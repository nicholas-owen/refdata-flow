#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Fetch one third-party resource, as described by its file in resources/.

Three things happen, in this order, and the order is the point:

  1. The version is worked out BEFORE downloading: the registry's fixed label, or
     the date the publisher last changed the file, from an HTTP HEAD request. If
     the vault already holds that version, nothing is fetched. ClinVar is 200 MB
     and replaced weekly; this is what stops a daily run from fetching it daily.
  2. The file is downloaded, and hashed as it arrives.
  3. If the publisher posts a checksum, the file is checked against it, and a
     mismatch fails the task - which the workflow retries, since a truncated
     transfer is the usual cause.

What comes out is the file exactly as published and a record of the fetch. It is
not yet known to belong to the genome it was requested for; bin/prepare_resource.py
decides that, against the genome's own sequence names and lengths.

Standard library only. This runs under the 'network' label, with whatever Python
the profile supplies.

The output directory always exists afterwards. It holds either the file and
<name>.fetch.json, or a single file called SKIP saying why nothing was fetched.
"""
import argparse
import hashlib
import json
import os
import re
import sys
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime

import resource_registry

TIMEOUT = 120
TAG_SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]")
USER_AGENT = "refdata-flow"


def open_url(url, method="GET"):
    """The one place the network is touched, so tests can stand in for it."""
    request = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=TIMEOUT)


def safe_tag(value):
    cleaned = TAG_SAFE_RE.sub("-", str(value)).strip("-")
    return cleaned or None


def modified_tag(headers):
    """A response's Last-Modified as YYYYMMDD, or None."""
    value = headers.get("Last-Modified")
    if not value:
        return None
    try:
        return parsedate_to_datetime(value).strftime("%Y%m%d")
    except (TypeError, ValueError):
        return None


def published_version(variant):
    """
    The version the publisher is offering now, as (tag, how), or (None, how).

    None means it could not be determined - the HEAD request failed, or the server
    sent no date. The caller decides what that means; it is not an error here.
    """
    declared = str(variant["version"])
    if declared != "last_modified":
        return safe_tag(declared), "fixed in the registry"
    try:
        with open_url(variant["url"], method="HEAD") as response:
            return modified_tag(response.headers), "the publisher's Last-Modified date"
    except Exception as e:
        print(f"       Could not ask {variant['url']} for its date: {e}", file=sys.stderr,
              flush=True)
        return None, "the publisher's Last-Modified date"


def download(url, destination):
    """Stream a URL to a file. Returns (sha256, md5, bytes, response headers)."""
    sha256, md5, size = hashlib.sha256(), hashlib.md5(), 0
    with open_url(url) as response, open(destination, "wb") as out:
        headers = response.headers
        for chunk in iter(lambda: response.read(1 << 20), b""):
            out.write(chunk)
            sha256.update(chunk)
            md5.update(chunk)
            size += len(chunk)
    return sha256.hexdigest(), md5.hexdigest(), size, headers


def published_checksum(checksum):
    """The checksum the publisher posts: the first hex string in the file."""
    with open_url(checksum["url"]) as response:
        text = response.read(4096).decode("utf-8", "replace")
    found = re.search(r"\b[0-9a-fA-F]{32,64}\b", text)
    if not found:
        raise ValueError(f"no checksum found in {checksum['url']}")
    return found.group(0).lower()


def skip(outdir, reason):
    with open(os.path.join(outdir, "SKIP"), "w") as fh:
        fh.write(reason + "\n")
    print(f"  SKIP   {reason}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Fetch one registered resource.")
    parser.add_argument("--registry", required=True, help="The resource's file from resources/")
    parser.add_argument("--variant", required=True, help="Which build's file to fetch")
    parser.add_argument("--genome", required=True, help="The vault genome it is for (for messages)")
    parser.add_argument("--held", default="",
                        help="Comma-separated versions of this resource the genome already has")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--no-release-check", action="store_true",
                        help="Do not ask the publisher for its current version: any "
                             "version already held satisfies the request")
    args = parser.parse_args()

    try:
        entry = resource_registry.load_entry(args.registry)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    name = entry["name"]
    variant = entry["variants"].get(args.variant)
    if variant is None:
        print(f"ERROR: resource '{name}' has no variant '{args.variant}'", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    held = {tag for tag in args.held.split(",") if tag}
    label = f"{name}:{args.variant} for {args.genome}"
    print(f"Resource {label}", flush=True)

    # --- 1. is the published version already in the vault? ----------------------
    if held and args.no_release_check:
        return skip(args.outdir, f"{label}: held ({', '.join(sorted(held))}); release check is off")
    tag, how = published_version(variant)
    if tag and tag in held:
        return skip(args.outdir, f"{label}: version {tag} is in the vault and is current")
    if tag is None and held:
        # Unreachable, or undated. Having some version, do not fetch again on a guess.
        return skip(args.outdir, f"{label}: held ({', '.join(sorted(held))}); the "
                                 "publisher's current version could not be checked")

    # --- 2. download ---------------------------------------------------------------
    filename = os.path.basename(variant["url"])
    destination = os.path.join(args.outdir, filename)
    print(f"  FETCH  {variant['url']}", flush=True)
    sha256, md5, size, headers = download(variant["url"], destination)
    print(f"         {size:,} bytes", flush=True)
    if size == 0:
        os.remove(destination)
        print(f"ERROR: {variant['url']} returned an empty file", file=sys.stderr)
        sys.exit(1)
    if tag is None:
        # The HEAD told us nothing; the download's own headers may.
        tag = (modified_tag(headers) if str(variant["version"]) == "last_modified" else None) \
            or "default"

    # --- 3. the publisher's checksum, where there is one -------------------------
    verified = "no checksum is published for this file"
    checksum = variant.get("checksum")
    if checksum:
        expected = published_checksum(checksum)
        actual = {"md5": md5, "sha256": sha256}.get(checksum.get("type"))
        if actual is None:
            print(f"ERROR: unsupported checksum type '{checksum.get('type')}'", file=sys.stderr)
            sys.exit(1)
        if actual != expected:
            os.remove(destination)
            print(f"ERROR: {filename} does not match the publisher's {checksum['type']} "
                  f"({actual} here, {expected} at {checksum['url']}). The transfer was "
                  "probably cut short; the file has been discarded.", file=sys.stderr)
            sys.exit(1)
        verified = f"{checksum['type']} matches {checksum['url']}"
        print(f"         {verified}", flush=True)

    record = {
        "resource": name,
        "variant": args.variant,
        "description": entry.get("description", ""),
        "url": variant["url"],
        "terms": entry.get("terms", ""),
        "version": tag,
        "version_from": how,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "file": filename,
        "file_bytes": size,
        "file_sha256": sha256,
        "file_md5": md5,
        "checksum_verified": verified,
    }
    with open(os.path.join(args.outdir, f"{name}.fetch.json"), "w") as fh:
        json.dump(record, fh, indent=4)
    print(f"         version {tag} ({how})", flush=True)


if __name__ == "__main__":
    main()
