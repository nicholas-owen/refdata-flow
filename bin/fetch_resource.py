#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Fetch one third-party resource, as described by its file in resources/.

Three things happen, in this order, and the order is the point:

  1. The version is worked out BEFORE downloading: the registry's fixed label, the
     date the publisher last changed the file (an HTTP HEAD request), the release
     the publisher states in a small file beside it, or - for a resource published
     per annotation release, such as a VEP cache - the release this genome's
     annotation is at. If the vault already holds that version, nothing is
     fetched. ClinVar is 200 MB and replaced weekly; this is what stops a daily
     run from fetching it daily.
  2. The file is downloaded, and hashed as it arrives. Before the first byte is
     written the size the publisher announces is compared with the space free
     where it is going, and afterwards with what actually arrived.
  3. If the publisher posts a checksum, in a file or with the download itself, the
     file is checked against it, and a mismatch fails the task - which the workflow
     retries, since a truncated transfer is the usual cause.

What comes out is the file exactly as published and a record of the fetch. It is
not yet known to belong to the genome it was requested for; bin/prepare_resource.py
decides that, against the genome's own sequence names and lengths. For a resource
whose sequence names have to be rewritten first (dbSNP), the table that says how -
an NCBI assembly report - is fetched here too, so that the step which uses it
needs no network.

Standard library only. This runs under the 'network' label, with whatever Python
the profile supplies. bin/fetch_database.py fetches a genome-less database with
the same code: fetch() below is the whole of it.

The output directory always exists afterwards. It holds either the file and
<name>.fetch.json, or a single file called SKIP saying why nothing was fetched.

Exit 3 means the request cannot be met and trying again would not help (no room
for the file, no annotation release to follow); the workflow does not retry it.
"""
import argparse
import base64
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime

import resource_registry

TIMEOUT = 120
TAG_SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]")
USER_AGENT = "refdata-flow"
EXIT_NO_RETRY = 3


class Refusal(Exception):
    """The request cannot be met as it stands; fetching again would change nothing."""


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


def as_text(body):
    """A small fetched file as text, whether or not it was gzipped (Pfam's is)."""
    if body[:2] == b"\x1f\x8b":
        body = gzip.decompress(body)
    return body.decode("utf-8", "replace")


def published_version(source, release=None):
    """
    The version the publisher is offering now, as (tag, how), or (None, how).

    None means it could not be determined - the request failed, or the answer had
    no date or release in it. The caller decides what that means; it is not an
    error here.
    """
    declared = source["version"]
    if declared == "annotation_release":
        return safe_tag(release or ""), "the release of this genome's annotation"
    if isinstance(declared, dict):
        how = f"the release stated in {declared['url']}"
        try:
            with open_url(declared["url"]) as response:
                found = re.search(declared["pattern"], as_text(response.read(1 << 16)))
            return (safe_tag(found.group(1)) if found else None), how
        except Exception as e:
            print(f"       Could not read {declared['url']} for the release: {e}",
                  file=sys.stderr, flush=True)
            return None, how
    if str(declared) != "last_modified":
        return safe_tag(declared), "fixed in the registry"
    try:
        with open_url(source["url"], method="HEAD") as response:
            return modified_tag(response.headers), "the publisher's Last-Modified date"
    except Exception as e:
        print(f"       Could not ask {source['url']} for its date: {e}", file=sys.stderr,
              flush=True)
        return None, "the publisher's Last-Modified date"


def announced_size(headers):
    value = str(headers.get("Content-Length") or "").strip()
    return int(value) if value.isdigit() else None


def download(url, destination):
    """Stream a URL to a file. Returns (sha256, md5, bytes, response headers)."""
    sha256, md5, size = hashlib.sha256(), hashlib.md5(), 0
    with open_url(url) as response:
        headers = response.headers
        expected = announced_size(headers)
        if expected is not None:
            # Asked before writing anything: a disk that fills part-way through a
            # 30 GB file takes everything else on it down with it.
            free = shutil.disk_usage(os.path.dirname(os.path.abspath(destination))).free
            if expected > free:
                raise Refusal(
                    f"{url} is {expected:,} bytes and only {free:,} are free in the work "
                    "directory. Nothing was downloaded. Note that a resource is held "
                    "there once as fetched, possibly once more as prepared, and then "
                    "copied into the vault")
        with open(destination, "wb") as out:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                out.write(chunk)
                sha256.update(chunk)
                md5.update(chunk)
                size += len(chunk)
    if expected is not None and size != expected:
        os.remove(destination)
        raise OSError(f"{url} announced {expected:,} bytes and {size:,} arrived; the "
                      "transfer was cut short and the file has been discarded")
    return sha256.hexdigest(), md5.hexdigest(), size, headers


def checksum_in_file(checksum, filename):
    """
    The checksum the publisher posts in a file: the hex string on the line that
    names our file, or the first one in it when no line does (a .md5 beside a
    single file often names a path from the publisher's own machine).
    """
    with open_url(checksum["url"]) as response:
        text = as_text(response.read(1 << 16))
    pattern = r"\b[0-9a-fA-F]{32,64}\b"
    # Named exactly: a list that has Pfam-A.hmm.gz also has Pfam-A.hmm.dat.gz.
    named = [line for line in text.splitlines()
             if re.search(pattern, line)
             and any(os.path.basename(word.lstrip("*")) == filename for word in line.split())]
    found = re.search(pattern, named[0] if named else text)
    if not found:
        raise ValueError(f"no checksum found in {checksum['url']}")
    return found.group(0).lower()


def checksum_in_header(checksum, headers):
    """
    A checksum sent with the download, such as Google Cloud Storage's
    `x-goog-hash: md5=<base64>`. None when the publisher sent none of that type:
    Cloud Storage omits md5 for an object that was uploaded in parts.
    """
    name = checksum["header"]
    values = headers.get_all(name) if hasattr(headers, "get_all") else [headers.get(name)]
    for part in ",".join(value for value in (values or []) if value).split(","):
        kind, separator, value = part.strip().partition("=")
        if separator and kind.strip().lower() == checksum["type"]:
            try:
                return base64.b64decode(value.strip()).hex()
            except ValueError:
                return None
    return None


def skip(outdir, reason):
    with open(os.path.join(outdir, "SKIP"), "w") as fh:
        fh.write(reason + "\n")
    print(f"  SKIP   {reason}", flush=True)


def fetch(source, name, label, held, outdir, no_release_check, release=None, what="Resource"):
    """
    Fetch what `source` describes into outdir, unless a version in `held` is
    already current. `source` is a resource's variant or a database's definition.

    Returns the part of the record that describes the download, or None when
    nothing was fetched (a SKIP file then says why). Raises Refusal when the
    request cannot be met, and exits 1 on a failed or corrupt transfer.
    """
    os.makedirs(outdir, exist_ok=True)
    print(f"{what} {label}", flush=True)

    # --- 1. is the published version already held? -------------------------------
    if held and no_release_check:
        return skip(outdir, f"{label}: held ({', '.join(sorted(held))}); release check is off")
    tag, how = published_version(source, release)
    if tag and tag in held:
        return skip(outdir, f"{label}: version {tag} is in the vault and is current")
    if tag is None and held:
        # Unreachable, or undated. Having some version, do not fetch again on a guess.
        return skip(outdir, f"{label}: held ({', '.join(sorted(held))}); the "
                            "publisher's current version could not be checked")

    # --- 2. download -----------------------------------------------------------------
    url = str(source["url"]).replace("{release}", tag or "")
    filename = os.path.basename(url)
    destination = os.path.join(outdir, filename)
    print(f"  FETCH  {url}", flush=True)
    if source.get("size"):
        print(f"         the registry says it is {source['size']}", flush=True)
    try:
        sha256, md5, size, headers = download(url, destination)
    except OSError as e:
        print(f"ERROR: could not fetch {url}: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"         {size:,} bytes", flush=True)
    if size == 0:
        os.remove(destination)
        print(f"ERROR: {url} returned an empty file", file=sys.stderr)
        sys.exit(1)
    if tag is None:
        # Asking beforehand told us nothing; the download's own headers may.
        tag = (modified_tag(headers) if source["version"] == "last_modified" else None) \
            or "default"

    # --- 3. the publisher's checksum, where there is one ---------------------------
    verified = "no checksum is published for this file"
    checksum = source.get("checksum")
    if checksum:
        actual = {"md5": md5, "sha256": sha256}.get(checksum.get("type"))
        if actual is None:
            print(f"ERROR: unsupported checksum type '{checksum.get('type')}'", file=sys.stderr)
            sys.exit(1)
        if checksum.get("header"):
            expected = checksum_in_header(checksum, headers)
            where = f"the publisher's {checksum['header']} header"
        else:
            expected = checksum_in_file(checksum, filename)
            where = checksum["url"]
        if expected is None:
            verified = f"the publisher sent no {checksum['type']} with this file"
        elif actual != expected:
            os.remove(destination)
            print(f"ERROR: {filename} does not match the publisher's {checksum['type']} "
                  f"({actual} here, {expected} at {where}). The transfer was "
                  "probably cut short; the file has been discarded.", file=sys.stderr)
            sys.exit(1)
        else:
            verified = f"{checksum['type']} matches {where}"
        print(f"         {verified}", flush=True)

    record = {
        "url": url,
        "version": tag,
        "version_from": how,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "file": filename,
        "file_bytes": size,
        "file_sha256": sha256,
        "file_md5": md5,
        "checksum_verified": verified,
    }

    # --- 4. the table its sequence names are rewritten through, if it needs one ---
    rename = source.get("rename")
    if rename:
        report = os.path.join(outdir, f"{name}.names.txt")
        print(f"  FETCH  {rename['report']}", flush=True)
        try:
            report_sha256, _, report_size, _ = download(rename["report"], report)
        except OSError as e:
            print(f"ERROR: could not fetch {rename['report']}: {e}", file=sys.stderr)
            sys.exit(1)
        if report_size == 0:
            print(f"ERROR: {rename['report']} returned an empty file", file=sys.stderr)
            sys.exit(1)
        record["names_report"] = {"url": rename["report"], "file": os.path.basename(report),
                                  "file_bytes": report_size, "file_sha256": report_sha256}

    print(f"         version {tag} ({how})", flush=True)
    return record


def write_record(outdir, name, record):
    with open(os.path.join(outdir, f"{name}.fetch.json"), "w") as fh:
        json.dump(record, fh, indent=4)


def main():
    parser = argparse.ArgumentParser(description="Fetch one registered resource.")
    parser.add_argument("--registry", required=True, help="The resource's file from resources/")
    parser.add_argument("--variant", required=True, help="Which build's file to fetch")
    parser.add_argument("--genome", required=True, help="The vault genome it is for (for messages)")
    parser.add_argument("--held", default="",
                        help="Comma-separated versions of this resource the genome already has")
    parser.add_argument("--annotation", default="",
                        help="The genome's annotation as asset:release, e.g. ensembl_gtf:116. "
                             "Needed only by a resource that follows the annotation's release.")
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

    held = {tag for tag in args.held.split(",") if tag}
    label = f"{name}:{args.variant} for {args.genome}"

    try:
        release = None
        if variant["version"] == "annotation_release":
            wanted = variant.get("annotation", "ensembl_gtf")
            asset, _, release = args.annotation.partition(":")
            if asset != wanted or not release:
                raise Refusal(
                    f"{name}:{args.variant} is published per release of a genome's "
                    f"'{wanted}' annotation, and {args.genome} was requested "
                    + (f"with '{args.annotation}'" if args.annotation else "without an annotation")
                    + f". Ask for it on a row whose annotation is '{wanted}'")
        download_record = fetch(variant, name, label, held, args.outdir,
                                args.no_release_check, release)
    except Refusal as e:
        print(f"ERROR: {e}.", file=sys.stderr)
        sys.exit(EXIT_NO_RETRY)
    if download_record is None:
        return

    record = {
        "resource": name,
        "variant": args.variant,
        "description": entry.get("description", ""),
        "terms": entry.get("terms", ""),
    }
    record.update(download_record)
    write_record(args.outdir, name, record)


if __name__ == "__main__":
    main()
