#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Decide which resolved requests still need downloading, by looking in the vault.

Before this existed, the only thing that stopped a finished genome being
downloaded again was Nextflow's resume cache in work/. That is a cache, not a
record: `run_pipeline.sh --clean` deletes it, and on an HPC work/ is not kept at
all, so every run fetched every genome from scratch. The vault's config is kept,
lives beside the data, and already says what has been built. So it is asked first.

A request is satisfied when the vault holds:

  * the genome, under its vault name (<assembly>_<provider>), with the sequence
    dictionary actually on disk; and
  * if annotation was requested, the provider's annotation asset at the release
    the request means:
      - pinned (Ensembl): the pinned release;
      - unpinned: the release the provider publishes NOW, worked out from the
        same source genomepy would download and by the same rules that tagged
        the asset at ingest (bin/download.py, release_from_source).

So an unpinned request tracks the provider: when a newer release appears, the row
is fetched again and the vault gains a tag; the older one stays reachable.

Two things this never does:

  * Fail because a provider is unreachable. If the current release cannot be
    determined and the vault holds some release of the annotation, the row is
    treated as satisfied. An offline run must still work, and re-downloading a
    genome that is already there is the wrong response to a network blip.
  * Guess in the other direction. Anything else it is unsure of is fetched:
    doing work twice is recoverable, and silently skipping a genome is not.

Writes one line per request that needs fetching to --fetch-keys, in the key
format update_refgenie.source_key() defines, and a report to stdout.
"""
import argparse
import csv
import json
import os
import sys

import yaml

import download
import update_refgenie as vault


def annotation_source(provider, assembly):
    """
    What genomepy would record as the 'annotation url' for a download made now.

    The same string its README.txt gets after a real download, obtained without
    downloading: a URL for the file-based providers, and for UCSC the table it
    would read, in the "UCSC MySQL database: ..." form release_from_source() parses.
    """
    from genomepy.providers import create

    instance = create(str(provider).strip().lower())
    link = instance.get_annotation_download_link(assembly)
    if instance.name.lower() == "ucsc":
        return f"UCSC MySQL database: {assembly}, table: {link}"
    return link


def current_annotation_release(provider, assembly):
    """
    The tag the provider's current annotation would be ingested under, or None.

    None for any failure at all. This asks the network a question whose answer is
    only an optimisation, so nothing that goes wrong here may stop the run.
    """
    try:
        source = annotation_source(provider, assembly)
        release, _kind = download.release_from_source(source, provider, f"[{assembly}]")
        return download.safe_tag(release) if release else None
    except Exception as e:
        print(f"[{assembly}]       Could not determine the current annotation release "
              f"from '{provider}': {e}", file=sys.stderr, flush=True)
        return None


def load_config(config_file):
    """The vault's config, or None when there is no vault yet (or it is unreadable)."""
    try:
        with open(config_file, "r") as f:
            return yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError):
        return None


def assess(row, config, config_file, current_release=current_annotation_release,
           check_releases=True):
    """
    Does this request need downloading? Returns (fetch, reason).

    `current_release` is a parameter so tests can answer for the provider.
    """
    assembly = (row.get("assembly") or "").strip()
    provider = (row.get("provider") or "").strip()
    annotation = (row.get("annotation") or "false").strip().lower() == "true"
    pinned = (row.get("release") or "").strip()

    if config is None:
        return True, "there is no vault yet"

    try:
        vault_name = vault.vault_genome_name(assembly, provider)
    except ValueError as e:
        return True, str(e)

    digest = None
    for candidate, data in (config.get("genomes") or {}).items():
        if vault_name in ((data or {}).get("aliases") or []):
            digest = candidate
            break
    if digest is None:
        return True, f"'{vault_name}' is not in the vault"

    has_dictionary = vault.fasta_has_dictionary(config_file, vault_name)
    if has_dictionary is None:
        return True, f"'{vault_name}' has no fasta asset"
    if has_dictionary is False:
        return True, f"'{vault_name}' has no sequence dictionary; it is built from the download"

    if not annotation:
        return False, "genome is in the vault"

    asset = vault.ANNOTATION_RECIPES.get(provider.lower())
    if asset is None:
        # The ingest reports this properly; let it, rather than hide the row here.
        return True, f"no annotation recipe for provider '{provider}'"

    # On disk as well as in the config: refgenie's config can list a tag whose
    # files were never written, and that is not a reason to skip a download.
    tags = vault.asset_tags_on_disk(config, config_file, digest, asset)
    if not tags:
        return True, f"genome is in the vault but its {asset} annotation is not"

    if pinned:
        if pinned in tags:
            return False, f"genome and {asset}:{pinned} are in the vault"
        return True, (f"release {pinned} of {asset} is not in the vault "
                      f"(it has {', '.join(sorted(tags))})")

    held = ", ".join(sorted(tags))
    if provider.lower() == "url":
        # An address is one file. There is no provider to ask what is current, and
        # no newer release to follow: a different file is a different address.
        return False, f"genome and {asset} ({held}) are in the vault, from the address given"
    if not check_releases:
        return False, f"genome and {asset} ({held}) are in the vault; release check is off"

    current = current_release(provider, assembly)
    if current is None:
        return False, (f"genome and {asset} ({held}) are in the vault; the provider's "
                       "current release could not be checked, so nothing is refetched")
    if current in tags:
        return False, f"genome and {asset}:{current} are in the vault, and {current} is current"
    return True, f"{provider} now publishes {asset} release {current}; the vault has {held}"


def address_conflict(row, config, config_file):
    """
    For a genome requested by address: is its name already in the vault, taken
    from a DIFFERENT address? Returns a description, or None.

    A named assembly from a provider is the same thing whoever asks for it. A
    genome requested by address is whatever file sits at that address, under
    whatever name the request gives it. So the same name with another address is
    a different genome - and the vault would not notice: the name is there, the
    preflight would skip the download, and the request would be "satisfied" by
    the file somebody else fetched. Fetching it would not help either, because
    the ingest keeps the genome a name already has.

    So it is neither skipped nor fetched. It is refused, here, before anything
    is downloaded, and the request has to use another name.

    Only the vault's own records can say where a genome came from. One ingested
    before the address was recorded has nothing to compare, and is let through.
    """
    provider = (row.get("provider") or "").strip()
    assembly = (row.get("assembly") or "").strip()
    wanted = (row.get("url") or "").strip()
    if provider.lower() != "url" or not wanted or config is None:
        return None
    try:
        vault_name = vault.vault_genome_name(assembly, provider)
    except ValueError:
        return None
    for digest, data in (config.get("genomes") or {}).items():
        if vault_name not in ((data or {}).get("aliases") or []):
            continue
        recorded = set()
        for tag in vault.asset_tags_on_disk(config, config_file, digest, "provenance"):
            path = vault.asset_file(config, config_file, digest, "provenance", tag, "provenance")
            try:
                with open(path) as fh:
                    address = (json.load(fh).get("url_requested") or "").strip()
            except (OSError, ValueError):
                continue
            if address:
                recorded.add(address)
        if recorded and wanted not in recorded:
            return (f"'{vault_name}' is already in the vault, taken from "
                    f"{', '.join(sorted(recorded))}. This request gives {wanted}. A name "
                    "identifies one genome: give this one a different name in the query "
                    "column")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="List the resolved requests the vault does not already satisfy."
    )
    parser.add_argument("resolved_csv", help="The resolved request CSV from bin/resolve.py")
    parser.add_argument("outdir", help="Output directory holding the refgenie vault")
    parser.add_argument("--fetch-keys", required=True, metavar="FILE",
                        help="Write the key of each request to download here, one per line")
    parser.add_argument(
        "--no-release-check", action="store_true",
        help="Do not ask providers for their current release. An unpinned request "
             "is then satisfied by any release of its annotation already in the vault.",
    )
    args = parser.parse_args()

    config_file = os.path.join(args.outdir, "refgenie_export.yaml")
    config = load_config(config_file)

    try:
        with open(args.resolved_csv, "r", newline="") as f:
            rows = list(csv.DictReader(f))
    except OSError as e:
        print(f"Preflight: cannot read '{args.resolved_csv}': {e}", file=sys.stderr)
        sys.exit(1)

    keys, conflicts = [], []
    print(f"Preflight: checking {len(rows)} request(s) against the vault at {args.outdir}",
          flush=True)
    for row in rows:
        assembly = (row.get("assembly") or "").strip()
        provider = (row.get("provider") or "").strip()
        release = (row.get("release") or "").strip()
        conflict = address_conflict(row, config, config_file)
        if conflict:
            print(f"  STOP   {assembly} from {provider}: {conflict}", flush=True)
            conflicts.append(conflict)
            continue
        try:
            fetch, reason = assess(row, config, config_file,
                                   check_releases=not args.no_release_check)
        except Exception as e:
            # A bug here must cost a download, never a genome.
            fetch, reason = True, f"the check itself failed ({e})"
        label = f"{assembly} from {provider}" + (f", release {release}" if release else "")
        print(f"  {'FETCH' if fetch else 'SKIP '}  {label}: {reason}", flush=True)
        if fetch:
            keys.append(vault.source_key(provider, assembly, release))

    with open(args.fetch_keys, "w") as f:
        f.write("".join(key + "\n" for key in keys))

    if conflicts:
        # Before any download. The run cannot give the user what the request
        # describes, and carrying on would report success for a genome it did
        # not fetch.
        print(f"Preflight: stopping. {len(conflicts)} request(s) name a genome that is in "
              "the vault from a different address (see STOP above). Nothing was "
              "downloaded.", file=sys.stderr, flush=True)
        sys.exit(1)

    print(f"Preflight: {len(keys)} to download, {len(rows) - len(keys)} already in the vault.",
          flush=True)


if __name__ == "__main__":
    main()
