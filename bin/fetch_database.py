#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Fetch one genome-less database, as described by its file in databases/.

Swiss-Prot, Pfam, SILVA, a Kraken2 index: reference data that belongs to no
assembly, and so has no genome to be checked against and no place in refgenie.
What can be done for it is what bin/fetch_resource.py does for a resource before
its gate, and this runs the same code: work out the published version first and
fetch only if the tree does not hold it, compare the announced size with the
space free, hash the file as it arrives, and verify it against the publisher's
checksum where there is one.

Standard library only, like fetch_resource.py, and for the same reason.

The output directory always exists afterwards. It holds either the file and
<name>.fetch.json, or a single file called SKIP saying why nothing was fetched.
bin/file_database.py puts the result in <outdir>/databases/<name>/<version>/.
"""
import argparse
import sys

import fetch_resource
import resource_registry


def main():
    parser = argparse.ArgumentParser(description="Fetch one registered database.")
    parser.add_argument("--definition", required=True, help="The database's file from databases/")
    parser.add_argument("--held", default="",
                        help="Comma-separated versions of this database already in the tree")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--no-release-check", action="store_true",
                        help="Do not ask the publisher for its current version: any "
                             "version already held satisfies the request")
    args = parser.parse_args()

    try:
        entry = resource_registry.load_database(args.definition)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    name = entry["name"]
    held = {tag for tag in args.held.split(",") if tag}

    try:
        download_record = fetch_resource.fetch(entry, name, name, held, args.outdir,
                                               args.no_release_check, what="Database")
    except fetch_resource.Refusal as e:
        print(f"ERROR: {e}.", file=sys.stderr)
        sys.exit(fetch_resource.EXIT_NO_RETRY)
    if download_record is None:
        return

    record = {
        "database": name,
        "description": entry.get("description", ""),
        "terms": entry.get("terms", ""),
        "unpack": entry.get("unpack"),
    }
    record.update(download_record)
    fetch_resource.write_record(args.outdir, name, record)


if __name__ == "__main__":
    main()
