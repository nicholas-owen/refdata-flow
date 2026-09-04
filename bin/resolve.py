#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import sys
import os
import csv
import genomepy

# The resolved CSV carries `query` so the cache can tell *which request* produced
# a given assembly. Without it the cache could only be keyed on species, which
# meant an edited query (e.g. GRCh38 -> GRCh37) silently returned the previously
# resolved assembly. Downstream consumers (main.nf, update_refgenie.py) read by
# column name, so the extra column is additive and harmless to them.
RESOLVED_FIELDNAMES = ["species", "assembly", "provider", "annotation", "query"]


def load_cache(output_file):
    """
    Load previous resolutions so an unchanged request can skip the slow network
    search and the interactive prompt.

    Keyed on (species, query) rather than species alone: a species may legitimately
    be needed at more than one assembly, and an edited query must not be answered
    from a stale entry.

    A resolved CSV written before the `query` column existed yields an empty query
    for every row, so its entries simply fail to match and are re-resolved. That is
    the safe direction - a one-off re-prompt, never a silently wrong assembly.
    """
    cache = {}
    if not os.path.exists(output_file):
        return cache
    try:
        with open(output_file, mode="r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                species = (row.get("species") or "").strip()
                query = (row.get("query") or "").strip()
                if species:
                    cache[(species, query)] = row
    except Exception as e:
        # Previously this was `except Exception: pass`, which made a corrupt cache
        # indistinguishable from an absent one. Re-resolving everything is still the
        # correct recovery, but the user is now told why it is happening.
        print(
            f"Warning: could not read existing resolutions from '{output_file}' ({e}).\n"
            "         Every request will be re-resolved.",
            file=sys.stderr,
        )
        return {}
    return cache


def main():
    if len(sys.argv) < 2:
        print("Usage: python bin/resolve.py <requests_raw.csv> [output_resolved.csv]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "requests_resolved.csv"

    resolved_cache = load_cache(output_file)

    resolved_rows = []
    # Rows that could not be turned into a concrete assembly. Collected rather than
    # silently skipped so the run can report them and exit non-zero.
    unresolved = []

    try:
        with open(input_file, mode="r", newline="") as infile:
            reader = csv.DictReader(infile)
            # start=2 so the number reported matches the line in the CSV (row 1 is the header).
            for lineno, row in enumerate(reader, start=2):
                species = (row.get("species") or "").strip()
                query = (row.get("query") or "").strip()
                provider_filter = (row.get("provider") or "").strip().lower()
                annotation = (row.get("annotation") or "false").strip()

                if not query:
                    unresolved.append(
                        (lineno, species or "(no species)", "(empty)", "no query given")
                    )
                    continue

                # ------------------------------------------------------------------
                # Cache lookup
                # ------------------------------------------------------------------
                cached_row = resolved_cache.get((species, query))
                if cached_row is not None:
                    cached_assembly = (cached_row.get("assembly") or "").strip()
                    cached_provider = (cached_row.get("provider") or "").strip()

                    # A provider filter added or changed since the cached run must not be
                    # ignored: the cached assembly may come from a provider the user has
                    # now explicitly ruled out.
                    if provider_filter and cached_provider.lower() != provider_filter:
                        print(
                            f"\n--- Cached resolution for '{species}' ({cached_assembly} from "
                            f"'{cached_provider}') does not match the requested provider "
                            f"'{provider_filter}'. Re-resolving. ---"
                        )
                    else:
                        print(
                            f"\n--- '{query}' for {species} already resolved to {cached_assembly}. "
                            "Skipping search. ---"
                        )
                        # Rebuild the row rather than mutating the cached one, so stale or
                        # unexpected columns in an older resolved CSV are not carried through.
                        resolved_rows.append({
                            "species": species,
                            "assembly": cached_assembly,
                            "provider": cached_provider,
                            "annotation": annotation,
                            "query": query,
                        })
                        continue

                # ------------------------------------------------------------------
                # Network search
                # ------------------------------------------------------------------
                print(f"\n--- Resolving query: '{query}' for {species} ---")

                # genomepy search returns a generator of lists
                results = list(genomepy.search(query))

                # Filter by provider if the user specified one
                if provider_filter:
                    results = [r for r in results if r[1].lower() == provider_filter]

                if not results:
                    reason = f"no matches found for '{query}'"
                    if provider_filter:
                        reason += f" with provider '{provider_filter}'"
                    # Not reason.capitalize(): that would also lower-case the rest of
                    # the string, mangling assembly names like 'GRCh38' into 'grch38'.
                    print(f"ERROR: {reason}. Skipping.", file=sys.stderr)
                    unresolved.append((lineno, species, query, reason))
                    continue

                if len(results) == 1:
                    match = results[0]
                    assembly = match[0]
                    provider = match[1]
                    print(f"Found exact match: {assembly} from {provider}.")
                else:
                    # Multiple matches, interactive prompt
                    print(f"Found {len(results)} matches for '{query}':")
                    for idx, match in enumerate(results):
                        name = match[0]
                        provider = match[1]
                        extra = " | ".join(str(x) for x in match[2:6] if x)
                        print(f"  [{idx + 1}] {name} (Provider: {provider}) | {extra}")

                    choice = 0
                    while choice < 1 or choice > len(results):
                        try:
                            choice_str = input(f"Select an option [1-{len(results)}]: ")
                            choice = int(choice_str)
                        except ValueError:
                            pass

                    selected_match = results[choice - 1]
                    assembly = selected_match[0]
                    provider = selected_match[1]
                    print(f"Selected: {assembly} from {provider}.")

                resolved_rows.append({
                    "species": species,
                    "assembly": assembly,
                    "provider": provider,
                    "annotation": annotation,
                    "query": query,
                })

    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Write whatever did resolve, so a partial run is not thrown away.
    with open(output_file, mode="w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=RESOLVED_FIELDNAMES)
        writer.writeheader()
        for row in resolved_rows:
            writer.writerow(row)

    print(f"\nResolved {len(resolved_rows)} request(s). Written to {output_file}.")

    # ----------------------------------------------------------------------
    # Report unresolved rows and fail. Previously these were skipped silently and
    # the script still exited 0, so a user who asked for five genomes and received
    # three had no way to notice. Exit code 2 distinguishes "some rows could not be
    # resolved" from a fatal error (1), letting run_pipeline.sh explain the difference.
    # ----------------------------------------------------------------------
    if unresolved:
        print(
            f"\nERROR: {len(unresolved)} request(s) could not be resolved:",
            file=sys.stderr,
        )
        for lineno, species, query, reason in unresolved:
            # ASCII only: this goes to stderr, and a Windows console in a legacy
            # code page renders an em-dash as a replacement character.
            print(f"  - line {lineno}: {species} / '{query}' -> {reason}", file=sys.stderr)
        print(
            "\nFix or remove these rows in the request CSV and re-run. Nothing was downloaded.",
            file=sys.stderr,
        )
        sys.exit(2)


if __name__ == "__main__":
    main()
