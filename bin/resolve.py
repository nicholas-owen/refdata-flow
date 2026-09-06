#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import csv
import os
import sys
import tempfile
from collections import namedtuple

import genomepy
from genomepy.providers import list_providers

# The resolved CSV carries `query` so the cache can tell *which request* produced
# a given assembly. Without it the cache could only be keyed on species, which
# meant an edited query (e.g. GRCh38 -> GRCh37) silently returned the previously
# resolved assembly. Downstream consumers (main.nf, update_refgenie.py) read by
# column name, so the extra column is additive and harmless to them.
RESOLVED_FIELDNAMES = ["species", "assembly", "provider", "annotation", "query", "release"]

# Which providers can actually be pinned to a release.
#
# This is deliberately a hard gate rather than a best-effort pass-through.
# genomepy's install_genome() takes **kwargs and forwards them to the provider,
# but only Ensembl reads a `version` from them - verified against 0.16.4, where
# `kwargs.get("version")` occurs exactly twice, both in providers/ensembl.py.
# GENCODE, UCSC and NCBI accept the argument and silently ignore it.
#
# Accepting `release` for those providers would therefore produce a request that
# looks pinned, reports as pinned, and is not - which is the failure class the
# annotation-recipe work in 0.9.1 was spent removing. Refuse instead.
#
# UCSC cannot be pinned even in principle: it patches assemblies in place and
# builds annotations live from MySQL, so there is no release to name. Its
# annotation is still *recorded*, by table modification date, at ingest time.
PINNABLE_PROVIDERS = {"ensembl"}

# How many candidates to show before truncating the interactive list. A query like
# "mouse" with no provider filter returns several hundred matches, nearly all of
# them gut metagenomes and viruses; printing them all is not a usable prompt.
DISPLAY_LIMIT = 12

# Exit codes. 2 means "some rows could not be resolved" as distinct from 1 for a
# fatal error; run_pipeline.sh treats them differently.
EXIT_UNRESOLVED = 2


# ---------------------------------------------------------------------------
# genomepy search adapter
# ---------------------------------------------------------------------------

SearchResult = namedtuple(
    "SearchResult",
    "assembly provider accession taxid has_annotation species other",
)

# genomepy/providers/__init__.py builds each row as:
#     list(row[:1]) + [provider_name] + list(row[1:])
# where `row` comes from that provider's own _genome_info_tuple(), giving:
#     [0] assembly  [1] provider  [2] accession  [3] taxid
#     [4] annotations  [5] species  [6] other
#
# Every provider returns those positions, but each implements _genome_info_tuple
# independently, so this is a convention rather than an enforced interface. Two
# things about it are worth pinning down:
#
#   * The TYPES are already inconsistent. `annotations` is a bool for Ensembl,
#     GENCODE and NCBI, but a list[bool] for UCSC (one flag per annotation source):
#       micMur2 (UCSC) | GCA_000165445.2 | 30608 | [True, False, True, False]
#     so treating index 4 as a boolean is wrong for UCSC today, not merely fragile.
#
#   * Passing size=True inserts `length` at index 5, shifting species to 6 and
#     other to 7. We always pass size=False explicitly - not just to keep the
#     positions stable, but because for NCBI and UCSC the size is fetched with one
#     network download per candidate, which on a 548-match query is unusable.
#     Genome size is captured after download instead, from the .fa.sizes file
#     genomepy leaves next to the FASTA.
_SEARCH_MIN_FIELDS = 6
_shape_verified = False


def _coerce_has_annotation(value):
    """Normalise the provider-dependent `annotations` field to a single bool."""
    if isinstance(value, (list, tuple, set)):
        return any(bool(v) for v in value)
    return bool(value)


def _verify_row_shape(row):
    """
    Sanity-check the first search row against the layout above, once per run.

    genomepy is pinned exactly in requirements.txt, so the shape is fixed today;
    this exists so that an upgrade which reorders the tuple fails loudly here
    instead of yielding confidently wrong assemblies downstream.
    """
    global _shape_verified
    if _shape_verified:
        return

    problems = []
    if len(row) < _SEARCH_MIN_FIELDS:
        problems.append(f"expected at least {_SEARCH_MIN_FIELDS} fields, got {len(row)}")
    else:
        if str(row[1]) not in list_providers():
            problems.append(
                f"field 1 was {row[1]!r}, which is not a known provider "
                f"({', '.join(list_providers())})"
            )
        try:
            int(row[3])
        except (TypeError, ValueError):
            problems.append(f"field 3 (taxid) was {row[3]!r}, which is not an integer")

    if problems:
        version = getattr(genomepy, "__version__", "unknown")
        raise RuntimeError(
            "genomepy.search() returned an unexpected row layout, so results cannot "
            "be interpreted safely.\n"
            f"  genomepy version: {version}\n"
            f"  observed row:     {row!r}\n"
            "  problems:         " + "; ".join(problems) + "\n"
            "This usually means genomepy changed its search result format. Check "
            "genomepy/providers/__init__.py and update the adapter in bin/resolve.py."
        )
    _shape_verified = True


def provider_reachable(provider_name):
    """
    Is a named provider answering right now? True / False / None if unknown.

    Needed because a provider being *offline* is otherwise indistinguishable from a
    query genuinely having *no match*. genomepy's online_providers() catches the
    ConnectionError for an unreachable provider, logs it as a warning, and carries
    on with the others; it only raises if every provider is down. Since this script
    searches all providers and filters afterwards, a single provider outage simply
    contributes no rows, and the filter then yields nothing.

    Observed for real: Ensembl was briefly unreachable, and the run reported
    "no matches found for 'GRCh38' with provider 'ensembl'" - sending the user to
    hunt for a typo in a query that was correct. The advice for the two cases is
    opposite (retry later vs fix the CSV), so they must not share a message.
    """
    try:
        from genomepy.providers import PROVIDERS
        provider = PROVIDERS.get(str(provider_name).strip().lower())
        if provider is None:
            return None
        return bool(provider.ping())
    except Exception:
        # Never let a diagnostic check fail the run; fall back to "unknown".
        return None


def unpinnable_reason(provider_name, release):
    """
    Why this provider cannot honour `release`, or None if it can.

    Cheap and offline, so it runs before the search rather than after it: a request
    that can never be satisfied should not cost a network round trip first.
    """
    name = str(provider_name).strip().lower()
    if not release or name in PINNABLE_PROVIDERS:
        return None
    return (
        f"release '{release}' was requested but provider '{name}' cannot be pinned "
        f"to a release (only {'/'.join(sorted(PINNABLE_PROVIDERS))} can). genomepy "
        "would accept the value and silently ignore it, so the request is refused "
        "rather than run unpinned. Remove the release, or use a provider that "
        "supports pinning. The annotation version is recorded either way at ingest"
    )


def validate_release(provider_name, release, assembly):
    """
    Check a release exists and actually contains this assembly. Reason, or None.

    Deferred to genomepy rather than reimplemented: EnsemblProvider.get_version()
    already checks both conditions and raises listing the valid alternatives, which
    is a better message than anything reconstructed here. The roadmap originally
    proposed rebuilding this with releases_with_assembly(); that would duplicate it.

    Called once the assembly is known - the second check is assembly-specific, since
    a release can exist without containing the assembly asked for.
    """
    if not release:
        return None

    name = str(provider_name).strip().lower()
    if name not in PINNABLE_PROVIDERS:
        return unpinnable_reason(name, release)

    try:
        from genomepy.providers import PROVIDERS
        provider = PROVIDERS.get(name)
        if provider is None:
            return None
        # PROVIDERS holds classes, not instances, and get_version needs the
        # provider's genome list to answer the assembly-specific half.
        provider().get_version(assembly, release)
        return None
    except Exception as e:
        return (
            f"release '{release}' is not usable for {assembly} on '{name}': {e}"
        )


def search_assemblies(query):
    """Run a genomepy search and return typed, shape-checked results."""
    rows = list(genomepy.search(query, size=False))
    if rows:
        _verify_row_shape(rows[0])
    results = []
    for r in rows:
        results.append(SearchResult(
            assembly=r[0],
            provider=str(r[1]),
            accession=r[2],
            taxid=r[3],
            has_annotation=_coerce_has_annotation(r[4]),
            species=r[5] if len(r) > 5 else None,
            other=r[6] if len(r) > 6 else None,
        ))
    return results


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------

def _normalise_species(value):
    """'Mus_musculus' and 'Mus musculus' both -> 'mus musculus'."""
    return str(value or "").replace("_", " ").strip().casefold()


def rank_results(results, species):
    """
    Order candidates by how well they match the requested species.

    A query like "mouse" returns 548 matches: GRCm39 and GRCm38, ~20 Ensembl strain
    assemblies, a mouse lemur and Mus pahari, then several hundred gut metagenomes,
    mouse-associated viruses and Candidatus Arthromitus. The species column in the
    request CSV already names what was wanted, but until now was used only as a
    directory name and never for matching.

    Ranking rather than filtering, deliberately: an exact-match filter would also
    drop the legitimate subspecies assemblies (CAST_EiJ_v3 is Mus musculus
    castaneus, PWK_PhJ_v3 is Mus musculus musculus), and silently hiding candidates
    is the failure mode this project keeps fixing. Nothing is removed; the useful
    entries simply come first.

    Tiers: exact species match, then same genus, then everything else. Sorting is
    stable, so provider order within a tier is whatever genomepy returned.
    """
    wanted = _normalise_species(species)
    if not wanted:
        return results, False

    wanted_genus = wanted.split(" ")[0]

    def tier(r):
        got = _normalise_species(r.species)
        if got and got == wanted:
            return 0
        if got and wanted_genus and got.split(" ")[0] == wanted_genus:
            return 1
        return 2

    ranked = sorted(results, key=tier)
    reordered = any(a is not b for a, b in zip(ranked, results))
    return ranked, reordered


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_resolved(output_file, rows):
    """
    Write the resolved CSV atomically.

    Previously this opened the output with mode="w" - over the same file the cache
    was read from - so an interruption mid-write lost every interactive choice made
    in the session and could leave a truncated cache behind. The temp file is
    created in the same directory because os.replace() is only atomic within a
    filesystem.

    Called after every resolution, not just at the end. Rewriting the whole file
    each time is O(n^2) writes, but n is a few dozen rows at most and the
    alternative (an append-only journal) is materially more complex for no benefit
    at this scale.
    """
    directory = os.path.dirname(os.path.abspath(output_file)) or "."
    # The output now lives in requests/, which normally exists because the request
    # CSV is read from it. Create it anyway: mkstemp fails with a bare FileNotFound
    # naming a hidden temp file, which says nothing about the real cause.
    os.makedirs(directory, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".resolve-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=RESOLVED_FIELDNAMES)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, output_file)
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def describe(result):
    """One-line summary of a candidate for the interactive list."""
    bits = [str(result.accession or "").strip(), str(result.taxid or "").strip()]
    bits.append("annotated" if result.has_annotation else "no annotation")
    if result.species:
        bits.append(str(result.species))
    return " | ".join(b for b in bits if b)


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Interactive selection
# ---------------------------------------------------------------------------

def choose_interactively(results, query, reordered, species):
    """
    Print the candidate list and ask the user to pick one.

    Only the first DISPLAY_LIMIT candidates are shown; 'a' lists the rest. Returns
    the chosen SearchResult.
    """
    total = len(results)
    shown = min(total, DISPLAY_LIMIT)
    print(f"Found {total} matches for '{query}':")
    if reordered:
        print(f"  (ordered by closeness to species '{species}')")

    def render(upto):
        for idx, r in enumerate(results[:upto], start=1):
            print(f"  [{idx}] {r.assembly} (Provider: {r.provider}) | {describe(r)}")

    render(shown)
    if total > shown:
        print(f"  ... showing {shown} of {total}. Enter 'a' to list them all.")

    choice = 0
    while choice < 1 or choice > total:
        # EOFError is left to propagate: the caller turns it into an unresolved row.
        # Catching it here would loop forever, since stdin stays closed.
        raw = input(f"Select an option [1-{total}], or 'a' for the full list: ").strip()
        if raw.lower() == "a":
            render(total)
            continue
        try:
            choice = int(raw)
        except ValueError:
            pass
    return results[choice - 1]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Resolve generic species queries into exact assemblies."
    )
    parser.add_argument("input_csv",
                        help="Request CSV (species, query, provider, annotation, aliases, release)")
    parser.add_argument(
        "output_csv", nargs="?", default="requests/requests_resolved.csv",
        help="Resolved CSV to write (default: requests/requests_resolved.csv)",
    )
    parser.add_argument(
        "--non-interactive", action="store_true",
        help="Never prompt. A query matching more than one assembly is reported as "
             "unresolved instead of asking, so unattended runs fail rather than hang.",
    )
    args = parser.parse_args()

    input_file = args.input_csv
    output_file = args.output_csv

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
                # Optional. Absent from older request CSVs, where DictReader yields
                # None - which must read as "unpinned", not crash.
                release = (row.get("release") or "").strip()

                if not query:
                    unresolved.append(
                        (lineno, species or "(no species)", "(empty)", "no query given")
                    )
                    continue

                # Refuse an unpinnable combination before searching. Only possible
                # here when the provider is stated in the request; when it is not,
                # the same check runs again once the selection determines it.
                if provider_filter:
                    reason = unpinnable_reason(provider_filter, release)
                    if reason:
                        print(f"ERROR: {reason}.", file=sys.stderr)
                        unresolved.append((lineno, species, query, reason))
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
                    cached_release = (cached_row.get("release") or "").strip()

                    if provider_filter and cached_provider.lower() != provider_filter:
                        print(
                            f"\n--- Cached resolution for '{species}' ({cached_assembly} from "
                            f"'{cached_provider}') does not match the requested provider "
                            f"'{provider_filter}'. Re-resolving. ---"
                        )
                    elif cached_release != release:
                        # Same reasoning as the provider check above: a release added,
                        # changed or removed since the cached run describes a different
                        # download, so the cached answer must not stand in for it.
                        print(
                            f"\n--- Cached resolution for '{species}' ({cached_assembly}, "
                            f"release '{cached_release or 'unpinned'}') does not match the "
                            f"requested release '{release or 'unpinned'}'. Re-resolving. ---"
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
                            "release": release,
                        })
                        write_resolved(output_file, resolved_rows)
                        continue

                # ------------------------------------------------------------------
                # Network search
                # ------------------------------------------------------------------
                print(f"\n--- Resolving query: '{query}' for {species} ---")

                results = search_assemblies(query)

                # Filter by provider if the user specified one
                if provider_filter:
                    results = [r for r in results if r.provider.lower() == provider_filter]

                if not results:
                    # Distinguish "the provider is down" from "the query matched
                    # nothing" before reporting - the two need opposite responses.
                    if provider_filter and provider_reachable(provider_filter) is False:
                        reason = (
                            f"provider '{provider_filter}' is unreachable, so '{query}' "
                            "could not be resolved. This is usually transient; re-run "
                            "when the provider is back. The request itself may be fine"
                        )
                    else:
                        reason = f"no matches found for '{query}'"
                        if provider_filter:
                            reason += f" with provider '{provider_filter}'"
                    # Not reason.capitalize(): that would also lower-case the rest of
                    # the string, mangling assembly names like 'GRCh38' into 'grch38'.
                    print(f"ERROR: {reason}. Skipping.", file=sys.stderr)
                    unresolved.append((lineno, species, query, reason))
                    continue

                results, reordered = rank_results(results, species)

                if len(results) == 1:
                    chosen = results[0]
                    print(f"Found exact match: {chosen.assembly} from {chosen.provider}.")
                elif args.non_interactive:
                    # Never auto-pick. Ranking makes the interactive list usable; it does
                    # not make automatic selection safe, and an unattended run quietly
                    # choosing a genome is the bug class fixed in 0.9.1. A hands-off run
                    # should already name its assembly.
                    top = ", ".join(f"{r.assembly} ({r.provider})" for r in results[:3])
                    reason = (
                        f"{len(results)} assemblies match '{query}' and --non-interactive "
                        f"was given, so no selection was made. Best candidates: {top}"
                    )
                    print(f"ERROR: {reason}.", file=sys.stderr)
                    unresolved.append((lineno, species, query, reason))
                    continue
                else:
                    try:
                        chosen = choose_interactively(results, query, reordered, species)
                    except EOFError:
                        # stdin is closed or not a terminal, so the prompt can never be
                        # answered. Treat it exactly like the --non-interactive case:
                        # report this row and carry on with the rest, rather than
                        # aborting the whole run on what is really a usage problem.
                        reason = (
                            f"{len(results)} assemblies match '{query}' and stdin is not "
                            "available to choose between them. Re-run attached to a "
                            "terminal, or pass --non-interactive and make the request "
                            "specific enough to resolve on its own"
                        )
                        print(f"ERROR: {reason}.", file=sys.stderr)
                        unresolved.append((lineno, species, query, reason))
                        continue
                    print(f"Selected: {chosen.assembly} from {chosen.provider}.")

                # Validate the release against what was actually selected. This is the
                # only point where both are known: the release must exist AND contain
                # this assembly, and when the request named no provider, whether the
                # provider can be pinned at all is only settled by the selection.
                if release:
                    reason = validate_release(chosen.provider, release, chosen.assembly)
                    if reason:
                        print(f"ERROR: {reason}.", file=sys.stderr)
                        unresolved.append((lineno, species, query, reason))
                        continue

                resolved_rows.append({
                    "species": species,
                    "assembly": chosen.assembly,
                    "provider": chosen.provider,
                    "annotation": annotation,
                    "query": query,
                    "release": release,
                })
                # Flush after every resolution so an interrupted interactive session
                # keeps the answers already given.
                write_resolved(output_file, resolved_rows)

    except Exception as e:
        print(f"Error resolving requests: {e}", file=sys.stderr)
        sys.exit(1)

    # Final write, so the file exists even when nothing resolved.
    write_resolved(output_file, resolved_rows)

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
        sys.exit(EXIT_UNRESOLVED)


if __name__ == "__main__":
    main()
