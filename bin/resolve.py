#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import csv
import os
import re
import sys
import tempfile
from collections import namedtuple

import genomepy
from genomepy.providers import PROVIDERS, list_providers
# genomepy logs through loguru's single global logger. Importing it here gives
# this script the same object, to restyle its output and to listen to it.
from loguru import logger as genomepy_logger

# The resolved CSV carries `query` so the cache can tell *which request* produced
# a given assembly. Without it the cache could only be keyed on species, which
# meant an edited query (e.g. GRCh38 -> GRCh37) silently returned the previously
# resolved assembly. Downstream consumers (main.nf, update_refgenie.py) read by
# column name, so the extra column is additive and harmless to them.
#
# `aliases` is carried through for the same reason, one step further on. The ingest
# used to join each resolved row back to the request CSV on (species, query) to find
# its aliases, which cannot tell apart two rows asking for the same assembly from two
# providers - Mus_musculus/GRCm39 from Ensembl and from GENCODE share that key, so
# both rows were handed the last row's aliases and the first genome ended up
# carrying the second's. A resolved row that carries its own aliases needs no join.
RESOLVED_FIELDNAMES = ["species", "assembly", "provider", "annotation", "query",
                       "release", "aliases"]

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


# ---------------------------------------------------------------------------
# genomepy's log output
# ---------------------------------------------------------------------------

# genomepy's own format, which colours the level name blue whatever the level:
#   <green>{time}</green> | <blue>{level}</blue> | {message}
# So "Ensembl appears to be offline" - the one line in a long search that changes
# what the user should do next - looked exactly like the progress lines around it.
# INFO keeps genomepy's styling; WARNING and above are red, level and message both.
_LOG_PREFIX = "<green>{time:HH:mm:ss}</green> <bold>|</bold> "
_LOG_INFO = _LOG_PREFIX + "<blue>{level}</blue> <bold>|</bold> {message}\n{exception}"
_LOG_WARNING = (_LOG_PREFIX + "<red><bold>{level}</bold></red> <bold>|</bold> "
                "<red>{message}</red>\n{exception}")

# Providers genomepy reported as offline during the search now running. See
# search_assemblies() for why this is recorded rather than checked afterwards.
_OFFLINE_RE = re.compile(r"^(.+?) appears to be offline")
_offline_during_search = set()


def _log_format(record):
    return _LOG_WARNING if record["level"].no >= 30 else _LOG_INFO


def _note_offline(message):
    """loguru sink: remember which providers genomepy says it could not reach."""
    found = _OFFLINE_RE.match(message.record["message"])
    if found:
        for name in found.group(1).split(","):
            _offline_during_search.add(name.strip().lower())


def configure_genomepy_logging(stream=None, colorize=None):
    """
    Replace genomepy's log handler with one that makes warnings stand out, and add
    a second that listens for providers going offline.

    `colorize=None` leaves the decision to loguru, which colours a terminal and
    writes plain text to a file or a pipe - so a log captured by a scheduler holds
    no escape codes.
    """
    genomepy_logger.remove()
    genomepy_logger.add(stream or sys.stderr, format=_log_format, level="INFO",
                        colorize=colorize)
    genomepy_logger.add(_note_offline, level="WARNING", format="{message}")


def search_assemblies(query, provider=None):
    """
    Run a genomepy search. Returns (typed, shape-checked results, offline providers).

    With `provider`, only that provider is searched. The resolver used to search
    all of them for every row and discard the rest, so a request naming Ensembl
    waited several minutes on NCBI's assembly summaries, and an outage at a
    provider nobody had asked for still appeared in the output.

    The second value is the set of providers genomepy could not reach WHILE this
    search ran, taken from its own warnings. Asking afterwards is not the same
    question. A search across all providers takes minutes - NCBI alone is several -
    and a provider that was down at the start can be back by the end. Seen for
    real: Ensembl was offline when the search reached it, so it contributed no
    rows; five minutes later a ping succeeded, and the run reported "no matches
    found for 'GRCh38' with provider 'ensembl'", sending the user to look for a
    typo in a request that was correct.
    """
    _offline_during_search.clear()
    try:
        rows = list(genomepy.search(query, provider=provider or None, size=False))
    except ConnectionError:
        # genomepy warns and carries on while at least one provider answers, and
        # raises only when none does - which, with one provider named, is whenever
        # that provider is down. Either way it is an outage to report against the
        # row, not a reason to abandon every other row in the request.
        rows = []
        _offline_during_search.update(
            [provider.lower()] if provider else [name.lower() for name in list_providers()]
        )
    offline = set(_offline_during_search)
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
    return results, offline


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

    Keyed on (species, query, release) rather than species alone: a species may
    legitimately be needed at more than one assembly, an edited query must not be
    answered from a stale entry, and two releases of one assembly are two requests.
    Release was previously compared after the lookup instead of being part of the
    key, so two rows differing only in release shared one slot and evicted each
    other on every run.

    Each key holds a LIST, because the key still does not identify a request on its
    own: the same query against two providers (GRCm39 from Ensembl and from GENCODE)
    is two rows with one key. cached_resolution() picks between them by provider.

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
                release = (row.get("release") or "").strip()
                if species:
                    cache.setdefault((species, query, release), []).append(row)
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


def cached_resolution(cache, species, query, release, provider_filter):
    """
    The one cached row that answers this request, or None to re-resolve.

    With a provider named in the request, only a cached row from that provider
    counts: a filter added or changed since the cached run must not be ignored, as
    the cached assembly may come from a provider the user has now ruled out.

    Anything other than exactly one candidate is a miss. Two candidates means the
    cache cannot say which was meant, and guessing is how a request ends up with
    another row's genome; re-resolving costs one search.
    """
    rows = cache.get((species, query, release), [])
    if provider_filter:
        rows = [r for r in rows
                if (r.get("provider") or "").strip().lower() == provider_filter]
    return rows[0] if len(rows) == 1 else None


# ---------------------------------------------------------------------------
# Request consistency
# ---------------------------------------------------------------------------

def wants_annotation(value):
    """The same test main.nf and update_refgenie.py apply to the annotation column."""
    return str(value or "").strip().lower() == "true"


def parse_aliases(value):
    """'yeast; sacCer' -> ['yeast', 'sacCer']."""
    return [a.strip() for a in str(value or "").split(";") if a.strip()]


class RequestLedger:
    """
    What the rows accepted so far have claimed, so a later row can be checked
    against them before anything is downloaded.

    Both checks are offline and both were previously left to fail late or not at all:

      * Two rows resolving to the same (provider, assembly, release) became two
        identical download tasks publishing into one directory. Wasteful when the
        rows agree; when they disagree on `annotation`, the provenance.json that
        survives can be the sequence-only task's, and the GTF is then ingested
        under the 'default' tag with its release lost.

      * An alias names exactly one genome. The ingest already refuses to repoint
        one, but only discovers the clash after every genome has been fetched.

    The same alias on two rows is fine when both rows are the same genome - two
    releases of one assembly, say - so the alias check compares (provider, assembly)
    and ignores release.
    """

    def __init__(self):
        self._targets = {}       # (provider, assembly, release) -> line number
        self._alias_owner = {}   # alias -> ((provider, assembly), line number)

    def conflict(self, lineno, provider, assembly, release, aliases):
        """Why this row cannot stand beside the accepted ones, or None. Claims on None."""
        genome = (str(provider).strip().lower(), assembly)
        target = genome + (release,)

        earlier = self._targets.get(target)
        if earlier is not None:
            pinned = f" release {release}" if release else ""
            return (
                f"resolves to {assembly} from '{provider}'{pinned}, which line "
                f"{earlier} already requests. Merge the two rows into one, putting "
                "both sets of aliases on it"
            )

        for alias in aliases:
            owner = self._alias_owner.get(alias)
            if owner is not None and owner[0] != genome:
                (other_provider, other_assembly), other_line = owner
                return (
                    f"alias '{alias}' is already given to {other_assembly} from "
                    f"'{other_provider}' on line {other_line}. An alias can name only "
                    "one genome; give this row a different one"
                )

        self._targets[target] = lineno
        for alias in aliases:
            self._alias_owner.setdefault(alias, (genome, lineno))
        return None


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

    configure_genomepy_logging()

    resolved_cache = load_cache(output_file)
    ledger = RequestLedger()

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
                aliases = parse_aliases(row.get("aliases"))

                if not query:
                    unresolved.append(
                        (lineno, species or "(no species)", "(empty)", "no query given")
                    )
                    continue

                # A provider genomepy has never heard of. This used to pass through
                # unnoticed - every provider was searched, nothing matched the
                # filter, and the row failed with "no matches found" - which hid a
                # typo behind a message about the query. Now that the name is handed
                # to genomepy, it would raise instead and end the whole run.
                if provider_filter and provider_filter not in PROVIDERS:
                    reason = (
                        f"unknown provider '{provider_filter}' "
                        f"(known: {', '.join(sorted(PROVIDERS))})"
                    )
                    print(f"ERROR: {reason}.", file=sys.stderr)
                    unresolved.append((lineno, species, query, reason))
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
                # A changed provider or release is a cache miss rather than a
                # mismatch to report: both are part of what identifies the request
                # (see load_cache and cached_resolution).
                cached_row = cached_resolution(
                    resolved_cache, species, query, release, provider_filter
                )
                if cached_row is not None:
                    cached_assembly = (cached_row.get("assembly") or "").strip()
                    cached_provider = (cached_row.get("provider") or "").strip()

                    print(
                        f"\n--- '{query}' for {species} already resolved to {cached_assembly}. "
                        "Skipping search. ---"
                    )
                    reason = ledger.conflict(
                        lineno, cached_provider, cached_assembly, release, aliases
                    )
                    if reason:
                        print(f"ERROR: {reason}.", file=sys.stderr)
                        unresolved.append((lineno, species, query, reason))
                        continue
                    # Rebuild the row rather than mutating the cached one, so stale or
                    # unexpected columns in an older resolved CSV are not carried through.
                    resolved_rows.append({
                        "species": species,
                        "assembly": cached_assembly,
                        "provider": cached_provider,
                        "annotation": annotation,
                        "query": query,
                        "release": release,
                        "aliases": ";".join(aliases),
                    })
                    write_resolved(output_file, resolved_rows)
                    continue

                # ------------------------------------------------------------------
                # Network search
                # ------------------------------------------------------------------
                print(f"\n--- Resolving query: '{query}' for {species} ---")

                results, offline = search_assemblies(query, provider_filter)

                # The search was already limited to the named provider; this keeps
                # the guarantee here, where the choice is made, rather than resting
                # it on how genomepy treats the argument.
                if provider_filter:
                    results = [r for r in results if r.provider.lower() == provider_filter]

                # With no provider named, the candidates below come from whichever
                # providers answered. Say when one did not: the list may be missing
                # the very assembly that was wanted, and a single remaining match
                # would otherwise be reported as "the exact match".
                if offline and not provider_filter:
                    warning = (
                        f"WARNING: {', '.join(sorted(offline))} could not be searched "
                        f"(offline), so the candidates for '{query}' may be incomplete."
                    )
                    # Red on a terminal, to match genomepy's warnings above it;
                    # plain when stderr is a file or a pipe.
                    if sys.stderr.isatty():
                        warning = f"\033[31m{warning}\033[0m"
                    print(warning, file=sys.stderr)

                if not results:
                    # Distinguish "the provider is down" from "the query matched
                    # nothing" before reporting - the two need opposite responses.
                    # Offline DURING the search is the evidence that counts (see
                    # search_assemblies); the ping covers a failure genomepy did not
                    # announce.
                    if provider_filter and (
                        provider_filter in offline
                        or provider_reachable(provider_filter) is False
                    ):
                        reason = (
                            f"provider '{provider_filter}' is unreachable, so '{query}' "
                            "could not be resolved. This is usually transient; re-run "
                            "when the provider is back. The request itself may be fine"
                        )
                    else:
                        reason = f"no matches found for '{query}'"
                        if provider_filter:
                            reason += f" with provider '{provider_filter}'"
                        elif offline:
                            reason += (
                                f" among the providers that answered "
                                f"({', '.join(sorted(offline))} could not be searched)"
                            )
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

                # An annotation the provider does not publish cannot be delivered, and
                # the search result already says so. Without this check the genome was
                # downloaded in full first - genomepy logs a warning and carries on -
                # and the request only failed at ingest, with "no GTF found".
                #
                # Not applied on a cache hit: the cached row records which assembly was
                # chosen, not whether it is annotated, and looking that up would cost
                # the network search the cache exists to avoid. A cached row switched
                # to annotation=true still fails at ingest, as before.
                if wants_annotation(annotation) and not chosen.has_annotation:
                    reason = (
                        f"annotation was requested, but genomepy reports none for "
                        f"{chosen.assembly} on '{chosen.provider}'. Set annotation to "
                        "false for this row, or choose an assembly or provider that "
                        "publishes one"
                    )
                    print(f"ERROR: {reason}.", file=sys.stderr)
                    unresolved.append((lineno, species, query, reason))
                    continue

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

                reason = ledger.conflict(
                    lineno, chosen.provider, chosen.assembly, release, aliases
                )
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
                    "aliases": ";".join(aliases),
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
