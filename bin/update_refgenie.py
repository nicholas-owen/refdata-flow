#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import sys
import os
import csv
import glob
import json
import logging
import subprocess
import shutil

import yaml


# Metadata files written into the download directory that must outlive it.
# Everything else genomepy leaves there is either the sequence itself or an index
# refgenie regenerates in the vault (*.fai, *.gzi, *.sizes, *.gaps.bed), so only
# these carry information that would otherwise be lost:
#   provenance.json     - written by bin/download.py: what was requested, when
#   README.txt          - written by genomepy: the source URLs
#   assembly_report.txt - written by genomepy: sequence-name mappings between
#                         UCSC / Ensembl / GenBank accessions, needed later to
#                         reconcile contig naming between a vault and an annotation
PROVENANCE_FILES = ("provenance.json", "README.txt", "assembly_report.txt")


def read_annotation_release(assembly_dir):
    """
    The annotation release recorded by bin/download.py, as (value, kind).

    Used as the refgenie asset tag. The genome digest is computed from the FASTA
    alone, so without a tag two annotation releases of the same assembly overwrite
    one another under 'default' and `refgenie seek` returns a confident path to
    whichever landed last.

    download.py identifies the release because that is where README.txt and the
    network are both available; this side only reads the result, so a vault can be
    rebuilt from an existing download directory without any network access.

    Returns (None, None) when the release could not be identified - the caller then
    keeps refgenie's default tag rather than failing an otherwise good ingest.
    """
    path = os.path.join(assembly_dir, "provenance.json")
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        # Absent for downloads made before release tagging existed, and unreadable
        # is treated the same way: fall back rather than block the build.
        return None, None

    release = data.get("annotation_release")
    return (release, data.get("annotation_release_kind")) if release else (None, None)

# Which refgenie recipe registers an annotation from which provider.
#
# These are provenance claims, not interchangeable containers, so the mapping is
# explicit and has no default. Previously the code read
#     "ensembl_gtf" if provider == "ensembl" else "gencode_gtf"
# which filed every non-Ensembl annotation under GENCODE - including UCSC yeast,
# despite GENCODE only ever publishing human and mouse. The vault then asserted a
# source that cannot exist.
#
# The two recipes also differ in what they do. `ensembl_gtf` derives
# <genome>_ensembl_TSS.bed and <genome>_ensembl_gene_body.bed, prefixing contigs
# with "chr" because Ensembl names them 1/2/3 and awk-ing Ensembl-specific
# attribute columns; run against an already-chr-prefixed UCSC GTF it would emit
# "chrchr1". `gencode_gtf` is a bare copy. So the old fallback was not merely
# mislabelled - it was the only one of the two that did not corrupt the output.
#
# refgenie 0.13 ships no ucsc_gtf or ncbi_gtf recipe, but `refgenie build --recipe`
# accepts a path to a JSON recipe file as well as a built-in name (see
# refgenie/refgenie.py, which json-loads a .json argument and takes the recipe's
# "name" from it). So providers without a built-in recipe get a passthrough recipe
# of our own under recipes/, named for the real source. Those are deliberately
# copy-only: the derivations in ensembl_gtf assume Ensembl contig naming and
# attribute ordering and would corrupt a UCSC or NCBI GTF.
ANNOTATION_RECIPES = {
    "ensembl": "ensembl_gtf",
    "gencode": "gencode_gtf",
    "ucsc":    "ucsc_gtf",
    "ncbi":    "ncbi_gtf",
}

# Recipes refgenie already knows by name; anything else is supplied as a file from
# recipes/. gencode_gtf is used as shipped - it is a passthrough, and GENCODE is
# chr-prefixed natively, so nothing needs correcting.
#
# ensembl_gtf is NOT used as shipped. The stock recipe pipes its derived BEDs
# through `sed 's/^/chr/'`, because refgenie's Ensembl recipes are written to feed
# pipelines that assume a UCSC- or GENCODE-named genome - the case where an Ensembl
# GTF is paired with a non-Ensembl FASTA. refdata-flow takes both from the same
# provider, so an Ensembl genome has contigs named 1/2/X/MT and that rewrite yields
# ensembl_tss / ensembl_gene_body BEDs referencing chr1/chr2, which match nothing in
# the vault. recipes/ensembl_gtf.json is the stock recipe with those two calls
# removed; the GTF asset itself is unchanged.
BUILTIN_RECIPES = {"gencode_gtf"}

# Custom recipe JSON lives beside the project, not inside bin/.
RECIPES_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, "recipes"
)


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(outdir: str) -> logging.Logger:
    """
    Configure a logger that writes to both stdout and a persistent log file
    in outdir. Console output is INFO level; the file captures full DEBUG
    detail for human review after a run.
    """
    logger = logging.getLogger("update_refgenie")
    logger.setLevel(logging.DEBUG)

    fmt = logging.Formatter(
        fmt="%(asctime)s [%(levelname)-8s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S"
    )

    # Stdout handler - visible in terminal and captured by SLURM/Nextflow logs
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(fmt)

    # File handler - append mode so multiple runs accumulate in one audit trail
    log_path = os.path.join(outdir, "update_refgenie.log")
    fh = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    logger.addHandler(console)
    logger.addHandler(fh)
    return logger


# ---------------------------------------------------------------------------
# Path safety
# ---------------------------------------------------------------------------

def safe_path_under(base: str, *parts: str) -> str:
    """
    Construct a path from parts and assert it resolves inside base.

    Both paths are canonicalised with os.path.realpath() to resolve symlinks
    and any '..' components before the containment check. Raises ValueError
    if the resolved path escapes the base directory, preventing path traversal
    from user-supplied CSV values reaching filesystem operations.
    """
    base_real = os.path.realpath(base)
    candidate = os.path.realpath(os.path.join(base, *parts))
    if not candidate.startswith(base_real + os.sep):
        raise ValueError(
            f"Path traversal detected: '{candidate}' is outside base directory '{base_real}'"
        )
    return candidate


# ---------------------------------------------------------------------------
# Safe directory removal
# ---------------------------------------------------------------------------

def logged_rmtree(path: str, logger: logging.Logger) -> bool:
    """
    Remove a directory tree, logging every individual failure via the onerror
    callback. Unlike ignore_errors=True, no failure is silently discarded -
    each problem is recorded with the affected path, failed operation, and
    exception detail.

    Returns True if the removal completed without errors, False otherwise.
    """
    errors = []

    def on_error(func, failed_path, excinfo):
        exc_type, exc_value, _ = excinfo
        errors.append((failed_path, exc_type.__name__, str(exc_value)))
        logger.error(
            "Cleanup failed | path='%s' | operation='%s' | error=%s: %s",
            failed_path, func.__name__, exc_type.__name__, exc_value
        )

    shutil.rmtree(path, onerror=on_error)

    if errors:
        log_path = os.path.join(os.path.dirname(path), "update_refgenie.log")
        logger.warning(
            "Cleanup of '%s' completed with %d error(s) - manual review required. "
            "Full details in: %s",
            path, len(errors), log_path
        )
        return False

    logger.info("Cleanup of '%s' completed successfully", path)
    return True


# ---------------------------------------------------------------------------
# Provenance preservation
# ---------------------------------------------------------------------------

def preserve_provenance(asm_dir: str, outdir: str, provider: str, species: str,
                        assembly: str, logger: logging.Logger) -> None:
    """
    Copy the run's provenance artefacts out of the download directory before that
    directory is removed.

    The files listed in PROVENANCE_FILES live in the published assembly directory
    and are not carried into the refgenie vault by `refgenie build`. Removing the
    directory therefore destroyed the only record of where the data came from -
    the audit trail the README advertises.

    They are copied to <outdir>/provenance/<provider>/<species>/<assembly>/, which
    is outside the download tree and is never removed by cleanup. Relocating these
    into the vault's per-digest asset directory is deferred until the digest is
    resolved earlier in the flow (roadmap Phase 4).
    """
    dest = safe_path_under(outdir, "provenance", provider, species, assembly)
    os.makedirs(dest, exist_ok=True)

    preserved = []
    for name in PROVENANCE_FILES:
        src = os.path.join(asm_dir, name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dest, name))
            preserved.append(name)

    if preserved:
        logger.info("Preserved provenance for '%s' (%s) in '%s'",
                    assembly, ", ".join(preserved), dest)
    else:
        logger.warning(
            "No provenance artefacts found in '%s' - nothing to preserve for '%s'",
            asm_dir, assembly
        )


# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------

def merge_aliases(assembly, assembly_safe, species, query, query_map,
                  custom_aliases_map, config_file, logger):
    """
    Merge this request's aliases into the digest's alias list in the refgenie config.

    Returns a list of problem descriptions (empty on success).

    Ordering matters and is the reason this is a separate step. refgenie links an
    asset into <vault>/alias/<alias>/<asset>/ only while a build actually executes,
    and only for the aliases the digest holds at that moment. A build that is
    skipped because its target already exists does not re-link. So if aliases are
    merged after the assets are built - as this script used to do - a freshly built
    annotation is linked under the primary alias alone, and no later run repairs it:
    `refgenie seek yeast/gencode_gtf` then returns a confident, well-formed path to
    a file that does not exist.

    Aliases are therefore merged immediately after the FASTA build (which is what
    creates the digest in the first place) and before anything else is built.
    """
    problems = []

    aliases_to_set = set([assembly_safe, assembly, assembly.split(".")[0]])

    # Look up by (species, query), so a species with more than one requested
    # assembly gets each row's own aliases rather than the last row's.
    key = (species, query)
    raw_query = query_map.get(key)
    if raw_query:
        aliases_to_set.add(raw_query)

    for ca in custom_aliases_map.get(key, []):
        aliases_to_set.add(ca)

    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)

        # Resolve the digest from this assembly's OWN identifiers rather than from
        # any overlap with aliases_to_set. Matching on overlap meant a custom alias
        # shared between two rows - 'yeast' on both the UCSC and the Ensembl yeast,
        # say - resolved to whichever digest was written first, and this row's
        # aliases were then merged onto the wrong genome.
        own_names = {assembly_safe, assembly, assembly.split(".")[0]}
        target_digest = None
        alias_owner = {}
        if "genomes" in config:
            for digest, data in config["genomes"].items():
                existing_aliases = data.get("aliases", []) or []
                for existing in existing_aliases:
                    alias_owner[existing] = digest
                if target_digest is None and any(a in existing_aliases for a in own_names):
                    target_digest = digest

        # An alias names exactly one genome. If one of ours already belongs to a
        # different digest, taking it would silently repoint it - so report it and
        # leave it where it is. The rest of this row's aliases still apply.
        if target_digest:
            conflicts = sorted(
                a for a in aliases_to_set
                if alias_owner.get(a) not in (None, target_digest)
            )
            if conflicts:
                logger.error(
                    "Alias(es) %s already name a different genome, so they were NOT "
                    "applied to '%s'. An alias cannot point at two assemblies; give "
                    "this row its own alias in the request CSV.",
                    ", ".join(f"'{c}'" for c in conflicts), assembly
                )
                problems.append(
                    f"{assembly}: alias(es) {', '.join(conflicts)} already belong to "
                    "another genome and were not applied"
                )
                aliases_to_set -= set(conflicts)

        if target_digest:
            current_aliases = config["genomes"][target_digest].get("aliases", [])
            merged_aliases  = sorted(list(set(current_aliases) | aliases_to_set))
            config["genomes"][target_digest]["aliases"] = merged_aliases
            logger.info(
                "Applying %d aliases to digest %s...: %s",
                len(merged_aliases), target_digest[:8], merged_aliases
            )
            with open(config_file, "w") as f:
                yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        else:
            logger.warning(
                "Could not resolve a digest for '%s' in config - skipping alias update",
                assembly
            )
            problems.append(f"{assembly}: no digest found in config, aliases not applied")
    except Exception as e:
        logger.error("Failed to update aliases via YAML for '%s': %s", assembly, e)
        problems.append(f"{assembly}: alias update failed ({e})")

    return problems


def relink_aliases(config_file, assembly, logger):
    """
    Refresh the alias symlink tree for every asset of this genome.

    refgenie links an asset into <vault>/alias/<alias>/<asset>/ only while a build
    actually executes, and only for the aliases the digest holds at that moment.
    The FASTA is necessarily built before its aliases exist - the digest it is keyed
    on is computed from that very FASTA - so its links start out covering the primary
    alias alone, and a later run does not repair them because the build is skipped.

    The obvious repair, re-running `refgenie build <genome>/fasta`, is actively
    harmful: the fasta recipe re-runs initialize_genome/set_genome_alias, which
    RESETS the alias list to the primary and undoes the merge. So the relink is done
    against refgenconf directly, which touches the symlinks and nothing else.

    Failure is logged, not treated as a problem: the data is in the vault either way,
    only the convenience shortcuts are affected.
    """
    try:
        from refgenconf import RefGenConf

        # RefGenConf.__init__ takes a dict of entries; from_yaml_file is the loader.
        rgc = RefGenConf.from_yaml_file(config_file)
        digest = rgc.get_genome_alias_digest(assembly, fallback=True)

        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
        assets = sorted(
            (config.get("genomes", {}).get(digest, {}).get("assets") or {}).keys()
        )

        for asset in assets:
            try:
                tag = rgc.get_default_tag(digest, asset)
            except Exception:
                tag = "default"
            rgc._symlink_alias(digest, asset, tag)

        logger.info(
            "Refreshed alias links for %s across assets: %s",
            assembly, ", ".join(assets) or "(none)"
        )
    except Exception as e:
        logger.warning(
            "Could not refresh alias links for '%s': %s. Assets are in the vault; some "
            "alias shortcuts may be missing.", assembly, e
        )


# ---------------------------------------------------------------------------
# Per-assembly processing
# ---------------------------------------------------------------------------

def process_assembly(row, outdir, config_file, refgenie_bin, query_map,
                     custom_aliases_map, do_cleanup, logger):
    """
    Register one resolved request into the refgenie vault.

    Returns a list of problem descriptions; an empty list means the row was fully
    ingested. Problems are collected rather than raised so one bad assembly does
    not abandon the rest of the batch, while still being reported and reflected in
    the process exit code.
    """
    problems = []

    species    = (row.get("species")    or "").strip()
    assembly   = (row.get("assembly")   or "").strip()
    provider   = (row.get("provider")   or "").strip()
    annotation = (row.get("annotation") or "false").strip().lower() == "true"
    # Identifies which raw request produced this row, so its aliases can be looked
    # up without colliding with another assembly of the same species. Absent from
    # resolved CSVs written before the column existed, which degrades to the old
    # behaviour for that row rather than failing.
    query      = (row.get("query")      or "").strip()

    logger.info("Processing %s...", assembly)

    # -------------------------------------------------------
    # Path traversal guard
    # Validate that the assembly directory resolves inside outdir
    # before any filesystem operation (glob, build, rmtree).
    # -------------------------------------------------------
    try:
        asm_dir = safe_path_under(outdir, provider, species, assembly)
    except ValueError as e:
        logger.error("Skipping %s - %s", assembly, e)
        return [f"{assembly}: {e}"]

    assembly_safe = assembly.replace(".", "_")

    # Track what actually made it into the vault. Cleanup is gated on these:
    # source data must never be removed on the strength of an ingest that only
    # partially succeeded.
    fasta_registered = False
    gtf_registered = False
    # Aliases are normally merged mid-flow, right after the FASTA build. This records
    # whether that happened so the fallback at the end does not repeat it.
    aliases_merged = False
    fasta_path = None

    if os.path.exists(asm_dir):
        # Find and register FASTA
        fasta_files = (
            glob.glob(os.path.join(asm_dir, "*.fa")) +
            glob.glob(os.path.join(asm_dir, "*.fa.gz"))
        )
        if fasta_files:
            fasta_path = os.path.abspath(fasta_files[0])
            logger.info("Building FASTA digest in refgenie for '%s'", assembly_safe)
            cmd = [
                sys.executable, refgenie_bin, "build",
                f"{assembly_safe}/fasta",
                "--files", f"fasta={fasta_path}",
                "-c", config_file
            ]
            try:
                subprocess.run(cmd, check=True)
                fasta_registered = True
            except subprocess.CalledProcessError as e:
                logger.error("FASTA build failed for '%s': %s", assembly, e)
                problems.append(f"{assembly}: FASTA build failed ({e})")
        else:
            logger.error("No FASTA found for '%s' in '%s'", assembly, asm_dir)
            problems.append(f"{assembly}: no FASTA found in '{asm_dir}'")

        # -------------------------------------------------------
        # Aliases, before any further asset is built.
        # The FASTA build above is what created the digest, so this is the first
        # moment the aliases can be attached - and it must happen before the
        # annotation build, or that asset is linked under the primary alias only.
        # See merge_aliases() for the full explanation.
        # -------------------------------------------------------
        if fasta_registered:
            problems.extend(merge_aliases(
                assembly, assembly_safe, species, query, query_map,
                custom_aliases_map, config_file, logger
            ))
            aliases_merged = True

        # Find and register annotation
        if annotation:
            recipe_name = ANNOTATION_RECIPES.get(provider.lower())
            gtf_files = (
                glob.glob(os.path.join(asm_dir, "*.gtf")) +
                glob.glob(os.path.join(asm_dir, "*.gtf.gz"))
            )
            recipe_file = None
            if recipe_name and recipe_name not in BUILTIN_RECIPES:
                recipe_file = os.path.join(RECIPES_DIR, f"{recipe_name}.json")
                if not os.path.isfile(recipe_file):
                    logger.error(
                        "Provider '%s' maps to custom recipe '%s' but '%s' is missing",
                        provider, recipe_name, recipe_file
                    )
                    problems.append(
                        f"{assembly}: custom recipe file not found ({recipe_file})"
                    )
                    recipe_name = None

            if recipe_name is None:
                # Refuse rather than assert a provenance that is not true: filing this
                # annotation under another provider's recipe would record a false source.
                logger.error(
                    "Annotation requested for '%s' but there is no annotation recipe for "
                    "provider '%s' (known: %s). The GTF was downloaded and is retained in "
                    "'%s', but it has NOT been registered in the vault. Set "
                    "annotation=false for this row, resolve it against a known provider, "
                    "or add a recipe for '%s' under recipes/.",
                    assembly, provider, ", ".join(sorted(ANNOTATION_RECIPES)), asm_dir,
                    provider.lower()
                )
                problems.append(
                    f"{assembly}: no annotation recipe for provider '{provider}' "
                    f"(known: {', '.join(sorted(ANNOTATION_RECIPES))})"
                )
            elif gtf_files:
                gtf_path = os.path.abspath(gtf_files[0])
                release, release_kind = read_annotation_release(asm_dir)
                # An unidentified release keeps refgenie's own default, which is the
                # behaviour that existed before release tagging. The warning was
                # already issued at download time, where the cause is known.
                target = f"{assembly_safe}/{recipe_name}"
                if release:
                    target += f":{release}"
                logger.info(
                    "Building GTF annotation digest in refgenie for '%s' (provider: %s, "
                    "recipe: %s%s, tag: %s%s)", assembly_safe, provider, recipe_name,
                    " [custom]" if recipe_file else "",
                    release or "default",
                    f" [{release_kind}]" if release_kind else " [release not identified]"
                )
                cmd = [
                    sys.executable, refgenie_bin, "build",
                    target,
                    "--files", f"{recipe_name}={gtf_path}",
                    "-c", config_file
                ]
                if recipe_file:
                    cmd += ["--recipe", os.path.abspath(recipe_file)]
                try:
                    subprocess.run(cmd, check=True)
                    gtf_registered = True
                except subprocess.CalledProcessError as e:
                    logger.error("GTF build failed for '%s': %s", assembly, e)
                    problems.append(f"{assembly}: GTF build failed ({e})")

                # Building a tag does NOT make it the default - verified against
                # refgenie 0.13.0, where a second tag registers alongside the first
                # and `refgenie seek <genome>/<asset>` keeps returning the older one.
                # Without this step a newly ingested release would be present in the
                # vault but invisible to every caller that does not name a tag, which
                # is the opposite of the intended "latest wins, older stays reachable"
                # behaviour. Not fatal: the asset is built and addressable by tag.
                if gtf_registered and release:
                    try:
                        subprocess.run(
                            [sys.executable, refgenie_bin, "tag",
                             f"{assembly_safe}/{recipe_name}:{release}",
                             "--default", "-c", config_file],
                            check=True
                        )
                    except subprocess.CalledProcessError as e:
                        logger.warning(
                            "Built '%s/%s:%s' but could not make it the default tag: %s. "
                            "The asset is addressable by tag; callers using no tag will "
                            "still resolve to the previous release.",
                            assembly_safe, recipe_name, release, e
                        )
            else:
                # Annotation was explicitly requested, so its absence is a failure of
                # the run, not a cosmetic warning.
                logger.error(
                    "Annotation requested for '%s' but no GTF found in '%s'",
                    assembly, asm_dir
                )
                problems.append(f"{assembly}: annotation requested but no GTF found")

        # Every asset for this genome now exists and the aliases are merged, so
        # refresh the symlink tree once, covering assets built before the aliases
        # were attached (the FASTA in particular).
        if aliases_merged:
            relink_aliases(config_file, assembly_safe, logger)

        # -------------------------------------------------------
        # Cleanup - opt-in, and only after a verified ingest
        #
        # Removing the downloaded source is a retention decision, not a tidy-up:
        # upstream providers archive or reorganise old releases, so a re-download is
        # not guaranteed to be possible later. It therefore happens only when the
        # user asks for it (--cleanup) AND everything that was supposed to reach the
        # vault actually did.
        # -------------------------------------------------------
        ingest_complete = fasta_registered and (gtf_registered or not annotation)

        if not do_cleanup:
            logger.info(
                "Retaining source directory '%s' (cleanup is opt-in - pass --cleanup "
                "to remove sources after a verified ingest)", asm_dir
            )
        elif not ingest_complete:
            logger.warning(
                "NOT removing '%s' - ingest incomplete "
                "(fasta_registered=%s, annotation_requested=%s, gtf_registered=%s). "
                "Source data is retained so nothing is lost.",
                asm_dir, fasta_registered, annotation, gtf_registered
            )
        else:
            preserve_provenance(asm_dir, outdir, provider, species, assembly, logger)
            logger.info("Cleaning up download directory: '%s'", asm_dir)
            if not logged_rmtree(asm_dir, logger):
                logger.warning(
                    "Directory '%s' may require manual cleanup. See '%s' for full details.",
                    asm_dir, os.path.join(outdir, "update_refgenie.log")
                )
                problems.append(f"{assembly}: cleanup did not complete cleanly")
    else:
        # This is the normal path on a re-run after a previous --cleanup, but it is
        # also what a provider/species case mismatch looks like, because the download
        # and ingest stages construct this path independently. The message no longer
        # asserts that the genome is present.
        logger.info(
            "Directory '%s' not found - nothing to ingest for '%s'. If this assembly was "
            "not registered by an earlier run, check that the provider and species in the "
            "resolved CSV match the published directory layout.", asm_dir, assembly
        )

    # -----------------------------------------------------------
    # Alias fallback
    #
    # The normal path merges aliases immediately after the FASTA build, so that
    # every asset built afterwards is linked against the full alias set. This only
    # runs when that did not happen - the download directory was absent (a re-run
    # after --cleanup), or the FASTA did not register - so an existing vault entry
    # still picks up aliases added to the request CSV since it was built.
    # -----------------------------------------------------------
    if not aliases_merged:
        problems.extend(merge_aliases(
            assembly, assembly_safe, species, query, query_map,
            custom_aliases_map, config_file, logger
        ))

    return problems


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Register downloaded reference data into a refgenie vault."
    )
    parser.add_argument("raw_csv", help="The original request CSV (source of queries and aliases)")
    parser.add_argument("resolved_csv", help="The resolved request CSV produced by bin/resolve.py")
    parser.add_argument("outdir", help="Output directory holding the refgenie vault")
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove each downloaded source directory after its contents have been "
             "verified into the vault. Off by default: upstream providers archive old "
             "releases, so a deleted download may not be re-fetchable.",
    )
    args = parser.parse_args()

    raw_csv      = args.raw_csv
    resolved_csv = args.resolved_csv
    outdir       = args.outdir

    # Ensure outdir exists before we try to create the log file inside it
    os.makedirs(outdir, exist_ok=True)

    logger      = setup_logger(outdir)
    config_file = os.path.join(outdir, "refgenie_export.yaml")
    refgenie_bin = os.path.join(os.path.dirname(sys.executable), "refgenie")

    logger.info("update_refgenie starting | outdir='%s' | cleanup=%s", outdir, args.cleanup)

    # Ensure refgenie config is initialized
    if not os.path.exists(config_file):
        logger.info("Initializing refgenie config at '%s'", config_file)
        subprocess.run(
            [sys.executable, refgenie_bin, "init", "-c", config_file],
            check=True
        )

    # Map each raw request to its custom aliases.
    #
    # Keyed by (species, query), NOT by species alone. Keying on species meant the
    # last row for a species overwrote every earlier one, so two assemblies of the
    # same species - Saccharomyces_cerevisiae as both UCSC sacCer3 and Ensembl
    # R64-1-1, say - collapsed into one entry. The consequences were silent and
    # severe: the second row's aliases were applied to the first row's digest, which
    # bound 'R64-1-1' to the UCSC genome; refgenie then resolved R64-1-1/fasta to
    # that digest, skipped the build because the flag already existed, and filed the
    # Ensembl annotation under a UCSC assembly. The run still exited 0.
    #
    # (species, query) is the right key because query is what distinguishes two rows
    # of the same species, and bin/resolve.py carries it into the resolved CSV - so
    # the row being ingested can look up exactly the request that produced it.
    query_map          = {}
    custom_aliases_map = {}
    try:
        with open(raw_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sp          = (row.get("species") or "").strip()
                q           = (row.get("query") or "").strip()
                aliases_str = (row.get("aliases") or "").strip()
                if sp and q:
                    key = (sp, q)
                    query_map[key] = q
                    if aliases_str:
                        custom_aliases_map[key] = [
                            a.strip() for a in aliases_str.split(";") if a.strip()
                        ]
    except Exception as e:
        logger.warning("Could not read '%s': %s", raw_csv, e)

    # Process each resolved request. Each row is isolated so a failure part-way
    # through the batch no longer abandons the assemblies after it.
    all_problems = []
    processed = 0
    try:
        with open(resolved_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            for lineno, row in enumerate(reader, start=2):
                assembly = (row.get("assembly") or "").strip() or f"row {lineno}"
                processed += 1
                try:
                    all_problems.extend(
                        process_assembly(
                            row, outdir, config_file, refgenie_bin,
                            query_map, custom_aliases_map, args.cleanup, logger
                        )
                    )
                except Exception as e:
                    logger.error("Unhandled error processing '%s': %s", assembly, e)
                    all_problems.append(f"{assembly}: unhandled error ({e})")
    except Exception as e:
        logger.error("Fatal error reading '%s': %s", resolved_csv, e)
        sys.exit(1)

    # ----------------------------------------------------------------------
    # Report. A run that failed to register something the user asked for now exits
    # non-zero rather than printing an unconditional success message.
    # ----------------------------------------------------------------------
    if all_problems:
        logger.error(
            "Refgenie update finished with %d problem(s) across %d request(s):",
            len(all_problems), processed
        )
        for problem in all_problems:
            logger.error("  - %s", problem)
        logger.error(
            "Source data for incomplete ingests has been retained. "
            "See '%s' for the full log.", os.path.join(outdir, "update_refgenie.log")
        )
        sys.exit(1)

    logger.info("Refgenie config updated successfully for %d request(s)!", processed)


if __name__ == "__main__":
    main()
