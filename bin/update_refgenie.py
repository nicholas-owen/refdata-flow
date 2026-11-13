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
import re
import subprocess
import shutil
import tempfile

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
    # An annotation given by address, for a genome given by address. No provider
    # stands behind it, and the asset name says that and no more.
    "url":     "url_gtf",
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

# The genome itself is built with our own fasta recipe too: refgenie's, plus a
# sequence dictionary. See recipes/fasta.json for why the dictionary belongs
# inside this asset. refgenie keys its genome-initialising behaviour on the
# recipe's NAME, which is still 'fasta', so supplying it from a file changes
# nothing about how the digest is computed.
FASTA_RECIPE = os.path.join(RECIPES_DIR, "fasta.json")


def fasta_has_dictionary(config_file, vault_name):
    """
    Does this genome's fasta asset include a sequence dictionary?

    True / False, or None when the genome or its fasta asset is not in the config.

    A genome ingested before recipes/fasta.json existed has a fasta asset built by
    refgenie's stock recipe, with no dictionary. `refgenie build` will not add one
    on its own: the asset's completion flag is already there, so the build is
    skipped whatever the recipe now says. The caller uses this to decide whether
    to force the rebuild.

    The answer comes from the FILE, not from the config's seek keys. A skipped
    build still writes the new recipe's seek keys into the config, so after one the
    config lists `dict` for an asset directory that holds no such file, and
    `refgenie seek <genome>/fasta.dict` returns a path to nothing. Trusting the
    seek key here would mean never repairing exactly the genomes that need it.
    """
    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError):
        return None
    # refgenie records where the vault is; a config that has been moved by hand
    # along with its vault may still carry the old location, so fall back to the
    # directory the config is in, which is where this pipeline always puts it.
    config_dir = os.path.dirname(os.path.abspath(config_file))
    roots = [config.get("genome_folder") or config_dir, config_dir]
    for digest, data in (config.get("genomes") or {}).items():
        if vault_name in ((data or {}).get("aliases") or []):
            fasta = (data.get("assets") or {}).get("fasta")
            if not fasta:
                return None
            tag = fasta.get("default_tag") or "default"
            seek_keys = ((fasta.get("tags") or {}).get(tag) or {}).get("seek_keys") or {}
            filename = seek_keys.get("dict")
            if not filename:
                return False
            return any(
                os.path.isfile(os.path.join(root, "data", digest, "fasta", tag, filename))
                for root in roots
            )
    return None


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
                        assembly: str, logger: logging.Logger,
                        release_parts: tuple = ()) -> None:
    """
    Copy the run's provenance artefacts out of the download directory before that
    directory is removed.

    The files listed in PROVENANCE_FILES live in the published assembly directory
    and are not carried into the refgenie vault by `refgenie build`. Removing the
    directory therefore destroyed the only record of where the data came from -
    the audit trail the README advertises.

    They are copied to <outdir>/provenance/<provider>/<species>/<assembly>/, which
    is outside the download tree and is never removed by cleanup. A pinned release
    adds the same release-<n>/ level the download tree has (`release_parts`), so
    two releases of one assembly do not overwrite each other's record here either. Relocating these
    into the vault's per-digest asset directory is deferred until the digest is
    resolved earlier in the flow (roadmap Phase 4).
    """
    dest = safe_path_under(outdir, "provenance", provider, species,
                           *release_parts, assembly)
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

def vault_genome_name(assembly, provider):
    """
    The name a genome is built under in refgenie: '<assembly>_<provider>'.

    It used to be the assembly name alone, which two providers can share - Ensembl
    and GENCODE both publish GRCm39 and GRCh38.p14. Reproduced against refgenie
    0.13.0 with two toy FASTAs of one name (contigs 1/2 and chr1/chr2):

      * `refgenie build <name>/fasta` for the second provider digests its FASTA and
        attaches the same name to the new digest, so one name now sits on two
        genomes;
      * `refgenie build <name>/gencode_gtf` then resolves that name to the FIRST
        digest, and files the chr-prefixed GENCODE annotation under the Ensembl
        genome whose contigs are 1/2/X.

    The run did exit non-zero, but for an alias conflict; nothing said an
    annotation had been attached to the wrong sequence, and
    `refgenie seek <name>/gencode_gtf` returned it without complaint.

    A build name unique to one (assembly, provider) means no build can resolve to
    another provider's genome. The bare assembly name is still offered as an
    alias, on a first-come basis (see merge_aliases).

    "." becomes "_" because refgenie's registry-path parser does not accept a dot
    in a genome name. The provider is lower-cased so 'Ensembl' from genomepy and
    'ensembl' from a hand-written CSV name one genome, not two.
    """
    provider_part = re.sub(r"[^a-z0-9]+", "_", str(provider).strip().lower()).strip("_")
    if not provider_part:
        raise ValueError(f"cannot name a vault genome for '{assembly}': no provider given")
    return f"{assembly.replace('.', '_')}_{provider_part}"


def parse_aliases(value):
    """'yeast; sacCer' -> ['yeast', 'sacCer']."""
    return [a.strip() for a in str(value or "").split(";") if a.strip()]


# ---------------------------------------------------------------------------
# Staged sources
#
# Run by hand, this script finds each download where main.nf publishes it, under
# <outdir>/<provider>/<species>/[release-<n>/]<assembly>/. Run as the REGISTER
# process it must not look there: Nextflow publishes in the background, so that
# copy may be half-written - or absent, when --cleanup turns publishing off. The
# process is handed each download as a staged input instead and names it here with
# --source, keyed on what makes a download unique.
# ---------------------------------------------------------------------------

def source_key(provider, assembly, release):
    """'ensembl/GRCh38.p14@116'; the release part is empty for an unpinned request."""
    return f"{str(provider).strip().lower()}/{assembly}@{release or ''}"


def parse_sources(values):
    """['ensembl/GRCm39@=dir', ...] -> {'ensembl/GRCm39@': '/abs/dir', ...}."""
    sources = {}
    for value in values or []:
        key, separator, path = value.partition("=")
        if not separator or not key or not path:
            raise ValueError(
                f"--source expects PROVIDER/ASSEMBLY@[RELEASE]=DIR, got '{value}'"
            )
        if key in sources:
            raise ValueError(f"--source names '{key}' twice")
        sources[key] = os.path.abspath(path)
    return sources


def source_mismatch(source_dir, provider, assembly, release):
    """
    Why this directory is not the download for this request, or None.

    A staged source is matched to its row by position in a list Nextflow builds, and
    nothing about a directory called 'GRCm39' says whose GRCm39 it is. Ingesting the
    wrong one files a genome under another's name - the failure this project keeps
    closing off - so the directory's own provenance.json is checked against the row
    before anything is built from it.

    Only fields the file actually has are compared: a download made before a field
    was recorded is accepted on the ones it does carry.
    """
    try:
        with open(os.path.join(source_dir, "provenance.json")) as fh:
            recorded = json.load(fh)
    except (OSError, ValueError):
        return None

    expected = {
        "assembly": assembly,
        "provider": str(provider).strip().lower(),
        "release_requested": release or "",
    }
    for field, want in expected.items():
        if field not in recorded:
            continue
        got = str(recorded[field] or "").strip()
        if field == "provider":
            got = got.lower()
        if got != want:
            return (
                f"the source directory '{source_dir}' records {field} "
                f"'{got or '(none)'}', not '{want or '(none)'}' - it is not the "
                "download for this request and has not been ingested"
            )
    return None


def vault_roots(config, config_file):
    """
    Directories that may hold the vault's data/ tree, most likely first.

    refgenie records the vault's location in the config; a config that was moved
    by hand along with its vault can still carry the old one. The directory the
    config sits in is where this pipeline always puts the vault, so it is tried too.
    """
    config_dir = os.path.dirname(os.path.abspath(config_file))
    recorded = config.get("genome_folder")
    return [recorded, config_dir] if recorded and recorded != config_dir else [config_dir]


def asset_tags_on_disk(config, config_file, digest, asset):
    """
    The tags of an asset that the config lists AND that have a directory in the
    vault.

    Both, because refgenie's config can run ahead of the disk (see
    fasta_has_dictionary), and a tag nobody can read is not a reason to skip work.
    """
    genome = (config.get("genomes") or {}).get(digest) or {}
    listed = ((genome.get("assets") or {}).get(asset) or {}).get("tags") or {}
    return {
        str(tag) for tag in listed
        if any(os.path.isdir(os.path.join(root, "data", digest, asset, str(tag)))
               for root in vault_roots(config, config_file))
    }


def asset_file(config, config_file, digest, asset, tag, seek_key):
    """
    The absolute path of one file of an asset in the vault's data/ tree, or ''.

    The data/ path, named by digest - not the alias path `refgenie seek` prints.
    It is what a derive process is handed: it does not move when an alias is
    added or removed. '' when the config does not list the file or it is not there.
    """
    genome = (config.get("genomes") or {}).get(digest) or {}
    tags = ((genome.get("assets") or {}).get(asset) or {}).get("tags") or {}
    # Compared as text: a tag like 116 can come back from YAML as a number or a
    # string depending on who wrote the file.
    tag_data = next((data for name, data in tags.items() if str(name) == str(tag)), None) or {}
    filename = (tag_data.get("seek_keys") or {}).get(seek_key)
    if not filename:
        return ""
    for root in vault_roots(config, config_file):
        candidate = os.path.join(root, "data", digest, asset, str(tag), filename)
        if os.path.isfile(candidate):
            return candidate
    return ""


MANIFEST_FIELDS = [
    "species", "assembly", "provider", "release", "vault_name", "digest",
    # What a derive process needs to build something FROM this genome:
    "fasta",             # the genome FASTA in the vault ('' if not there)
    "annotation_asset",  # e.g. ensembl_gtf; '' when no annotation was requested
    "annotation_tag",    # the release this row means: its pin, else the default tag
    "gtf",               # that release's GTF in the vault ('' if not there)
    "fasta_txome_tags",  # ';'-joined releases a transcriptome already exists for
]


def write_manifest(path, rows, config_file):
    """
    Write a TSV of the requested genomes that are in the vault.

    This is what a later stage needs to build something FROM a registered genome:
    the name to register under, the digest that identifies the sequence, and where
    the genome and the requested annotation release are. A row that did not reach
    the vault is left out, so nothing downstream is started on a genome that is not
    there; a file that is not on disk is an empty field, for the same reason.

    `annotation_tag` is what the row means by "the annotation": the pinned release
    if it has one, otherwise the asset's default tag, which the ingest keeps on the
    newest release held.
    """
    with open(config_file, "r") as f:
        config = yaml.safe_load(f) or {}
    digest_of = {}
    for digest, data in (config.get("genomes") or {}).items():
        for alias in ((data or {}).get("aliases") or []):
            digest_of[alias] = digest

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            assembly = (row.get("assembly") or "").strip()
            provider = (row.get("provider") or "").strip()
            pinned = (row.get("release") or "").strip()
            wants_annotation = (row.get("annotation") or "false").strip().lower() == "true"
            try:
                vault_name = vault_genome_name(assembly, provider)
            except ValueError:
                continue
            if vault_name not in digest_of:
                continue
            digest = digest_of[vault_name]

            annotation_asset = annotation_tag = gtf = ""
            recipe = ANNOTATION_RECIPES.get(provider.lower())
            if wants_annotation and recipe:
                assets = (config["genomes"][digest].get("assets") or {})
                tag = pinned or (assets.get(recipe) or {}).get("default_tag")
                if tag:
                    annotation_asset, annotation_tag = recipe, str(tag)
                    gtf = asset_file(config, config_file, digest, recipe, tag, recipe)

            writer.writerow({
                "species": (row.get("species") or "").strip(),
                "assembly": assembly,
                "provider": provider,
                "release": pinned,
                "vault_name": vault_name,
                "digest": digest,
                "fasta": asset_file(config, config_file, digest, "fasta", "default", "fasta"),
                "annotation_asset": annotation_asset,
                "annotation_tag": annotation_tag,
                "gtf": gtf,
                "fasta_txome_tags": ";".join(sorted(
                    asset_tags_on_disk(config, config_file, digest, "fasta_txome"))),
            })


def snapshot_aliases(config_file):
    """
    Every genome's alias list as it stands now: {digest: [aliases]}.

    Taken before a FASTA build, because the build destroys it. `refgenie build
    <name>/fasta` ends by setting the genome's aliases to <name> ALONE - even when
    the asset already exists and nothing is rebuilt. Reproduced against refgenie
    0.13.0: a genome ingested with the alias 'first_name', then ingested again by
    a request that asked for 'second_name', was left with only the second. The
    first was gone from the config, its symlink tree still on disk, and
    `refgenie seek first_name/fasta` failed. Nothing reported it.

    merge_aliases() puts back whatever this recorded for the genome it is working
    on, so an alias, once granted, stays until someone removes it on purpose.
    Keyed by digest, not by name, because the name a row builds under may be new
    to a genome that has been in the vault for a long time.
    """
    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError):
        return {}
    return {
        digest: list((data or {}).get("aliases") or [])
        for digest, data in (config.get("genomes") or {}).items()
    }


def merge_aliases(vault_name, assembly, auto_names, requested_aliases,
                  config_file, logger, aliases_before=None):
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

    The digest is found by `vault_name` and nothing else. That name is unique to one
    (assembly, provider), so the lookup cannot land on another provider's copy of
    the same assembly - which matching on the bare assembly name could, and did.

    Two kinds of alias, treated differently when one is already taken by another
    genome:

      * `requested_aliases` come from the `aliases` column. The user asked for them
        by name, so being unable to grant one is a problem and fails the run.
      * `auto_names` are derived - the assembly name, its stem, the query. Nobody
        asked for 'GRCm39' to mean this genome in particular, so when another
        provider's GRCm39 already holds it, this one is simply reachable under its
        vault name instead. Logged, not failed: requesting an assembly from two
        providers is legitimate and must not exit non-zero.
    """
    problems = []

    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}

        alias_owner = {}
        for digest, data in (config.get("genomes") or {}).items():
            for existing in ((data or {}).get("aliases") or []):
                alias_owner[existing] = digest

        target_digest = alias_owner.get(vault_name)
        if not target_digest:
            logger.warning(
                "No genome named '%s' in the config - skipping alias update for '%s'. "
                "Vault genomes are named <assembly>_<provider>; an entry built before "
                "that naming was introduced is picked up by re-ingesting its download.",
                vault_name, assembly
            )
            problems.append(
                f"{assembly}: no genome named '{vault_name}' in config, aliases not applied"
            )
            return problems

        def taken(alias):
            return alias_owner.get(alias) not in (None, target_digest)

        requested = set(requested_aliases)
        refused = sorted(a for a in requested if taken(a))
        if refused:
            logger.error(
                "Alias(es) %s already name a different genome, so they were NOT "
                "applied to '%s'. An alias cannot point at two assemblies; give "
                "this row its own alias in the request CSV.",
                ", ".join(f"'{c}'" for c in refused), assembly
            )
            problems.append(
                f"{assembly}: alias(es) {', '.join(refused)} already belong to "
                "another genome and were not applied"
            )

        derived = set(auto_names) - requested
        unavailable = sorted(a for a in derived if taken(a))
        if unavailable:
            logger.info(
                "Name(s) %s already refer to another genome in the vault (the same "
                "assembly from a different provider, usually) and stay with it. "
                "This genome is reachable as '%s'.",
                ", ".join(f"'{c}'" for c in unavailable), vault_name
            )

        # What this genome held before the build that has just reset its list (see
        # snapshot_aliases). They were this digest's a moment ago and nothing else
        # has run since, so none can have been taken by another genome; the check
        # is there for the case that someone edits the config between the two.
        restored = {
            a for a in (aliases_before or {}).get(target_digest, []) if not taken(a)
        }

        aliases_to_set = (
            {vault_name} | restored
            | (requested - set(refused)) | (derived - set(unavailable))
        )

        entry = config["genomes"][target_digest]
        merged_aliases = sorted(set(entry.get("aliases") or []) | aliases_to_set)
        entry["aliases"] = merged_aliases
        logger.info(
            "Applying %d aliases to digest %s...: %s",
            len(merged_aliases), target_digest[:8], merged_aliases
        )
        with open(config_file, "w") as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    except Exception as e:
        logger.error("Failed to update aliases via YAML for '%s': %s", assembly, e)
        problems.append(f"{assembly}: alias update failed ({e})")

    return problems


# ---------------------------------------------------------------------------
# Default tag
# ---------------------------------------------------------------------------

def current_default_tag(config_file, vault_name, asset):
    """
    The tag `refgenie seek <genome>/<asset>` resolves to right now, or None.

    None covers every "nothing to compare against" case alike - no such genome, no
    such asset yet, an unreadable config - because the caller's response to all of
    them is the same: let the newly built tag become the default.
    """
    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
        for data in (config.get("genomes") or {}).values():
            if vault_name in ((data or {}).get("aliases") or []):
                asset_data = (data.get("assets") or {}).get(asset) or {}
                return asset_data.get("default_tag")
    except Exception:
        pass
    return None


def release_sort_key(tag):
    """
    Order release tags the way the releases themselves are ordered, or None.

    Every tag download.py derives is numeric at heart: an Ensembl release ('116'),
    a GENCODE release ('44', or 'M33' for mouse), or an upstream modification date
    ('20240312'). Comparing the digits as integers orders all three correctly;
    comparing the strings would put '99' above '116'. A tag is only ever compared
    with another tag of the same asset on the same genome, so the three vocabularies
    never meet.

    None for a tag with no digits - refgenie's own 'default', in practice.
    """
    numbers = re.findall(r"\d+", str(tag or ""))
    return tuple(int(n) for n in numbers) if numbers else None


def should_become_default(current, new):
    """
    Should a newly built tag replace the current default?

    Yes unless it is demonstrably older. The ingest used to make every tag it built
    the default, so pinning an older release for one analysis - 110, with 116
    already in the vault - silently moved `refgenie seek <genome>/<asset>` back to
    110 for every caller that names no tag.

    An unparseable current tag ('default', from before release tagging) is replaced:
    a release that has been identified beats one that has not.
    """
    current_key, new_key = release_sort_key(current), release_sort_key(new)
    if current_key is None or new_key is None:
        return True
    return new_key >= current_key


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

# ---------------------------------------------------------------------------
# The download's record, in the vault
# ---------------------------------------------------------------------------

def build_passthrough(refgenie_bin, config_file, vault_name, asset, tag, source_file,
                      logger, force=False):
    """
    Build a one-file asset from recipes/<asset>.json. Returns a problem, or None.

    For the small files that describe a download rather than being one. The asset
    is built once and then left: `refgenie build` skips a tag that exists, which
    is right, because the vault's FASTA and GTF are also kept from the download
    that first supplied them, and a record should describe what is there. `force`
    is for the one case where the caller has just replaced what is there.

    The file is checked on disk afterwards. A skipped build still lists its seek
    keys in the config (see fasta_has_dictionary), so the config alone cannot say
    the asset exists; if the file is missing, the build is run again in earnest.
    """
    target = f"{vault_name}/{asset}:{tag}"
    recipe = os.path.join(RECIPES_DIR, f"{asset}.json")
    command = [
        sys.executable, refgenie_bin, "build", target,
        "--files", f"{asset}={os.path.abspath(source_file)}",
        "--recipe", os.path.abspath(recipe),
        "-c", config_file,
    ]
    previous_default = current_default_tag(config_file, vault_name, asset)

    def on_disk():
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
        for digest, data in (config.get("genomes") or {}).items():
            if vault_name in ((data or {}).get("aliases") or []):
                return bool(asset_file(config, config_file, digest, asset, tag, asset))
        return False

    try:
        subprocess.run(command + (["-N"] if force else []), check=True)
        if not on_disk():
            subprocess.run(command + ["-N"], check=True)
    except subprocess.CalledProcessError as e:
        logger.error("Could not build '%s': %s", target, e)
        return f"'{target}' was not built ({e})"
    if not on_disk():
        logger.error("'%s' is listed in the config but its file is not in the vault", target)
        return f"'{target}' is not in the vault after building"

    if previous_default not in (None, tag) and should_become_default(previous_default, tag):
        try:
            subprocess.run([sys.executable, refgenie_bin, "tag", target, "--default",
                            "-c", config_file], check=True)
        except subprocess.CalledProcessError as e:
            logger.warning("Built '%s' but could not make it the default tag: %s", target, e)
    return None


def register_download_record(asm_dir, vault_name, assembly, refgenie_bin, config_file,
                             logger, force=False):
    """
    File the download's own record in the vault, beside what it describes.

    Two assets, both passthroughs of files the download directory already holds:

      provenance       provenance.json. Tagged with the annotation release the
                       download carried, or 'default' for sequence alone.
      assembly_report  assembly_report.txt, where genomepy could fetch one: the
                       table of each sequence's name in every naming scheme.

    Until now these reached only <outdir>/provenance/<provider>/<species>/..., a
    path built from the request's own spelling, with no route to it from a digest
    and nothing in `refgenie list` to say it was there - and only when --cleanup
    ran. In the vault they are keyed by digest, linked under every alias, and
    found with `refgenie seek <genome>/provenance`.

    A failure here is reported as a problem. The genome is in the vault either
    way, but a reference manager that quietly drops its own records is the thing
    this pipeline exists not to be.
    """
    problems = []

    record = os.path.join(asm_dir, "provenance.json")
    if os.path.isfile(record):
        release, _kind = read_annotation_release(asm_dir)
        problem = build_passthrough(refgenie_bin, config_file, vault_name, "provenance",
                                    release or "default", record, logger, force)
        if problem:
            problems.append(f"{assembly}: {problem}")
    else:
        logger.warning("No provenance.json in '%s'; '%s' has no download record in the "
                       "vault.", asm_dir, vault_name)

    report = os.path.join(asm_dir, "assembly_report.txt")
    if os.path.isfile(report) and os.path.getsize(report) > 0:
        problem = build_passthrough(refgenie_bin, config_file, vault_name, "assembly_report",
                                    "default", report, logger, force)
        if problem:
            problems.append(f"{assembly}: {problem}")

    return problems


INDEX_FIELDS = [
    "digest", "aliases", "provider", "species", "assembly", "annotation_release",
    "assembly_accession", "mask", "retrieved_utc",
    "genome_content_sha256", "annotation_content_sha256", "assets",
]


def write_vault_index(outdir, config_file, logger):
    """
    Write <outdir>/provenance/index.tsv: what is in the vault, one line per record.

    The question asked most often of a shared vault is not "tell me about this
    digest" but "what is in here?", and the only answer so far was to parse the
    config's YAML. This is a plain table for someone with grep and no refgenie:
    each genome, the names it goes by, where and when it came from, the checksums
    of what was delivered, and every asset and tag it holds.

    One line per download record; a genome ingested before records were kept in
    the vault gets a single line with those columns empty.

    Rewritten in full on every run from the vault itself, never appended to. An
    appended log would drift from the vault the first time anything was removed
    or rebuilt by hand, and then it would be a second, wrong, source of truth.
    """
    with open(config_file, "r") as f:
        config = yaml.safe_load(f) or {}

    rows = []
    for digest, data in sorted((config.get("genomes") or {}).items()):
        data = data or {}
        assets = []
        for asset in sorted(data.get("assets") or {}):
            tags = sorted(asset_tags_on_disk(config, config_file, digest, asset))
            if tags:
                assets.append(f"{asset}:{','.join(tags)}")
        base = {
            "digest": digest,
            "aliases": ";".join(sorted(data.get("aliases") or [])),
            "assets": ";".join(assets),
        }
        records = []
        for tag in sorted(asset_tags_on_disk(config, config_file, digest, "provenance")):
            path = asset_file(config, config_file, digest, "provenance", tag, "provenance")
            try:
                with open(path) as fh:
                    records.append(json.load(fh))
            except (OSError, ValueError):
                continue
        if not records:
            rows.append(base)
        for record in records:
            rows.append(dict(
                base,
                provider=record.get("provider", ""),
                species=record.get("species", ""),
                assembly=record.get("assembly", ""),
                annotation_release=record.get("annotation_release", ""),
                assembly_accession=record.get("assembly_accession", ""),
                mask=(record.get("install_arguments") or {}).get("mask", ""),
                retrieved_utc=record.get("retrieved_utc", ""),
                genome_content_sha256=(record.get("genome") or {}).get("content_sha256", ""),
                annotation_content_sha256=(record.get("annotation") or {}).get(
                    "content_sha256", ""),
            ))

    directory = os.path.join(outdir, "provenance")
    os.makedirs(directory, exist_ok=True)
    # Written beside its destination and moved into place, so a reader never sees
    # half a table.
    fd, tmp_path = tempfile.mkstemp(dir=directory, prefix=".index-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=INDEX_FIELDS, delimiter="\t",
                                    restval="")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_path, os.path.join(directory, "index.tsv"))
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    logger.info("Wrote the vault index (%d line(s)): %s", len(rows),
                os.path.join(directory, "index.tsv"))


def process_assembly(row, outdir, config_file, refgenie_bin, do_cleanup, logger,
                     sources=None):
    """
    Register one resolved request into the refgenie vault.

    `sources` maps source_key() to a staged download directory. A row found there
    is ingested from that directory; any other row is looked for in the published
    layout under outdir, as before.

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
    query      = (row.get("query")      or "").strip()
    # The release the request PINNED, if any - distinct from the annotation release
    # read back from provenance.json further down, which exists for unpinned
    # requests too. This one only decides which directory the download is in.
    release_pin = (row.get("release")   or "").strip()
    # Carried on the resolved row by bin/resolve.py. They used to be looked up in
    # the request CSV by (species, query), a key two providers' copies of one
    # assembly share, so each received the other's aliases.
    requested_aliases = parse_aliases(row.get("aliases"))

    logger.info("Processing %s...", assembly)

    # -------------------------------------------------------
    # Path traversal guard
    # Validate that the assembly directory resolves inside outdir
    # before any filesystem operation (glob, build, rmtree).
    # -------------------------------------------------------
    #
    # A pinned release sits one level down, in release-<n>/, matching the
    # publishDir in main.nf. The two build this path independently; keep them in
    # step.
    release_parts = (f"release-{release_pin}",) if release_pin else ()
    try:
        layout_dir = safe_path_under(outdir, provider, species, *release_parts, assembly)
        vault_name = vault_genome_name(assembly, provider)
    except ValueError as e:
        logger.error("Skipping %s - %s", assembly, e)
        return [f"{assembly}: {e}"]

    # Where to read the download from. layout_dir is where it is PUBLISHED, and is
    # the only directory cleanup may ever remove; a staged source lives in
    # Nextflow's work directory and is not this script's to delete.
    staged_dir = (sources or {}).get(source_key(provider, assembly, release_pin))
    asm_dir = staged_dir or layout_dir
    if staged_dir:
        # A source that was named must be usable. Falling through to "nothing to
        # ingest" is right for a row that was cleaned up long ago, and wrong here.
        mismatch = (
            source_mismatch(staged_dir, provider, assembly, release_pin)
            if os.path.isdir(staged_dir)
            else f"the source directory '{staged_dir}' does not exist"
        )
        if mismatch:
            logger.error("Skipping %s - %s", assembly, mismatch)
            return [f"{assembly}: {mismatch}"]

    # Names offered as aliases without anyone asking for them by name. The first
    # genome to arrive keeps each; see merge_aliases.
    auto_names = {assembly, assembly.replace(".", "_"), assembly.split(".")[0]}
    if query:
        auto_names.add(query)

    # Track what actually made it into the vault. Cleanup is gated on these:
    # source data must never be removed on the strength of an ingest that only
    # partially succeeded.
    fasta_registered = False
    gtf_registered = False
    # Aliases are normally merged mid-flow, right after the FASTA build. This records
    # whether that happened so the fallback at the end does not repeat it.
    aliases_merged = False
    fasta_path = None
    # Filled in just before the FASTA build; see snapshot_aliases.
    aliases_before = None
    # True when the vault's FASTA was replaced from this download, which is the one
    # time the genome's existing download record should be replaced with it.
    fasta_rebuilt = False

    if os.path.exists(asm_dir):
        # Find and register FASTA
        fasta_files = (
            glob.glob(os.path.join(asm_dir, "*.fa")) +
            glob.glob(os.path.join(asm_dir, "*.fa.gz"))
        )
        if fasta_files:
            fasta_path = os.path.abspath(fasta_files[0])
            # Before the build, which resets this genome's alias list.
            aliases_before = snapshot_aliases(config_file)
            logger.info("Building FASTA digest in refgenie for '%s'", vault_name)
            cmd = [
                sys.executable, refgenie_bin, "build",
                f"{vault_name}/fasta",
                "--files", f"fasta={fasta_path}",
                "--recipe", os.path.abspath(FASTA_RECIPE),
                "-c", config_file
            ]
            try:
                subprocess.run(cmd, check=True)
                fasta_registered = True
            except subprocess.CalledProcessError as e:
                logger.error("FASTA build failed for '%s': %s", assembly, e)
                problems.append(f"{assembly}: FASTA build failed ({e})")

            # A genome that was already in the vault was not rebuilt by the call
            # above, so one ingested before the dictionary was part of the recipe
            # still has none. Checked AFTER the build rather than before it,
            # because only now is this genome sure to be findable by its vault
            # name. -N (pypiper's new-start) reruns the whole recipe from the
            # source FASTA, which is here; the sequence and its digest are the
            # same, and the alias merge that follows restores the alias list.
            if fasta_registered and fasta_has_dictionary(config_file, vault_name) is False:
                logger.info(
                    "'%s' is in the vault without a sequence dictionary (it was "
                    "ingested before one was built). Rebuilding its fasta asset "
                    "to add it.", vault_name
                )
                try:
                    subprocess.run(cmd + ["-N"], check=True)
                    fasta_rebuilt = True
                except subprocess.CalledProcessError as e:
                    logger.error("FASTA rebuild failed for '%s': %s", assembly, e)
                    problems.append(f"{assembly}: FASTA rebuild for the dictionary failed ({e})")
                    fasta_registered = False
                if fasta_registered and not fasta_has_dictionary(config_file, vault_name):
                    logger.error(
                        "'%s' still has no sequence dictionary after rebuilding", vault_name
                    )
                    problems.append(f"{assembly}: sequence dictionary was not built")
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
                vault_name, assembly, auto_names, requested_aliases,
                config_file, logger, aliases_before
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
                target = f"{vault_name}/{recipe_name}"
                if release:
                    target += f":{release}"
                # Read BEFORE the build: once the new tag exists there is no telling
                # whether refgenie made it the default (it does for an asset's first
                # tag) or left an older one in place.
                previous_default = current_default_tag(config_file, vault_name, recipe_name)
                logger.info(
                    "Building GTF annotation digest in refgenie for '%s' (provider: %s, "
                    "recipe: %s%s, tag: %s%s)", vault_name, provider, recipe_name,
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
                # is the opposite of the intended "newest wins, older stays reachable"
                # behaviour. Not fatal: the asset is built and addressable by tag.
                #
                # "Newest", not "latest ingested": an older release pinned after a
                # newer one is already in the vault must not take the default from
                # it. See should_become_default.
                if gtf_registered and release and not should_become_default(
                        previous_default, release):
                    logger.info(
                        "Built '%s/%s:%s'. The default tag stays at '%s', which is a "
                        "newer release; name the tag to use this one.",
                        vault_name, recipe_name, release, previous_default
                    )
                elif gtf_registered and release:
                    try:
                        subprocess.run(
                            [sys.executable, refgenie_bin, "tag",
                             f"{vault_name}/{recipe_name}:{release}",
                             "--default", "-c", config_file],
                            check=True
                        )
                    except subprocess.CalledProcessError as e:
                        logger.warning(
                            "Built '%s/%s:%s' but could not make it the default tag: %s. "
                            "The asset is addressable by tag; callers using no tag will "
                            "still resolve to the previous release.",
                            vault_name, recipe_name, release, e
                        )
            else:
                # Annotation was explicitly requested, so its absence is a failure of
                # the run, not a cosmetic warning.
                logger.error(
                    "Annotation requested for '%s' but no GTF found in '%s'",
                    assembly, asm_dir
                )
                problems.append(f"{assembly}: annotation requested but no GTF found")

        # The download's own record, into the vault. After the aliases are merged,
        # so it is built against the full alias set, and before the relink below,
        # so that covers it too.
        if fasta_registered:
            problems.extend(register_download_record(
                asm_dir, vault_name, assembly, refgenie_bin, config_file, logger,
                force=fasta_rebuilt
            ))

        # Every asset for this genome now exists and the aliases are merged, so
        # refresh the symlink tree once, covering assets built before the aliases
        # were attached (the FASTA in particular).
        if aliases_merged:
            relink_aliases(config_file, vault_name, logger)

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
            preserve_provenance(asm_dir, outdir, provider, species, assembly, logger,
                                release_parts)
            # layout_dir, not asm_dir. They are the same directory when this script
            # is run by hand. Under Nextflow asm_dir is a staged input, and the
            # published copy exists only if an earlier run without --cleanup left
            # one behind - which is exactly the copy --cleanup is asked to remove.
            if os.path.exists(layout_dir):
                logger.info("Cleaning up download directory: '%s'", layout_dir)
                if not logged_rmtree(layout_dir, logger):
                    logger.warning(
                        "Directory '%s' may require manual cleanup. See '%s' for full "
                        "details.", layout_dir, os.path.join(outdir, "update_refgenie.log")
                    )
                    problems.append(f"{assembly}: cleanup did not complete cleanly")
    else:
        # This is the normal path on a re-run after a previous --cleanup, but it is
        # also what a provider/species case mismatch looks like, because the download
        # and ingest stages construct this path independently. The message no longer
        # asserts that the genome is present.
        # A genome the vault already holds is the ordinary case since the preflight:
        # nothing was downloaded for it because nothing needed to be. Say that, and
        # keep the longer hint for the case where the genome is NOT there, which is
        # the one that needs explaining.
        if fasta_has_dictionary(config_file, vault_name):
            logger.info(
                "'%s' is already in the vault and no download was supplied for it; "
                "checking its aliases only.", vault_name
            )
        else:
            logger.info(
                "Directory '%s' not found - nothing to ingest for '%s'. If this assembly "
                "was not registered by an earlier run, check that the provider and "
                "species in the resolved CSV match the published directory layout.",
                asm_dir, assembly
            )
        # The dictionary is built from the source FASTA, so with no source there is
        # nothing to build it from. Say so: otherwise `refgenie seek <genome>/fasta.dict`
        # fails later with no hint that the fix is a re-download.
        if fasta_has_dictionary(config_file, vault_name) is False:
            logger.warning(
                "'%s' is in the vault without a sequence dictionary, and its download "
                "is not here to build one from. It will be added the next time this "
                "genome is ingested from a download.", vault_name
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
            vault_name, assembly, auto_names, requested_aliases,
            config_file, logger, aliases_before
        ))

    return problems


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Register downloaded reference data into a refgenie vault."
    )
    parser.add_argument("resolved_csv", help="The resolved request CSV produced by bin/resolve.py")
    parser.add_argument("outdir", help="Output directory holding the refgenie vault")
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove each downloaded source directory after its contents have been "
             "verified into the vault. Off by default: upstream providers archive old "
             "releases, so a deleted download may not be re-fetchable.",
    )
    parser.add_argument(
        "--source", action="append", default=[], metavar="PROVIDER/ASSEMBLY@[RELEASE]=DIR",
        help="Ingest this request from DIR instead of the published layout under "
             "outdir. Repeatable. Used by the REGISTER process in main.nf, which is "
             "handed each download as a staged input.",
    )
    parser.add_argument(
        "--manifest", metavar="TSV",
        help="Write the requested genomes that are in the vault to this file, with "
             "their vault names and digests.",
    )
    args = parser.parse_args()

    resolved_csv = args.resolved_csv
    outdir       = args.outdir

    try:
        sources = parse_sources(args.source)
    except ValueError as e:
        parser.error(str(e))

    # Ensure outdir exists before we try to create the log file inside it
    os.makedirs(outdir, exist_ok=True)

    logger      = setup_logger(outdir)
    config_file = os.path.join(outdir, "refgenie_export.yaml")
    refgenie_bin = os.path.join(os.path.dirname(sys.executable), "refgenie")

    logger.info("update_refgenie starting | outdir='%s' | cleanup=%s", outdir, args.cleanup)

    # Every genome is built with this recipe, so without it nothing can be ingested.
    # Checked once here: otherwise each row fails separately with refgenie's own,
    # less direct, complaint about the path.
    if not os.path.isfile(FASTA_RECIPE):
        logger.error("The fasta recipe '%s' is missing; nothing can be ingested.",
                     os.path.abspath(FASTA_RECIPE))
        sys.exit(1)

    # Ensure refgenie config is initialized
    if not os.path.exists(config_file):
        logger.info("Initializing refgenie config at '%s'", config_file)
        subprocess.run(
            [sys.executable, refgenie_bin, "init", "-c", config_file],
            check=True
        )

    # Process each resolved request. Each row is isolated so a failure part-way
    # through the batch no longer abandons the assemblies after it.
    all_problems = []
    processed = 0
    rows = []
    used_sources = set()
    try:
        with open(resolved_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            for lineno, row in enumerate(reader, start=2):
                assembly = (row.get("assembly") or "").strip() or f"row {lineno}"
                processed += 1
                rows.append(row)
                used_sources.add(source_key(
                    row.get("provider") or "", (row.get("assembly") or "").strip(),
                    (row.get("release") or "").strip()
                ))
                try:
                    all_problems.extend(
                        process_assembly(
                            row, outdir, config_file, refgenie_bin,
                            args.cleanup, logger, sources
                        )
                    )
                except Exception as e:
                    logger.error("Unhandled error processing '%s': %s", assembly, e)
                    all_problems.append(f"{assembly}: unhandled error ({e})")
    except Exception as e:
        logger.error("Fatal error reading '%s': %s", resolved_csv, e)
        sys.exit(1)

    # A source that matched no row was downloaded and then ignored. That is a
    # disagreement between the workflow and the request file, not something to
    # pass over: the genome it holds is not in the vault.
    for key in sorted(set(sources) - used_sources):
        logger.error("--source '%s' matches no row in '%s'", key, resolved_csv)
        all_problems.append(f"source '{key}' matches no request and was not ingested")

    # The plain-text index of the vault, rebuilt from the vault itself.
    try:
        write_vault_index(outdir, config_file, logger)
    except Exception as e:
        logger.error("Could not write the vault index: %s", e)
        all_problems.append(f"vault index not written ({e})")

    # Written whether or not there were problems: it lists what IS in the vault,
    # and the exit code below says whether that is everything that was asked for.
    if args.manifest:
        try:
            write_manifest(args.manifest, rows, config_file)
        except Exception as e:
            logger.error("Could not write the manifest '%s': %s", args.manifest, e)
            all_problems.append(f"manifest not written ({e})")

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
