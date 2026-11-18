#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
The registries of third-party data: what can be requested, and from where.

A genome and its annotation come from a provider that genomepy can search. What
is described here is different in kind: a dataset somebody else curates and
versions on their own schedule, which no provider search finds. There are two
sorts, each described by one file per dataset, and the description is the whole
of what the pipeline knows about it.

RESOURCES (resources/<name>.json) are written in an assembly's coordinates -
ClinVar, the ENCODE blacklist, the GATK known-sites files, dbSNP, a VEP cache -
and are registered against a genome in the vault:

    name         also the refgenie asset it is registered as, and the name of its
                 recipe in recipes/
    kind         'vcf' or 'bed': how its contigs and coordinates are read; or
                 'vep_cache', an Ensembl VEP cache tarball, whose sequence
                 directories are compared with the genome instead
    index        'tabix', or null to store it as fetched
    terms        where its conditions of use are stated
    min_matched  the share that must fit the genome before it may be registered
                 (see bin/prepare_resource.py)
    variants     one entry per genome build it is published for

A request names a resource and a variant together - `clinvar:GRCh38` - because
which build a file is for cannot be worked out from the genome. Two builds of one
species can share every contig name; only the person asking knows which they
hold. The check that follows the download is a safeguard, not a substitute.

DATABASES (databases/<name>.json) belong to no genome - Swiss-Prot, Pfam, SILVA,
a Kraken2 index. refgenie cannot hold them, since every asset there hangs off a
genome, so they go in a tree of their own beside the vault (bin/file_database.py).
A database's file has the fields of a resource's variant at its top level, plus
`unpack`: 'tar' for an archive to be unpacked, or null.

WHAT A VARIANT (OR A DATABASE) SAYS

    url        where the file is. For a variant whose version is
               'annotation_release' it contains {release}.
    version    how the file is versioned, one of:
                 a fixed label ("v2")    for a file published once;
                 "last_modified"         for one replaced in place: its version
                                         is the date the publisher last changed
                                         it, the same rule NCBI and UCSC
                                         annotations are tagged by;
                 {"url", "pattern"}      for one whose publisher states the
                                         release in a small file beside it: the
                                         pattern's first group is the version;
                 "annotation_release"    (resources only) for one published per
                                         release of the genome's own annotation,
                                         such as a VEP cache: the version is the
                                         release the genome's annotation is at.
    annotation which annotation asset 'annotation_release' reads the release
               from (default ensembl_gtf).
    checksum   {"type": "md5"|"sha256", "url": ...} for a checksum posted in a
               file, or {"type": "md5", "header": "x-goog-hash"} for one sent
               with the download itself. Optional.
    rename     (vcf and bed only) {"report": address of an NCBI assembly report}:
               the file's sequence names are rewritten to the genome's own
               through that table. For a file such as dbSNP, which names its
               sequences by RefSeq accession, as no genome in a vault does.
    size       roughly how big it is, in words, said before it is fetched.
    naming     what naming it uses and what it fits, for the refusal message.

JSON, not YAML, on purpose: this module is used by the fetch steps, which run
with whatever Python the profile supplies and nothing but the standard library.
"""
import json
import os
import re

REGISTRY_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "resources")
DATABASES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "databases")

# A resource name becomes a refgenie asset name, a file name and part of shell
# commands; a variant becomes part of shell commands. A database name becomes a
# directory name.
NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
KINDS = ("vcf", "bed", "vep_cache")
CHECKSUM_TYPES = ("md5", "sha256")
UNPACK = (None, "tar")


def is_address(value):
    return str(value or "").startswith(("https://", "http://"))


def source_fault(source, for_genome):
    """
    What is wrong with one place a file is fetched from, or None.

    `source` is a resource's variant or a database's whole definition: the two
    describe a download in the same words. `for_genome` says which it is, since
    only a resource has a genome whose annotation release it can follow.
    """
    if not is_address(source.get("url")):
        return "needs an http(s) 'url'"
    version = source.get("version")
    if isinstance(version, dict):
        if not is_address(version.get("url")) or not version.get("pattern"):
            return "a 'version' given as an object needs a 'url' and a 'pattern'"
        try:
            if re.compile(version["pattern"]).groups < 1:
                return "the 'version' pattern needs one group, which is the version"
        except re.error as e:
            return f"the 'version' pattern is not a regular expression: {e}"
    elif not version:
        return "needs a 'version' (a label, 'last_modified', or a url and pattern)"
    follows_annotation = version == "annotation_release"
    if follows_annotation and not for_genome:
        return "'annotation_release' is only for a resource, which has a genome"
    if ("{release}" in str(source["url"])) != follows_annotation:
        return ("'{release}' in the url and the version 'annotation_release' go "
                "together: one was given without the other")
    checksum = source.get("checksum")
    if checksum is not None:
        if not isinstance(checksum, dict) or checksum.get("type") not in CHECKSUM_TYPES:
            return f"a 'checksum' needs a 'type', one of {', '.join(CHECKSUM_TYPES)}"
        if bool(checksum.get("url")) == bool(checksum.get("header")):
            return "a 'checksum' needs either a 'url' or a 'header'"
        if checksum.get("url") and not is_address(checksum["url"]):
            return "a checksum 'url' must be an http(s) address"
    rename = source.get("rename")
    if rename is not None:
        if not for_genome:
            return "'rename' is only for a resource"
        if not isinstance(rename, dict) or not is_address(rename.get("report")):
            return "'rename' needs a 'report': the address of an NCBI assembly report"
    return None


def read_definition(path, what):
    try:
        with open(path) as fh:
            entry = json.load(fh)
    except (OSError, ValueError) as e:
        raise ValueError(f"cannot read {what} definition '{path}': {e}")
    if not isinstance(entry, dict):
        raise ValueError(f"{what} definition '{path}': must be a JSON object")
    name = entry.get("name")
    if not isinstance(name, str) or not NAME_RE.match(name):
        raise ValueError(f"{what} definition '{path}': 'name' must be letters, digits, "
                         "'.', '_' and '-'")
    if os.path.splitext(os.path.basename(path))[0] != name:
        raise ValueError(f"{what} definition '{path}': is for '{name}', so the file must "
                         f"be called '{name}.json'")
    return entry


def load_entry(path):
    """One resource definition, checked. Raises ValueError naming the file and the fault."""
    entry = read_definition(path, "resource")

    def fault(message):
        return ValueError(f"resource definition '{path}': {message}")

    if entry.get("kind") not in KINDS:
        raise fault(f"'kind' must be one of {', '.join(KINDS)}")
    if entry.get("index") not in (None, "tabix"):
        raise fault("'index' must be 'tabix' or null")
    if entry["kind"] == "vep_cache" and entry.get("index"):
        raise fault("a 'vep_cache' is a directory and cannot have an 'index'")
    variants = entry.get("variants")
    if not isinstance(variants, dict) or not variants:
        raise fault("needs at least one entry under 'variants'")
    for label, variant in variants.items():
        if not NAME_RE.match(label):
            raise fault(f"variant '{label}' must be letters, digits, '.', '_' and '-'")
        problem = source_fault(variant, for_genome=True)
        if problem:
            raise fault(f"variant '{label}' {problem}")
        if variant.get("rename") and entry["kind"] == "vep_cache":
            raise fault(f"variant '{label}': a 'vep_cache' cannot be renamed")
    return entry


def load_registry(directory=REGISTRY_DIR):
    """Every resource definition in a directory, as {name: entry}."""
    registry = {}
    if not os.path.isdir(directory):
        return registry
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".json"):
            entry = load_entry(os.path.join(directory, filename))
            registry[entry["name"]] = entry
    return registry


def load_database(path):
    """One database definition, checked. Raises ValueError naming the file and the fault."""
    entry = read_definition(path, "database")
    problem = source_fault(entry, for_genome=False)
    if problem:
        raise ValueError(f"database definition '{path}': {problem}")
    if entry.get("unpack") not in UNPACK:
        raise ValueError(f"database definition '{path}': 'unpack' must be 'tar' or null")
    return entry


def load_databases(directory=DATABASES_DIR):
    """Every database definition in a directory, as {name: entry}."""
    databases = {}
    if not os.path.isdir(directory):
        return databases
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".json"):
            entry = load_database(os.path.join(directory, filename))
            databases[entry["name"]] = entry
    return databases


def parse_request(value, registry):
    """
    'clinvar:GRCh38; blacklist:hg38' -> ([('clinvar', 'GRCh38'), ...], None)

    or ([], reason) when the request cannot be met. The reason names what IS
    available, since a mistyped name is the usual cause.
    """
    pairs = []
    for item in str(value or "").split(";"):
        item = item.strip()
        if not item:
            continue
        name, separator, variant = item.partition(":")
        name, variant = name.strip(), variant.strip()
        if not separator or not name or not variant:
            return [], (f"resource '{item}' must be written as name:build, for example "
                        "clinvar:GRCh38")
        if name not in registry:
            known = ", ".join(sorted(registry)) or "none are defined"
            return [], f"unknown resource '{name}' (known: {known})"
        if variant not in registry[name]["variants"]:
            known = ", ".join(sorted(registry[name]["variants"]))
            return [], (f"resource '{name}' is not published for '{variant}' "
                        f"(it has: {known})")
        if (name, variant) in pairs:
            return [], f"resource '{item}' is listed twice"
        if any(other == name for other, _ in pairs):
            return [], (f"resource '{name}' is requested for two builds on one row; a "
                        "genome is one build")
        pairs.append((name, variant))
    return pairs, None


def format_request(pairs):
    """[('clinvar', 'GRCh38')] -> 'clinvar:GRCh38'"""
    return ";".join(f"{name}:{variant}" for name, variant in pairs)
