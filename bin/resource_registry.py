#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
The registry of third-party resources: what can be requested, and from where.

A genome and its annotation come from a provider that genomepy can search. A
resource is different in kind: a dataset somebody else curates and versions on
their own schedule - ClinVar, the ENCODE blacklist - that merely happens to be
written in an assembly's coordinates. Nothing derives it and no provider search
finds it. So each one is described by a file in resources/, and the description
is the whole of what the pipeline knows about it:

    name         also the refgenie asset it is registered as, and the name of its
                 recipe in recipes/
    kind         'vcf' or 'bed': how its contigs and coordinates are read
    index        'tabix', or null to store it as fetched
    terms        where its conditions of use are stated
    min_matched  the share of its records that must fall on contigs the genome
                 has before it may be registered (see bin/prepare_resource.py)
    variants     one entry per genome build it is published for, each with a url,
                 a version, and optionally a published checksum

A request names a resource and a variant together - `clinvar:GRCh38` - because
which build a file is for cannot be worked out from the genome. Two builds of one
species can share every contig name; only the person asking knows which they
hold. The check that follows the download is a safeguard, not a substitute.

`version` is either a fixed label ("v2"), for a file that is published once, or
"last_modified", for one that is replaced in place: its version is then the date
the publisher last changed it, the same rule NCBI and UCSC annotations are tagged
by (bin/download.py).

JSON, not YAML, on purpose: this module is used by the fetch step, which runs
with whatever Python the profile supplies and nothing but the standard library.
"""
import json
import os
import re

REGISTRY_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "resources")

# A resource name becomes a refgenie asset name, a file name and part of shell
# commands; a variant becomes part of shell commands.
NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
KINDS = ("vcf", "bed")


def load_entry(path):
    """One registry file, checked. Raises ValueError naming the file and the fault."""
    try:
        with open(path) as fh:
            entry = json.load(fh)
    except (OSError, ValueError) as e:
        raise ValueError(f"cannot read resource definition '{path}': {e}")

    def fault(message):
        return ValueError(f"resource definition '{path}': {message}")

    name = entry.get("name")
    if not isinstance(name, str) or not NAME_RE.match(name):
        raise fault("'name' must be letters, digits, '.', '_' and '-'")
    if os.path.splitext(os.path.basename(path))[0] != name:
        raise fault(f"is for '{name}', so the file must be called '{name}.json'")
    if entry.get("kind") not in KINDS:
        raise fault(f"'kind' must be one of {', '.join(KINDS)}")
    if entry.get("index") not in (None, "tabix"):
        raise fault("'index' must be 'tabix' or null")
    variants = entry.get("variants")
    if not isinstance(variants, dict) or not variants:
        raise fault("needs at least one entry under 'variants'")
    for label, variant in variants.items():
        if not NAME_RE.match(label):
            raise fault(f"variant '{label}' must be letters, digits, '.', '_' and '-'")
        if not str(variant.get("url", "")).startswith(("https://", "http://")):
            raise fault(f"variant '{label}' needs an http(s) 'url'")
        if not variant.get("version"):
            raise fault(f"variant '{label}' needs a 'version' (a label, or 'last_modified')")
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
