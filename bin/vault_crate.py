#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Describe the vault as an RO-Crate: <outdir>/ro-crate-metadata.json.

The vault already records where everything came from, in a provenance.json beside
each asset. Those files are written for this pipeline and say exactly what it
knows; nothing else can read them without being told how. RO-Crate is the same
information in a vocabulary other tools already speak - schema.org, as JSON-LD -
so that a repository, a FAIR assessment or a colleague's script can find out what
the vault holds, how big each file is, and what made it, without knowing anything
about refdata-flow.

This does not replace the provenance files. They stay where they are, they are
listed in the crate as files that describe their asset, and they remain the only
place for what has no standard home: the summary of an annotation, the result of
checking a resource against a genome.

WHAT IS IN THE CRATE

    the vault                 the root dataset
      each genome             a Dataset, data/<digest>/, with its aliases
        each asset and tag    a Dataset, data/<digest>/<asset>/<tag>/
          each file           a File, with its size and format

    each download, each       a CreateAction: what ran (and which version), what
    derivation, each fetch    it read, what it produced, and when it finished

A derivation is recognised by its record, not by a list of asset names kept here:
any asset whose provenance file names itself and a tool. What it was derived from
is refgenie's own record of the asset's parents. So a new derived asset appears
in the crate with no change to this file.

WHICH STANDARD, AND WHY THESE VERSIONS

RO-Crate 1.1, with the Process Run Crate profile 0.5 for the actions. Both have
newer editions (1.3 and 0.6, checked 2026-11-13). These are the ones the reference
validator (roc-validator 0.13.1) can check a crate against, and a claim of
conformance that nothing can check is not worth making. The two version strings
are constants below.

HOW IT IS KEPT TRUE

Like the index, it is rebuilt from the vault itself after every registration and
never edited. Three rules follow from that:

  * A file is listed only if it is on disk now.
  * A checksum is given only where the pipeline recorded one that provably applies
    to the stored bytes. Hashing a three-gigabyte genome on every run to fill in a
    field would be the wrong trade; an absent `sha256` means "not recorded", and
    the size is always there.
  * Nothing in it depends on when it was written. `datePublished` is the latest
    retrieval or derivation among the vault's records, so two runs over an
    unchanged vault produce the same bytes.
"""
import json
import os
import tempfile

import yaml

import update_refgenie as vault

RO_CRATE_VERSION = "1.1"
PROCESS_RUN_VERSION = "0.5"

# Two contexts: RO-Crate's own, and the run profiles' terms, which is where
# `sha256` is defined - RO-Crate 1.1 has no word for a checksum (checked 2026-11-14).
CONTEXT = [f"https://w3id.org/ro/crate/{RO_CRATE_VERSION}/context",
           "https://w3id.org/ro/terms/workflow-run/context"]
RO_CRATE_SPEC = f"https://w3id.org/ro/crate/{RO_CRATE_VERSION}"
PROCESS_RUN_PROFILE = f"https://w3id.org/ro/wfrun/process/{PROCESS_RUN_VERSION}"
# Written as a link, as every example in the RO-Crate specification writes it.
# roc-validator 0.13.1 compares it with a text value at its "recommended" level
# and so reports it there; the specification's own form is kept.
COMPLETED = {"@id": "http://schema.org/CompletedActionStatus"}

PIPELINE_URL = "https://github.com/nicholas-owen/refdata-flow"
KENT_URL = "https://github.com/ucscGenomeBrowser/kent"
SOFTWARE_URLS = {
    "genomepy": "https://github.com/vanheeringen-lab/genomepy",
    "gffread": "https://github.com/gpertea/gffread",
    "htslib": "https://github.com/samtools/htslib",
    "gtfToGenePred": KENT_URL,
    "genePredToBed": KENT_URL,
    "faToTwoBit": KENT_URL,
    "refdata-flow": PIPELINE_URL,
}

# What a derivation did, where one sentence says more than the generic one.
DERIVATIONS = {
    "fasta_txome": "{tool} cut every transcript in the annotation out of the genome; "
                   "sequence was then written in upper case.",
}

# EDAM format identifiers, checked against EDAM.tsv on 2026-11-13 and -14.
EDAM = {
    "FASTA": "http://edamontology.org/format_1929",
    "GTF": "http://edamontology.org/format_2306",
    "VCF": "http://edamontology.org/format_3016",
    "BED": "http://edamontology.org/format_3003",
    "tabix": "http://edamontology.org/format_3616",
    "2bit": "http://edamontology.org/format_3009",
    "TSV": "http://edamontology.org/format_3475",
}

# (suffix, media types, EDAM name or None). First match wins, so longer first.
FORMATS = [
    (".gz.tbi", [], "tabix"),
    (".fa.fai", ["text/tab-separated-values"], None),
    (".chrom.sizes", ["text/tab-separated-values"], None),
    (".vcf.gz", ["application/gzip"], "VCF"),
    (".gtf.gz", ["application/gzip"], "GTF"),
    (".bed.gz", ["application/gzip"], "BED"),
    (".bed", ["text/plain"], "BED"),
    (".fa", ["text/plain"], "FASTA"),
    (".2bit", [], "2bit"),
    (".tsv", ["text/tab-separated-values"], "TSV"),
    (".refFlat.txt", ["text/tab-separated-values"], None),
    (".json", ["application/json"], None),
    (".dict", ["text/plain"], None),
    (".txt", ["text/plain"], None),
]


def ref(identifier):
    return {"@id": identifier}


class Crate:
    """The graph under construction: entities keyed by @id, so none is added twice."""

    def __init__(self):
        self.entities = {}

    def add(self, entity):
        self.entities.setdefault(entity["@id"], entity)
        return ref(entity["@id"])

    def has(self, identifier):
        return identifier in self.entities

    def software(self, name, version):
        """A tool at a version, as something an action can name as its instrument."""
        version = str(version or "").strip()
        entity = {
            "@id": f"#{name}-{version}" if version else f"#{name}",
            "@type": "SoftwareApplication",
            "name": name,
            "url": SOFTWARE_URLS.get(name, PIPELINE_URL),
        }
        if version:
            entity["version"] = version
        return self.add(entity)

    def source(self, address, description):
        """
        Something that was read from outside the vault: a publisher's file, or for
        UCSC annotations a database table, which has no address of its own.
        """
        address = str(address or "").strip()
        if not address:
            return None
        if address.startswith(("http://", "https://", "ftp://")):
            return self.add({"@id": address, "@type": "CreativeWork",
                             "name": address.rsplit("/", 1)[-1] or address,
                             "description": description, "url": address})
        local = "#source-" + "".join(c if c.isalnum() else "-" for c in address).strip("-")
        return self.add({"@id": local, "@type": "CreativeWork", "name": address,
                         "description": description})


def encoding_format(crate, filename):
    for suffix, media_types, edam in FORMATS:
        if filename.endswith(suffix):
            formats = list(media_types)
            if edam:
                formats.append(crate.add({"@id": EDAM[edam], "@type": "WebSite", "name": edam}))
            return formats[0] if len(formats) == 1 else formats
    return None


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def tag_directory(digest, asset, tag):
    return f"data/{digest}/{asset}/{tag}/"


def build_crate(outdir, config_file):
    """The whole crate as a dict, built from the vault's config and what is on disk."""
    with open(config_file, "r") as f:
        config = yaml.safe_load(f) or {}
    genomes = config.get("genomes") or {}
    # Identifiers in a crate are paths relative to the crate's own directory, and
    # the pipeline always keeps the vault's data/ tree beside its config there.
    root = os.path.dirname(os.path.abspath(config_file))

    crate = Crate()
    actions, dates, genome_refs = [], [], []
    # {(digest, asset, tag): {seek key: File @id}} for files that are on disk.
    files_of = {}

    # ---- pass 1: genomes, assets, files ----------------------------------------
    for digest in sorted(genomes):
        data = genomes[digest] or {}
        aliases = sorted(data.get("aliases") or [])
        asset_refs = []
        for asset in sorted(data.get("assets") or {}):
            asset_data = data["assets"][asset] or {}
            for tag in sorted(vault.asset_tags_on_disk(config, config_file, digest, asset)):
                tag_data = next((v for k, v in (asset_data.get("tags") or {}).items()
                                 if str(k) == tag), None) or {}
                directory = tag_directory(digest, asset, tag)
                parts, held = [], {}
                for key, filename in sorted((tag_data.get("seek_keys") or {}).items()):
                    path = vault.asset_file(config, config_file, digest, asset, tag, key)
                    if not path:
                        continue
                    entity = {
                        "@id": directory + str(filename),
                        "@type": "File",
                        "name": str(filename),
                        "contentSize": str(os.path.getsize(path)),
                    }
                    formats = encoding_format(crate, str(filename))
                    if formats:
                        entity["encodingFormat"] = formats
                    held[key] = entity["@id"]
                    parts.append(crate.add(entity))
                if not parts:
                    continue
                files_of[(digest, asset, tag)] = held
                asset_refs.append(crate.add({
                    "@id": directory,
                    "@type": "Dataset",
                    "name": f"{asset}:{tag}",
                    "description": str(asset_data.get("asset_description") or
                                       f"The {asset} asset at tag {tag}."),
                    "version": tag,
                    "hasPart": parts,
                }))
        if asset_refs:
            genome_refs.append(crate.add({
                "@id": f"data/{digest}/",
                "@type": "Dataset",
                "name": aliases[0] if aliases else digest,
                "description": "A reference genome and the assets held for it. The "
                               "identifier is refgenie's digest of its sequences; the "
                               "alternate names are the aliases it can be asked for by.",
                "identifier": digest,
                "alternateName": aliases,
                "hasPart": asset_refs,
            }))

    def file_ref(digest, asset, tag, key):
        identifier = files_of.get((digest, asset, str(tag)), {}).get(key)
        return ref(identifier) if identifier else None

    # ---- pass 2: lineage between assets, and what made each one -----------------
    for digest in sorted(genomes):
        data = genomes[digest] or {}
        short = digest[:12]
        for (d, asset, tag), held in sorted(files_of.items()):
            if d != digest:
                continue
            directory = tag_directory(digest, asset, tag)
            dataset = crate.entities[directory]
            tag_data = next((v for k, v in ((data["assets"][asset] or {}).get("tags") or {}).items()
                             if str(k) == tag), None) or {}

            # refgenie's own record of what an asset was built from.
            parents, parent_keys = [], []
            for parent in tag_data.get("asset_parents") or []:
                parent_digest, _, rest = str(parent).partition("/")
                parent_asset, _, parent_tag = rest.partition(":")
                parent_dir = tag_directory(parent_digest, parent_asset, parent_tag)
                if crate.has(parent_dir):
                    parents.append(ref(parent_dir))
                    parent_keys.append((parent_digest, parent_asset, parent_tag))
            if parents:
                dataset["isBasedOn"] = parents

            record = load_json(os.path.join(root, held["provenance"])) if "provenance" in held else None
            if record is not None:
                crate.entities[held["provenance"]].update({
                    "description": "The pipeline's own record of where this asset came "
                                   "from, in full. The crate gives the standard part of it.",
                    "about": ref(directory),
                })
            if not record:
                continue

            if asset == "provenance":
                # A download: the genome, and its annotation if it had one.
                dataset["about"] = ref(f"data/{digest}/")
                if record.get("retrieved_utc"):
                    dates.append(record["retrieved_utc"])
                release = record.get("annotation_release")
                # Which annotation asset this download filled: the provider's, or
                # for a record too old to name its provider, whichever annotation
                # asset of this genome holds the release.
                named = vault.ANNOTATION_RECIPES.get(str(record.get("provider", "")).lower())
                recipe = next((r for r in [named] + sorted(set(vault.ANNOTATION_RECIPES.values()))
                               if r and (digest, r, str(release)) in files_of), None)
                genome_file = file_ref(digest, "fasta", "default", "fasta")
                gtf_file = file_ref(digest, recipe, release, recipe) if recipe and release else None
                if genome_file and (record.get("genome") or {}).get("content_sha256"):
                    # The vault keeps the FASTA decompressed, so the hash of the
                    # download's content is the hash of this file.
                    crate.entities[genome_file["@id"]].setdefault(
                        "sha256", record["genome"]["content_sha256"])
                if gtf_file and (record.get("annotation") or {}).get("file_sha256"):
                    # The annotation recipes copy the delivered file unchanged.
                    crate.entities[gtf_file["@id"]].setdefault(
                        "sha256", record["annotation"]["file_sha256"])
                what = f"{record.get('assembly', 'a genome')} from {record.get('provider', 'its provider')}"
                action = {
                    "@id": f"#download-{short}-{tag}",
                    "@type": "CreateAction",
                    "name": f"Download {what}",
                    "description": f"genomepy fetched {what}"
                                   + (f", with annotation release {release}" if release else "")
                                   + ", filtered and compressed it. See the record for the "
                                     "settings and for what genomepy changed.",
                    "instrument": crate.software(
                        "genomepy", (record.get("tool_versions") or {}).get("genomepy")),
                    "actionStatus": COMPLETED,
                }
                objects = [crate.source(record.get("genome_url"), "Genome sequence as published."),
                           crate.source(record.get("annotation_url"),
                                        "Gene annotation as published.") if release else None]
                results = [genome_file, gtf_file, ref(held["provenance"])]
            elif record.get("resource"):
                # A third-party resource fetched for this genome.
                if record.get("retrieved_utc"):
                    dates.append(record["retrieved_utc"])
                if record.get("terms"):
                    dataset["license"] = crate.add({
                        "@id": record["terms"], "@type": "CreativeWork",
                        "name": f"Terms of use for {record['resource']}",
                        "description": "Where the publisher states the conditions of use. "
                                       "Linked, not evaluated.",
                        "url": record["terms"]})
                data_file = file_ref(digest, asset, tag, asset)
                if data_file and record.get("stored_as_published") and record.get("file_sha256"):
                    crate.entities[data_file["@id"]].setdefault("sha256", record["file_sha256"])
                what = f"{record['resource']}:{record.get('variant', '')}"
                action = {
                    "@id": f"#fetch-{asset}-{short}-{tag}",
                    "@type": "CreateAction",
                    "name": f"Fetch {what}",
                    "description": f"refdata-flow fetched {what} at version "
                                   f"{record.get('version', tag)}, checked it against the "
                                   "genome's sequence names and lengths, and registered it. "
                                   "See the record for what the check found.",
                    "instrument": crate.software("refdata-flow", ""),
                    "actionStatus": COMPLETED,
                }
                objects = [crate.source(record.get("url"), record.get("description") or
                                        "The resource as published."),
                           file_ref(digest, "fasta", "default", "chrom_sizes")]
                results = [data_file, file_ref(digest, asset, tag, "tbi"), ref(held["provenance"])]
            elif record.get("asset") == asset and record.get("tool"):
                # Something derived from what the vault already held: a
                # transcriptome, a tx2gene table, a .2bit. What it was derived
                # from is refgenie's own record of the asset's parents, and the
                # file read from each parent is the one named after it.
                if record.get("derived_utc"):
                    dates.append(record["derived_utc"])
                tool = str(record["tool"])
                sources = ", ".join(f"{a}:{t}" for _, a, t in parent_keys) or "the vault"
                action = {
                    "@id": f"#derive-{asset}-{short}-{tag}",
                    "@type": "CreateAction",
                    "name": f"Derive {asset} from {sources}",
                    "description": DERIVATIONS.get(asset, "").format(tool=tool) or
                                   f"Derived from {sources} with {tool}. The record beside "
                                   "the files gives the command.",
                    "instrument": crate.software(tool, record.get("tool_version")),
                    "actionStatus": COMPLETED,
                }
                objects = [file_ref(d, a, t, a) for d, a, t in parent_keys]
                results = [file_ref(digest, asset, tag, asset), file_ref(digest, asset, tag, "tbi"),
                           ref(held["provenance"])]
            else:
                continue

            finished = record.get("retrieved_utc") or record.get("derived_utc")
            if finished:
                action["endTime"] = finished
            objects = [item for item in objects if item]
            if objects:
                action["object"] = objects
            action["result"] = [item for item in results if item]
            actions.append(crate.add(action))

    # ---- the root, and the two entities every crate has ------------------------
    licence = crate.add({
        "@id": "#terms",
        "@type": "CreativeWork",
        "name": "Terms vary by component",
        "description": "This vault gathers files from several publishers. Each remains "
                       "under its own publisher's terms. Where the pipeline recorded where "
                       "those are stated, the asset's entry links to them.",
    })
    profile = crate.add({
        "@id": PROCESS_RUN_PROFILE,
        "@type": "CreativeWork",
        "name": "Process Run Crate",
        "version": PROCESS_RUN_VERSION,
    })
    # The last time anything in the vault was retrieved or derived; the config's own
    # date for a vault whose records predate those times.
    if dates:
        published = max(dates)
    else:
        from datetime import datetime, timezone
        published = datetime.fromtimestamp(os.path.getmtime(config_file),
                                           timezone.utc).date().isoformat()
    root_entity = {
        "@id": "./",
        "@type": "Dataset",
        "name": "Reference data vault built by refdata-flow",
        "description": "Reference genomes, their annotations, transcriptomes derived "
                       "from them and third-party resources checked against them, held "
                       "in a refgenie vault. Each asset has a provenance.json beside it "
                       "with the full record of where it came from; this crate carries "
                       "the part of that record that has a standard vocabulary.",
        "datePublished": published,
        "license": licence,
        "conformsTo": profile,
        "hasPart": genome_refs,
    }
    if actions:
        root_entity["mentions"] = actions
    descriptor = {
        "@id": "ro-crate-metadata.json",
        "@type": "CreativeWork",
        "conformsTo": ref(RO_CRATE_SPEC),
        "about": ref("./"),
    }
    graph = [descriptor, root_entity] + [crate.entities[key] for key in sorted(crate.entities)]
    # RO-Crate asks for a lone value to be written as the value, not as a list of
    # one. Done here, once, so everything above can build lists without caring.
    for entity in graph:
        for key, value in entity.items():
            if isinstance(value, list) and len(value) == 1:
                entity[key] = value[0]
    return {"@context": CONTEXT, "@graph": graph}


def write_vault_crate(outdir, config_file, logger):
    """Write <outdir>/ro-crate-metadata.json, replacing it in one step."""
    crate = build_crate(outdir, config_file)
    fd, tmp_path = tempfile.mkstemp(dir=outdir, prefix=".ro-crate-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(crate, f, indent=2)
            f.write("\n")
        os.replace(tmp_path, os.path.join(outdir, "ro-crate-metadata.json"))
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    logger.info("Wrote the vault's RO-Crate (%d entities): %s", len(crate["@graph"]),
                os.path.join(outdir, "ro-crate-metadata.json"))
