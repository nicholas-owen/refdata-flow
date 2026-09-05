#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import sys
import json
import os
import re
import shutil
from datetime import datetime
import genomepy

# UCSC command-line tools that genomepy shells out to when converting gene
# annotations between GTF / GFF3 / BED / genePred. They ship with the Bioconda
# genomepy package but NOT with the pip one, so a pip-only environment can
# download sequence perfectly well and then fail on every annotated request.
#
# GATE_TOOL is the single tool genomepy's own check_ucsc_tools() tests
# (genomepy/utils.py). We fail on exactly that one so this preflight can never
# reject a run genomepy would have accepted; the others are reported as a warning
# because a partial install is a latent failure worth knowing about.
UCSC_ANNOTATION_TOOLS = (
    "genePredToGtf",
    "gtfToGenePred",
    "genePredToBed",
    "bedToGenePred",
    "gff3ToGenePred",
)
GATE_TOOL = "genePredToGtf"

# Providers whose downloads can be pinned to a release. Kept in step with
# PINNABLE_PROVIDERS in bin/resolve.py, which is where the request is validated;
# this copy exists so download.py cannot be made to ignore a release silently when
# invoked directly. Verified against genomepy 0.16.4: only providers/ensembl.py
# reads kwargs["version"].
PINNABLE_PROVIDERS = {"ensembl"}

# Exit code reserved for "a required external tool is missing". Distinct from 1
# (the download itself failed) because a missing binary is not worth retrying -
# useful once main.nf grows an errorStrategy.
EXIT_MISSING_TOOL = 3


def check_annotation_tools(tag):
    """
    Verify the UCSC tools are present before anything is downloaded.

    genomepy performs this check itself, but only once it reaches the annotation
    step - by which point the genome FASTA has already been fetched and filtered.
    For a human assembly that is several GB and many minutes of work thrown away,
    so the same check is made here, up front, and the message names the fix.
    """
    missing = [t for t in UCSC_ANNOTATION_TOOLS if shutil.which(t) is None]
    if not missing:
        return

    if GATE_TOOL in missing:
        print(
            f"{tag} ERROR: annotation was requested but the UCSC tools genomepy needs "
            f"are not on PATH (missing: {', '.join(missing)}).\n"
            f"{tag}        These ship with the Bioconda genomepy package but not the pip one.\n"
            f"{tag}        Install them with:\n"
            f"{tag}          conda install -c bioconda "
            + " ".join("ucsc-" + t.lower() for t in UCSC_ANNOTATION_TOOLS) + "\n"
            f"{tag}        or see https://github.com/vanheeringen-lab/genomepy#pip for direct\n"
            f"{tag}        download links. Re-run without --annotation to fetch sequence only.",
            file=sys.stderr, flush=True,
        )
        sys.exit(EXIT_MISSING_TOOL)

    print(
        f"{tag} WARNING: some UCSC annotation tools are missing ({', '.join(missing)}). "
        f"'{GATE_TOOL}' is present so genomepy will proceed, but conversions needing the "
        "others will fail.",
        file=sys.stderr, flush=True,
    )


# --- annotation release identification (roadmap 1b, Half B) ------------------
#
# The genome digest covers the FASTA and nothing else, so two vaults with the
# same digest can hold annotations from different releases and `refgenie seek`
# returns a confident path either way. The fix is to tag the annotation asset,
# which first requires knowing what release was actually downloaded.
#
# genomepy renames its downloads to '<localname>.annotation.gtf.gz', so the
# release is NOT recoverable from the filename. It IS recorded in README.txt,
# which genomepy writes with the URL it actually used:
#
#   Ensembl  annotation url: http://ftp.ensembl.org/pub/release-116/gtf/...
#   UCSC     annotation url: UCSC MySQL database: sacCer3, table: ncbiRefSeq
#
# Each provider therefore needs its own rule, because each publishes a different
# kind of identifier. We record the most specific *stable* one available rather
# than forcing a single vocabulary: an Ensembl release number says more than the
# date that release happened to be published.
README_ANNOTATION_RE = re.compile(r"^annotation url:\s*(.+)$", re.MULTILINE)
ENSEMBL_RELEASE_RE = re.compile(r"/release-(\d+)/")
GENCODE_RELEASE_RE = re.compile(r"/release_(M?\d+)")
UCSC_MYSQL_RE = re.compile(r"UCSC MySQL database:\s*(\w+),\s*table:\s*(\w+)")

# refgenie parses 'genome/asset:tag', so a tag may not contain a colon. Dates are
# therefore YYYYMMDD rather than ISO timestamps. Same-day collisions are not a
# concern: a tag is scoped to one (genome, asset), so a collision would need the
# same provider to republish the same annotation for the same assembly twice
# within a day.
TAG_SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]")


def _ucsc_table_updated(database, table, tag):
    """
    When did UCSC last change the annotation table genomepy actually read?

    UCSC has no release number - it patches assemblies in place - so the honest
    identifier is the modification time of the source itself. Note this must come
    from MySQL and not from the published file dump under goldenPath/<db>/database/:
    the two disagree (for sacCer3.sgdGene, 2011-08-29 in MySQL versus 2011-10-06
    for the dump), and MySQL is what genomepy reads.

    Returns YYYYMMDD, or None if the server is unreachable or reports no time.
    """
    try:
        import mysql.connector

        conn = mysql.connector.connect(
            host="genome-mysql.soe.ucsc.edu", user="genome",
            database=database, connection_timeout=30,
        )
        try:
            cur = conn.cursor()
            cur.execute(
                "SELECT UPDATE_TIME, CREATE_TIME FROM information_schema.tables "
                "WHERE TABLE_SCHEMA=%s AND TABLE_NAME=%s",
                (database, table),
            )
            row = cur.fetchone()
            cur.close()
        finally:
            conn.close()
    except Exception as e:
        print(f"{tag}       Could not reach UCSC MySQL for a table date: {e}",
              file=sys.stderr, flush=True)
        return None

    if not row:
        return None
    # UPDATE_TIME is storage-engine dependent: MyISAM persists it, InnoDB can
    # report NULL or reset it on restart. CREATE_TIME is the documented fallback,
    # and for UCSC's static tables the two are usually identical anyway.
    stamp = row[0] or row[1]
    return stamp.strftime("%Y%m%d") if stamp else None


def _http_last_modified(url, tag):
    """Upstream's Last-Modified for `url`, as YYYYMMDD, or None."""
    try:
        import requests
        from email.utils import parsedate_to_datetime

        resp = requests.head(url, timeout=30, allow_redirects=True)
        resp.raise_for_status()
        header = resp.headers.get("Last-Modified")
        return parsedate_to_datetime(header).strftime("%Y%m%d") if header else None
    except Exception as e:
        print(f"{tag}       Could not read Last-Modified from {url}: {e}",
              file=sys.stderr, flush=True)
        return None


def derive_annotation_release(assembly_dir, provider, tag):
    """
    Identify the annotation release, as (value, kind), or (None, None).

    Deliberately never raises. This runs after a successful download, and losing
    a provenance label is not a reason to discard gigabytes of correctly fetched
    sequence. An unidentified release means the asset is tagged 'default', which
    is exactly the behaviour before this existed.
    """
    readme = os.path.join(assembly_dir, "README.txt")
    try:
        with open(readme) as fh:
            text = fh.read()
    except OSError as e:
        print(f"{tag}       No README.txt to read a release from: {e}",
              file=sys.stderr, flush=True)
        return None, None

    match = README_ANNOTATION_RE.search(text)
    if not match:
        return None, None
    source = match.group(1).strip()

    name = str(provider).strip().lower()

    if name == "ensembl":
        found = ENSEMBL_RELEASE_RE.search(source)
        return (found.group(1), "ensembl_release") if found else (None, None)

    if name == "gencode":
        # Not pinnable through genomepy, but the release is still in the path and
        # is a better identifier than the date it was published.
        found = GENCODE_RELEASE_RE.search(source)
        return (found.group(1), "gencode_release") if found else (None, None)

    if name == "ucsc":
        found = UCSC_MYSQL_RE.search(source)
        if not found:
            return None, None
        stamp = _ucsc_table_updated(found.group(1), found.group(2), tag)
        return (stamp, f"ucsc_table_update:{found.group(2)}") if stamp else (None, None)

    if source.startswith("http"):
        # NCBI and anything else file-based. The accession version pins the
        # assembly, not the annotation, so the file's own date is what we have.
        stamp = _http_last_modified(source, tag)
        return (stamp, "http_last_modified") if stamp else (None, None)

    return None, None


def safe_tag(value):
    """Reduce a release identifier to something refgenie will accept as a tag."""
    cleaned = TAG_SAFE_RE.sub("-", str(value)).strip("-")
    return cleaned or None


def main():
    parser = argparse.ArgumentParser(description="Wrapper for genomepy for Nextflow pipeline.")
    parser.add_argument("--species", required=True, help="Species name (e.g., Homo_sapiens)")
    parser.add_argument("--assembly", required=True, help="Assembly name (e.g., GRCh38)")
    parser.add_argument("--provider", required=True, help="Provider (e.g., ensembl, ncbi, ucsc)")
    parser.add_argument("--annotation", action="store_true", help="Download annotation")
    parser.add_argument("--release", default=None,
                        help="Provider release to pin (Ensembl only)")

    args = parser.parse_args()

    # We download into the current working directory.
    # genomepy will automatically create a directory named after the assembly (e.g. GRCh38).
    local_dir = "."

    # Progress messages are printed to stdout with flush=True so they stream live
    # to the Nextflow console (the process sets `debug true`). Each line is tagged
    # with the assembly so parallel FETCH_GENOME tasks remain attributable when
    # their output interleaves. genomepy's own byte-level progress bar goes to
    # stderr and is captured in the task's .command.log (tail it for detail).
    tag = f"[{args.assembly}]"
    artefacts = "FASTA + annotation (GTF)" if args.annotation else "FASTA"

    # Fail before a single byte is downloaded if the annotation toolchain is
    # incomplete. genomepy only checks once it reaches the annotation step, which
    # for a large assembly is minutes of downloading and filtering already spent.
    if args.annotation:
        check_annotation_tools(tag)

    print(f"{tag} [1/3] Starting: {args.species} {args.assembly} "
          f"from {args.provider} - fetching {artefacts}", flush=True)

    kwargs = {
        "name": args.assembly,
        "provider": args.provider,
        "genomes_dir": local_dir,
        "annotation": args.annotation,
        "bgzip": True
    }

    # Only Ensembl reads `version` out of install_genome's **kwargs; the other
    # providers accept it and discard it without complaint. bin/resolve.py already
    # refuses that combination, so reaching here with one would be a bug - guard at
    # the call site anyway, because the failure it prevents is silent.
    release = (args.release or "").strip()
    if release:
        if args.provider.strip().lower() in PINNABLE_PROVIDERS:
            kwargs["version"] = release
            print(f"{tag}       Pinning to {args.provider} release {release}.", flush=True)
        else:
            print(f"{tag} Refusing to run: release '{release}' was passed for provider "
                  f"'{args.provider}', which cannot be pinned. genomepy would ignore it "
                  "and download whatever is current, recording a release that was never "
                  "requested.", file=sys.stderr, flush=True)
            sys.exit(1)

    try:
        # Install the genome using genomepy's Python API. This is the slow step:
        # it downloads the sequence (and annotation, if requested) and bgzip-
        # compresses it. Detailed progress is in this task's .command.log.
        print(f"{tag} [2/3] Downloading and bgzip-compressing via genomepy "
              f"(slow step - see .command.log for the detailed progress bar)...", flush=True)
        genomepy.install_genome(**kwargs)
        print(f"{tag}       Download + compression complete.", flush=True)

        # Provenance tracking: capture metadata of the download session
        provenance = {
            "timestamp": datetime.now().isoformat(),
            "species": args.species,
            "assembly": args.assembly,
            "provider": args.provider,
            "annotation_requested": args.annotation,
            "tool": "genomepy"
        }

        # Identify which annotation release this actually is, so the ingest stage
        # can tag the asset with it. Recorded here rather than at ingest because
        # this is where README.txt and the network are both to hand, and because
        # provenance.json is already preserved by update_refgenie.py's
        # preserve_provenance() before any cleanup.
        if args.annotation:
            release, kind = derive_annotation_release(args.assembly, args.provider, tag)
            release = safe_tag(release) if release else None
            if release:
                provenance["annotation_release"] = release
                provenance["annotation_release_kind"] = kind
                print(f"{tag}       Annotation release identified: {release} ({kind})",
                      flush=True)
            else:
                # Not fatal. The asset falls back to the 'default' tag, which is
                # the behaviour that existed before release tagging.
                print(f"{tag}       WARNING: could not identify the annotation release "
                      f"for provider '{args.provider}'; the asset will be tagged "
                      "'default' and two releases could later collide under it.",
                      file=sys.stderr, flush=True)

        # Write the provenance.json into the created assembly directory
        print(f"{tag} [3/3] Writing provenance.json", flush=True)
        prov_path = os.path.join(args.assembly, "provenance.json")
        with open(prov_path, "w") as f:
            json.dump(provenance, f, indent=4)

        print(f"{tag} Done.", flush=True)

    except Exception as e:
        print(f"{tag} Error downloading genome: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
