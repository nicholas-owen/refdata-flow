#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import sys
import json
import os
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

def main():
    parser = argparse.ArgumentParser(description="Wrapper for genomepy for Nextflow pipeline.")
    parser.add_argument("--species", required=True, help="Species name (e.g., Homo_sapiens)")
    parser.add_argument("--assembly", required=True, help="Assembly name (e.g., GRCh38)")
    parser.add_argument("--provider", required=True, help="Provider (e.g., ensembl, ncbi, ucsc)")
    parser.add_argument("--annotation", action="store_true", help="Download annotation")
    
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
