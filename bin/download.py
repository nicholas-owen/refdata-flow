#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
import argparse
import gzip
import hashlib
import platform
import subprocess
import sys
import json
import os
import re
import shutil
from datetime import datetime, timezone
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
    return release_from_source(match.group(1).strip(), provider, tag)


def release_from_source(source, provider, tag):
    """
    The release an annotation source denotes, as (value, kind), or (None, None).

    `source` is what genomepy records as the 'annotation url': a download URL, or
    for UCSC "UCSC MySQL database: <db>, table: <table>".

    Split out from derive_annotation_release() so the same rules can be applied
    BEFORE a download, to a source genomepy has only been asked about.
    bin/preflight.py does that to decide whether the vault already holds the
    release a provider is currently publishing. One implementation, so the tag the
    preflight looks for is by construction the tag the ingest will have written.
    """
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


# --- what was downloaded: the provenance record -------------------------------
#
# provenance.json used to hold the four request fields, a timestamp and the word
# "genomepy". It could say that a download happened, not what was downloaded: no
# source URL, no checksum, no tool version, nothing about which variant of the
# assembly. Everything below fills that in from three places - the arguments this
# script passed, the README.txt genomepy wrote, and the files themselves.
#
# None of it may fail a download. It runs after gigabytes have been fetched
# correctly, and a record that is missing a field is better than a genome thrown
# away for want of one. Each part is collected separately, and a part that fails
# is named under "provenance_errors" so its absence is not mistaken for nothing
# to report.

# Bumped when the record's shape changes, so a reader can tell an old record from
# an incomplete one. 1 was the request fields only (and was not written down).
PROVENANCE_SCHEMA = 2

GENE_ID_RE = re.compile(rb'gene_id "([^"]+)"')
TRANSCRIPT_ID_RE = re.compile(rb'transcript_id "([^"]+)"')
BIOTYPE_RE = re.compile(rb'gene_(?:bio)?type "([^"]+)"')


def sha256_of_file(path):
    """SHA-256 of a file's bytes as stored, and its size."""
    digest, size = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def describe_fasta(path, sizes_path=None):
    """
    Identify a delivered FASTA two ways, because the two answer different questions.

    file_sha256     the .fa.gz exactly as downloaded and bgzip-compressed here. It
                    changes if the file is recompressed, so it identifies this
                    copy and nothing more.
    content_sha256  the sequence text inside it. Stable under recompression, and
                    the one that can be checked against the vault: refgenie stores
                    the FASTA decompressed, so `sha256sum` of the vault's .fa must
                    give this value.

    Neither is the vault's genome digest, which is computed from sequence names,
    lengths and content and is refgenie's own.
    """
    file_sha256, file_bytes = sha256_of_file(path)
    content, content_bytes = hashlib.sha256(), 0
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            content.update(chunk)
            content_bytes += len(chunk)
    record = {
        "file": os.path.basename(path),
        "file_sha256": file_sha256,
        "file_bytes": file_bytes,
        "content_sha256": content.hexdigest(),
        "content_bytes": content_bytes,
    }
    if sizes_path and os.path.exists(sizes_path):
        with open(sizes_path) as fh:
            lengths = [int(line.split("\t")[1]) for line in fh if line.strip()]
        record["sequences"] = len(lengths)
        record["total_length"] = sum(lengths)
    return record


def describe_gtf(path):
    """
    Identify and summarise a delivered GTF in one pass over its decompressed text.

    The hashes are as for describe_fasta. Read them with one caveat that the FASTA
    does not have: this is the GTF genomepy DELIVERED, not the file the provider
    published. genomepy drops the header lines and rewrites contig names to match
    the genome, and for UCSC there is no upstream file at all - the GTF is
    generated here from database tables. So content_sha256 identifies what is in
    the vault and can be compared between two vaults; it cannot be checked against
    Ensembl's CHECKSUMS or GENCODE's MD5SUMS.

    That is why the summary is recorded beside it. Two annotations that differ
    will have different hashes and nothing more to say; the counts say how they
    differ, which is what someone comparing them needs.
    """
    file_sha256, file_bytes = sha256_of_file(path)
    content, content_bytes, lines = hashlib.sha256(), 0, 0
    features, contigs, biotypes = {}, set(), set()
    genes, transcripts, headers = set(), set(), []

    with gzip.open(path, "rb") as fh:
        for line in fh:
            content.update(line)
            content_bytes += len(line)
            if line.startswith(b"#"):
                if line.startswith(b"#!"):
                    headers.append(line.decode("utf-8", "replace").strip())
                continue
            fields = line.split(b"\t", 8)
            if len(fields) < 9:
                continue
            lines += 1
            contigs.add(fields[0])
            features[fields[2]] = features.get(fields[2], 0) + 1
            attributes = fields[8]
            found = GENE_ID_RE.search(attributes)
            if found:
                genes.add(found.group(1))
            found = TRANSCRIPT_ID_RE.search(attributes)
            if found:
                transcripts.add(found.group(1))
            found = BIOTYPE_RE.search(attributes)
            if found:
                biotypes.add(found.group(1))

    return {
        "file": os.path.basename(path),
        "file_sha256": file_sha256,
        "file_bytes": file_bytes,
        "content_sha256": content.hexdigest(),
        "content_bytes": content_bytes,
        "records": lines,
        "feature_counts": {k.decode(): v for k, v in sorted(features.items())},
        "genes": len(genes),
        "transcripts": len(transcripts),
        "contigs": sorted(c.decode() for c in contigs),
        "biotypes": sorted(b.decode() for b in biotypes),
        # Ensembl and GENCODE state provider and release in '#!' lines. genomepy
        # 0.16.4 strips them, so this is normally empty; kept in case that changes.
        "header_lines": headers,
    }


def parse_genomepy_readme(path):
    """
    genomepy's README.txt as a dict, plus its free text.

    genomepy already records where everything came from - the source URLs, the
    assembly accession, the taxon, the masking, and what it changed - but as text
    nobody can query. 'key: value' lines become fields with the spaces in their
    keys turned to underscores; anything else (for instance "No contigs were
    removed.") is kept under 'notes'.
    """
    fields, notes = {}, []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            key, separator, value = line.partition(": ")
            if separator and re.fullmatch(r"[A-Za-z][A-Za-z _]*", key):
                fields[key.strip().lower().replace(" ", "_")] = value.strip()
            else:
                notes.append(line)
    if notes:
        fields["notes"] = notes
    return fields


def tool_versions():
    """
    Versions of what shaped the delivered files.

    bgzip and samtools are asked directly because they are whichever ones are on
    PATH in this task, not a property of this script. The UCSC converters that
    build an annotation from UCSC's tables have no version flag to ask.
    """
    versions = {
        "python": platform.python_version(),
        "genomepy": getattr(genomepy, "__version__", "unknown"),
    }
    for tool in ("bgzip", "samtools"):
        try:
            out = subprocess.run([tool, "--version"], capture_output=True, text=True,
                                 timeout=30).stdout.splitlines()
            if out:
                versions[tool] = out[0].strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return versions


def describe_download(assembly_dir, assembly, annotation, tag):
    """
    Everything the record says about the files, as (fields, errors).

    Each part is attempted on its own; see the note at the top of this section.
    """
    fields, errors = {}, []

    def attempt(name, build):
        try:
            build()
        except Exception as e:
            errors.append(f"{name}: {e}")
            print(f"{tag}       WARNING: could not record {name} in provenance.json: {e}",
                  file=sys.stderr, flush=True)

    def readme():
        parsed = parse_genomepy_readme(os.path.join(assembly_dir, "README.txt"))
        fields["genomepy_readme"] = parsed
        # Hoisted, because these are what a reader looks for first.
        for key in ("assembly_accession", "tax_id", "genome_url", "annotation_url"):
            if parsed.get(key):
                fields[key] = parsed[key]

    def genome():
        fields["genome"] = describe_fasta(
            os.path.join(assembly_dir, f"{assembly}.fa.gz"),
            os.path.join(assembly_dir, f"{assembly}.fa.sizes"))

    def gtf():
        fields["annotation"] = describe_gtf(
            os.path.join(assembly_dir, f"{assembly}.annotation.gtf.gz"))

    attempt("genomepy's README.txt", readme)
    attempt("the genome's checksums", genome)
    if annotation:
        attempt("the annotation's checksums and summary", gtf)
    attempt("tool versions", lambda: fields.__setitem__("tool_versions", tool_versions()))
    return fields, errors


def main():
    parser = argparse.ArgumentParser(description="Wrapper for genomepy for Nextflow pipeline.")
    parser.add_argument("--species", required=True, help="Species name (e.g., Homo_sapiens)")
    parser.add_argument("--assembly", required=True, help="Assembly name (e.g., GRCh38)")
    parser.add_argument("--provider", required=True, help="Provider (e.g., ensembl, ncbi, ucsc)")
    parser.add_argument("--annotation", action="store_true", help="Download annotation")
    parser.add_argument("--release", default=None,
                        help="Provider release to pin (Ensembl only)")
    parser.add_argument("--pipeline-version", default="",
                        help="refdata-flow version, recorded in provenance.json")
    parser.add_argument("--url", default="",
                        help="Address of the genome FASTA (provider 'url' only)")
    parser.add_argument("--annotation-url", default="",
                        help="Address of the annotation (provider 'url' only)")

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

    # The arguments that decide WHICH variant of the assembly is delivered, and so
    # what the vault's digest will be. They are genomepy's defaults, passed
    # explicitly: the record below is then a copy of what was asked for, not a
    # belief about what a default was in the genomepy version that happened to run.
    #
    #   mask          soft-masked sequence where the provider offers it
    #   keep_alt      False: alternative haplotypes and patches are filtered out
    #   invert_match  no regex is given, so nothing else is filtered
    install_arguments = {
        "mask": "soft",
        "keep_alt": False,
        "invert_match": False,
        "bgzip": True,
    }
    kwargs = {
        "name": args.assembly,
        "provider": args.provider,
        "genomes_dir": local_dir,
        "annotation": args.annotation,
        **install_arguments,
    }

    # A genome requested by address. genomepy's 'url' provider takes the address
    # where every other provider takes a name, and the name the genome is to be
    # known by as `localname`; the directory and files it writes are then called
    # after that name, exactly as for a named assembly, so nothing downstream
    # needs to know the difference. With no provider to ask, an annotation must be
    # given by address too: left to itself genomepy goes looking in the genome's
    # directory for a file with a similar name, which is a guess.
    is_url = args.provider.strip().lower() == "url"
    if is_url:
        if not args.url:
            print(f"{tag} Refusing to run: provider 'url' needs --url.",
                  file=sys.stderr, flush=True)
            sys.exit(1)
        if args.annotation and not args.annotation_url:
            print(f"{tag} Refusing to run: annotation was requested for a genome given "
                  "by address, and no --annotation-url says where it is.",
                  file=sys.stderr, flush=True)
            sys.exit(1)
        kwargs["name"] = args.url
        kwargs["localname"] = args.assembly
        if args.annotation:
            kwargs["to_annotation"] = args.annotation_url
        print(f"{tag}       Taking the genome from {args.url}", flush=True)
    elif args.url or args.annotation_url:
        print(f"{tag} Refusing to run: --url and --annotation-url are only for provider "
              f"'url', not '{args.provider}'. The address would be ignored and the "
              "provider's own file downloaded in its place.", file=sys.stderr, flush=True)
        sys.exit(1)

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
            "schema": PROVENANCE_SCHEMA,
            # UTC, and said to be. This was naive local time, which on a vault
            # shared between machines could not be placed. For UCSC it is also
            # the nearest thing to a version the sequence has.
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "pipeline_version": args.pipeline_version,
            "species": args.species,
            "assembly": args.assembly,
            "provider": args.provider,
            "annotation_requested": args.annotation,
            # The release the request pinned, empty when it took the current one.
            # With assembly and provider this identifies which request a download
            # directory answers; bin/update_refgenie.py checks all three before
            # ingesting a directory it was handed by position.
            "release_requested": release,
            "tool": "genomepy",
            "install_arguments": install_arguments,
        }
        # For a genome requested by address, the address IS the request. Recorded
        # as asked for, separately from the 'genome_url' genomepy reports, so the
        # preflight can tell when a name in the vault came from somewhere else.
        if is_url:
            provenance["url_requested"] = args.url
            provenance["annotation_url_requested"] = args.annotation_url if args.annotation else ""

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

        # What the files are: checksums, a summary of the annotation, genomepy's
        # own record, tool versions. A full read of both files, so for a large
        # genome this is a minute or two; it cannot fail the download.
        print(f"{tag} [3/3] Checksumming and writing provenance.json", flush=True)
        described, errors = describe_download(args.assembly, args.assembly,
                                              args.annotation, tag)
        provenance.update(described)
        if errors:
            provenance["provenance_errors"] = errors

        prov_path = os.path.join(args.assembly, "provenance.json")
        with open(prov_path, "w") as f:
            json.dump(provenance, f, indent=4)

        print(f"{tag} Done.", flush=True)

    except Exception as e:
        print(f"{tag} Error downloading genome: {e}", file=sys.stderr, flush=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
