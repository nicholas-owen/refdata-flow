#!/usr/bin/env bash
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Derive the small files that pipelines ask for alongside a GTF, from one
# annotation in the vault (the DERIVE_ANNOTATION process).
#
#   derive_annotation.sh GENOME ANNOTATION_ASSET RELEASE GTF OUTDIR [PIPELINE_VERSION]
#
# Four assets come out, each with a provenance record beside it:
#
#   tx2gene       tx2gene.tsv            transcript_id, gene_id, gene_name; no header
#   gtf_indexed   sorted.gtf.gz + .tbi   the same GTF, position-sorted, bgzipped, tabix-indexed
#   bed12         transcripts.bed        one BED12 line per transcript
#   refflat       refFlat.txt            the same gene models as Picard reads them
#
# Everything is taken from the GTF alone; the genome is not read. So the files say
# exactly what the annotation says, for every provider, and a transcript here is a
# transcript in the vault's transcriptome under the same name: both use the GTF's
# transcript_id as written.
#
# DECISIONS WORTH KNOWING
#
#   * tx2gene's third column falls back to the gene_id where a gene has no name.
#     In Ensembl human 116 that is 35,483 of 78,941 genes, and an empty field
#     breaks more readers than a repeated one.
#
#   * refFlat's first column is the gene_id, not the gene name. Picard groups
#     transcripts into genes by that column, and gene names are not unique: the
#     same release has 758 genes called Y_RNA, across the chromosomes, and 484
#     names in all that more than one gene carries. A gene whose transcripts sit
#     on different chromosomes or strands is one Picard cannot build.
#     tx2gene.tsv maps the id to the name for anything that wants one.
#
#   * BED12 comes from UCSC's gtfToGenePred and genePredToBed, not from
#     `gffread --bed`. Checked on Ensembl human 116: gffread marks the whole of a
#     non-coding transcript as coding sequence (thickStart..thickEnd), where the
#     UCSC convention - and what RSeQC reads - is an empty range for a transcript
#     with no CDS.
#
#   * The sort is stable and by sequence name then start only, so lines that start
#     together keep the order the provider wrote them in: a gene stays ahead of
#     its transcripts, and a transcript ahead of its exons.
#
# Needs gtfToGenePred, genePredToBed, bgzip and tabix on PATH. The conda profile's
# default environment has all four, because genomepy needs them too.
set -euo pipefail
export LC_ALL=C

if [ "$#" -lt 5 ]; then
    echo "usage: derive_annotation.sh GENOME ANNOTATION_ASSET RELEASE GTF OUTDIR [PIPELINE_VERSION]" >&2
    exit 2
fi
genome=$1
annotation_asset=$2
release=$3
gtf=$4
out=$5
pipeline_version=${6:-}
tab=$(printf '\t')

for tool in gtfToGenePred genePredToBed bgzip tabix; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[$genome] ERROR: $tool is not on PATH. The conda profile supplies it; outside conda, install it with: conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed htslib" >&2
        exit 3
    fi
done

# The version of a conda package, read from the environment the tool was found
# in. The UCSC tools have no --version to ask.
conda_version() {  # package, a binary it provides
    local bin meta found
    bin=$(command -v "$2")
    meta="$(dirname "$bin")/../conda-meta"
    found=$(ls "$meta"/"$1"-[0-9]*.json 2>/dev/null | head -n 1 || true)
    if [ -z "$found" ]; then
        echo unknown
    else
        basename "$found" .json | sed -E "s/^$1-//; s/-[^-]+\$//"
    fi
}

mkdir -p "$out"
gzip -dcf "$gtf" > annotation.gtf

# ---- tx2gene -----------------------------------------------------------------
# One line per transcript_id, in the order the GTF first names it. Read from every
# line that carries one, not only `transcript` lines: a GTF built from UCSC's
# tables may have none.
awk -F'\t' '
    function attribute(text, key,    value) {
        if (!match(text, "(^|;)[ ]*" key " \"[^\"]+\"")) return ""
        value = substr(text, RSTART, RLENGTH)
        sub(/^[^"]*"/, "", value)
        sub(/"$/, "", value)
        return value
    }
    /^#/ || NF < 9 { next }
    {
        transcript = attribute($9, "transcript_id")
        if (transcript == "" || (transcript in seen)) next
        seen[transcript] = 1
        gene = attribute($9, "gene_id")
        name = attribute($9, "gene_name")
        print transcript "\t" gene "\t" (name == "" ? gene : name)
    }' annotation.gtf > "$out/tx2gene.tsv"

transcripts=$(wc -l < "$out/tx2gene.tsv" | tr -d ' ')
genes=$(cut -f 2 "$out/tx2gene.tsv" | sort -u | wc -l | tr -d ' ')
if [ "$transcripts" -eq 0 ]; then
    echo "[$genome] ERROR: no line of $annotation_asset:$release carries a transcript_id, so there is nothing to derive from it." >&2
    exit 1
fi
if [ "$(awk -F'\t' '$2 == ""' "$out/tx2gene.tsv" | wc -l | tr -d ' ')" -ne 0 ]; then
    echo "[$genome] ERROR: $annotation_asset:$release has transcripts with no gene_id." >&2
    exit 1
fi

# ---- sorted, indexed GTF -------------------------------------------------------
{ grep '^#' annotation.gtf || true; } > sorted.gtf
{ grep -v -e '^#' -e '^$' annotation.gtf || true; } \
    | sort -s -t "$tab" -k1,1 -k4,4n -S 1G -T . >> sorted.gtf
bgzip -c sorted.gtf > "$out/sorted.gtf.gz"
tabix -p gff "$out/sorted.gtf.gz"

records=$(grep -c -v '^#' sorted.gtf || true)
sequences=$(grep -v '^#' sorted.gtf | cut -f 1 | uniq | sort -u | wc -l | tr -d ' ')
indexed=$(tabix -l "$out/sorted.gtf.gz" | wc -l | tr -d ' ')
if [ "$indexed" -ne "$sequences" ]; then
    echo "[$genome] ERROR: the tabix index lists $indexed sequence(s) but the GTF has records on $sequences." >&2
    exit 1
fi

# ---- gene models: BED12 and refFlat ---------------------------------------------
# -genePredExt adds the gene_id as column 12. -ignoreGroupsWithoutExons: a
# transcript with no exon line has no model to write, and without the flag one
# such transcript stops the whole conversion.
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons annotation.gtf annotation.genePred
genePredToBed annotation.genePred unsorted.bed
sort -t "$tab" -k1,1 -k2,2n -k3,3n -k4,4 unsorted.bed > "$out/transcripts.bed"
awk -F'\t' 'BEGIN { OFS = "\t" } { print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' annotation.genePred \
    | sort -t "$tab" -k3,3 -k5,5n -k6,6n -k2,2 > "$out/refFlat.txt"

models=$(wc -l < "$out/transcripts.bed" | tr -d ' ')
if [ "$models" -eq 0 ]; then
    echo "[$genome] ERROR: gtfToGenePred wrote no gene models from $annotation_asset:$release." >&2
    exit 1
fi
if [ "$models" -gt "$transcripts" ]; then
    echo "[$genome] ERROR: $models gene models from $transcripts transcripts in $annotation_asset:$release." >&2
    exit 1
fi
if [ "$models" -ne "$transcripts" ]; then
    echo "[$genome] WARNING: $((transcripts - models)) of $transcripts transcripts in $annotation_asset:$release have no exon line and so no gene model. They are in tx2gene and absent from bed12 and refflat." >&2
fi

# ---- the records ---------------------------------------------------------------
ucsc_version=$(conda_version ucsc-gtftogenepred gtfToGenePred)
htslib_version=$(bgzip --version 2>/dev/null | head -n 1 | awk '{ print $NF }' || true)
now=$(date -u +%Y-%m-%dT%H:%M:%S+00:00)

record() {  # asset, tool, tool version, command, further JSON fields
    printf '{"asset": "%s", "genome": "%s", "annotation": "%s:%s", "tool": "%s", "tool_version": "%s", "command": "%s", "derived_utc": "%s", %s}\n' \
        "$1" "$genome" "$annotation_asset" "$release" "$2" "$3" "$4" "$now" "$5" > "$out/$1.provenance.json"
}
record tx2gene refdata-flow "$pipeline_version" \
    "bin/derive_annotation.sh: transcript_id, gene_id and gene_name of each transcript, in GTF order" \
    "\"columns\": [\"transcript_id\", \"gene_id\", \"gene_name\"], \"header\": false, \"gene_name_fallback\": \"gene_id\", \"transcripts\": $transcripts, \"genes\": $genes"
record gtf_indexed htslib "${htslib_version:-unknown}" \
    "sort -s -k1,1 -k4,4n ANNOTATION.gtf | bgzip; tabix -p gff" \
    "\"records\": $records, \"sequences\": $sequences"
record bed12 genePredToBed "$ucsc_version" \
    "gtfToGenePred -genePredExt -ignoreGroupsWithoutExons ANNOTATION.gtf A.genePred; genePredToBed A.genePred; sort -k1,1 -k2,2n" \
    "\"transcripts\": $models, \"transcripts_in_annotation\": $transcripts"
record refflat gtfToGenePred "$ucsc_version" \
    "gtfToGenePred -genePredExt -ignoreGroupsWithoutExons ANNOTATION.gtf A.genePred; columns 12,1-10; sort -k3,3 -k5,5n" \
    "\"gene_column\": \"gene_id\", \"transcripts\": $models, \"transcripts_in_annotation\": $transcripts"

echo "[$genome] $annotation_asset:$release -> $transcripts transcripts of $genes genes; $models gene models; $records GTF records on $sequences sequences"
