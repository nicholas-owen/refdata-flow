#!/usr/bin/env bash
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Write a genome in UCSC's .2bit format, from the FASTA in the vault (the
# DERIVE_TWOBIT process).
#
#   derive_twobit.sh GENOME FASTA FAI OUT_2BIT OUT_PROVENANCE
#
# .2bit is the genome as deepTools, BLAT, liftOver's chain tools and rtracklayer
# read it: two bits a base, random access by sequence name. Soft-masking is kept:
# lower-case sequence in the FASTA is recorded as masked in the .2bit.
#
# The file is not trusted because faToTwoBit exited 0. The names and lengths in it
# are read back with twoBitInfo and must be the genome's own, as its .fai has them.
#
# Needs faToTwoBit and twoBitInfo on PATH; the conda profile gives this process an
# environment with both.
set -euo pipefail
export LC_ALL=C

if [ "$#" -ne 5 ]; then
    echo "usage: derive_twobit.sh GENOME FASTA FAI OUT_2BIT OUT_PROVENANCE" >&2
    exit 2
fi
genome=$1
fasta=$2
fai=$3
out=$4
provenance=$5

for tool in faToTwoBit twoBitInfo; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[$genome] ERROR: $tool is not on PATH. The conda profile supplies it; outside conda, install it with: conda install -c bioconda ucsc-fatotwobit ucsc-twobitinfo" >&2
        exit 3
    fi
done

sequences=$(wc -l < "$fai" | tr -d ' ')
bases=$(awk -F'\t' '{ total += $2 } END { printf "%.0f", total }' "$fai")

# The format's index holds 32-bit offsets; a genome past four gigabases needs the
# 64-bit variant, which older readers cannot open, so it is used only then.
long=""
if [ "$bases" -gt 4000000000 ]; then
    long="-long"
fi
# shellcheck disable=SC2086  # $long is one flag or nothing
faToTwoBit $long "$fasta" "$out"

twoBitInfo "$out" stdout | sort > twobit.sizes
cut -f 1,2 "$fai" | sort > genome.sizes
if ! cmp -s twobit.sizes genome.sizes; then
    echo "[$genome] ERROR: the .2bit does not hold the genome's sequences. Names or lengths that differ (genome '<', .2bit '>'):" >&2
    diff genome.sizes twobit.sizes | grep '^[<>]' | head -n 10 >&2 || true
    exit 1
fi

bin=$(command -v faToTwoBit)
found=$(ls "$(dirname "$bin")/../conda-meta"/ucsc-fatotwobit-[0-9]*.json 2>/dev/null | head -n 1 || true)
if [ -z "$found" ]; then
    version=unknown
else
    version=$(basename "$found" .json | sed -E 's/^ucsc-fatotwobit-//; s/-[^-]+$//')
fi

printf '{"asset": "twobit", "genome": "%s", "tool": "faToTwoBit", "tool_version": "%s", "command": "faToTwoBit %sGENOME.fa GENOME.2bit; twoBitInfo GENOME.2bit, compared with GENOME.fa.fai", "derived_utc": "%s", "soft_masking": "kept", "sequences": %s, "bases": %s}\n' \
    "$genome" "$version" "${long:+$long }" "$(date -u +%Y-%m-%dT%H:%M:%S+00:00)" "$sequences" "$bases" > "$provenance"

echo "[$genome] 2bit: $sequences sequences, $bases bases, names and lengths match the genome"
