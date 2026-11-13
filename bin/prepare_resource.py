#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Decide whether a fetched resource belongs to the genome it was requested for,
and if it does, make it ready to register.

THE GATE

A resource is a list of positions on named sequences. Whether those names and
positions mean anything depends entirely on the genome they are read against,
and nothing in a tool that takes the two together will say when they do not
match. A known-sites VCF that calls its contigs chr1 beside a genome that calls
them 1 is not an error to GATK: it finds no known site anywhere, and recalibrates
base qualities as though none existed. The run completes.

So nothing reaches the vault until it has been compared with the genome's own
chrom.sizes - the names and lengths of the sequences actually held. Three checks,
each of which refuses the resource outright when it fails:

  names        The share of the resource's records that sit on a sequence the
               genome has must reach the resource's `min_matched`. A
               differently-named file scores zero. A file with a few records on
               unplaced scaffolds the genome leaves out scores just under one,
               and passes, with those sequences and their counts recorded.

  coordinates  No record may end beyond the end of its sequence. This is what
               catches the right names on the wrong build: GRCh37 and GRCh38
               both have a chromosome 1, of different lengths.

  lengths      Where the file declares sequence lengths itself (a VCF's ##contig
               lines), they must equal the genome's.

None of these proves a match. Two builds can agree on every name and length that
a small file happens to touch. They are what can be checked, they are cheap, and
each one that fails is a certainty.

AFTERWARDS

A resource whose registry entry asks for a tabix index gets one, built here by
pysam's bundled htslib, after re-compressing to BGZF if it was published as plain
gzip. One that asks for none is stored as fetched.

Writes, to --outdir: the file, as <asset>.data.gz; its .tbi if indexed;
provenance.json, which is the fetch record plus what the gate found; and TAG,
the version to register it under.
"""
import argparse
import gzip
import json
import os
import shutil
import sys

import resource_registry


def read_chrom_sizes(path):
    """{sequence name: length} from a chrom.sizes file."""
    sizes = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0]:
                sizes[fields[0]] = int(fields[1])
    return sizes


def inspect(path, kind):
    """
    One pass over a gzipped VCF or BED. Returns what the gate needs:

        records           number of data lines
        per_contig        {name: [records, furthest end position]}
        declared_lengths  {name: length} from ##contig lines (VCF only)
        header            the VCF's ## lines other than INFO/FILTER/FORMAT/contig,
                          which is where a file says what it is

    A BED interval's end is its third column. A VCF record is taken to end at its
    position: enough to tell a wrong build, without parsing REF or END.
    """
    records, per_contig, declared, header = 0, {}, {}, []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                if kind == "vcf" and line.startswith("##"):
                    if line.startswith("##contig=<"):
                        body = dict(
                            part.split("=", 1) for part in
                            line.strip()[len("##contig=<"):-1].split(",") if "=" in part
                        )
                        if body.get("ID") and str(body.get("length", "")).isdigit():
                            declared[body["ID"]] = int(body["length"])
                    elif not line.startswith(("##INFO", "##FILTER", "##FORMAT", "##ALT")):
                        header.append(line.strip())
                continue
            if kind == "bed" and line.startswith(("track", "browser")):
                continue
            fields = line.split("\t", 3)
            try:
                end = int(fields[2]) if kind == "bed" else int(fields[1])
            except (IndexError, ValueError):
                continue
            records += 1
            entry = per_contig.setdefault(fields[0], [0, 0])
            entry[0] += 1
            if end > entry[1]:
                entry[1] = end
    return {"records": records, "per_contig": per_contig,
            "declared_lengths": declared, "header": header}


def gate(found, sizes, min_matched):
    """
    Compare what a resource contains with the genome's sequences.

    Returns (reasons it cannot be registered, summary). An empty list means it
    passed; the summary is recorded either way.
    """
    reasons = []
    per_contig = found["per_contig"]
    total = found["records"]
    matched = sum(count for name, (count, _) in per_contig.items() if name in sizes)
    absent = {name: count for name, (count, _) in per_contig.items() if name not in sizes}
    share = (matched / total) if total else 0.0

    summary = {
        "records": total,
        "contigs_in_resource": len(per_contig),
        "records_on_genome_contigs": matched,
        "share_matched": round(share, 6),
        "min_matched": min_matched,
        "contigs_absent_from_genome": dict(sorted(absent.items())),
    }

    if total == 0:
        reasons.append("it contains no records")
        return reasons, summary

    if share < min_matched:
        ours = ", ".join(list(sizes)[:4])
        theirs = ", ".join(sorted(per_contig, key=lambda n: -per_contig[n][0])[:4])
        reasons.append(
            f"only {share:.1%} of its {total:,} records are on sequences this genome has "
            f"(at least {min_matched:.0%} is required). The genome names its sequences "
            f"{ours}, ...; the resource names them {theirs}, .... They use different "
            "naming, or are for different assemblies"
        )

    beyond = {
        name: (end, sizes[name]) for name, (_, end) in per_contig.items()
        if name in sizes and end > sizes[name]
    }
    if beyond:
        name, (end, length) = sorted(beyond.items())[0]
        reasons.append(
            f"it has positions beyond the end of {len(beyond)} sequence(s): for example "
            f"{name} reaches {end:,}, and in this genome {name} is {length:,} long. The "
            "names match but the coordinates do not, so it is for a different build"
        )
    summary["sequences_with_positions_beyond_their_end"] = len(beyond)

    wrong = {
        name: (length, sizes[name]) for name, length in found["declared_lengths"].items()
        if name in sizes and length != sizes[name]
    }
    if wrong:
        name, (theirs, ours) = sorted(wrong.items())[0]
        reasons.append(
            f"it declares different lengths for {len(wrong)} sequence(s): for example "
            f"{name} as {theirs:,}, which in this genome is {ours:,}"
        )
    summary["declared_lengths_differing"] = len(wrong)

    return reasons, summary


def is_bgzf(path):
    """BGZF is gzip with a 'BC' extra field; tabix needs it, and plain gzip is not it."""
    with open(path, "rb") as fh:
        head = fh.read(18)
    return len(head) >= 14 and head[:4] == b"\x1f\x8b\x08\x04" and head[12:14] == b"BC"


def tabix(path, kind):
    """Index a gzipped VCF or BED in place, recompressing to BGZF first if needed."""
    import pysam

    if not is_bgzf(path):
        plain = path + ".plain"
        with gzip.open(path, "rb") as source, open(plain, "wb") as out:
            shutil.copyfileobj(source, out)
        pysam.tabix_compress(plain, path, force=True)
        os.remove(plain)
    pysam.tabix_index(path, preset=kind, force=True)
    return path + ".tbi"


def main():
    parser = argparse.ArgumentParser(
        description="Check a fetched resource against a genome and prepare it for the vault."
    )
    parser.add_argument("--registry", required=True, help="The resource's file from resources/")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--genome", required=True, help="The vault genome it was requested for")
    parser.add_argument("--fetched", required=True, help="Directory written by fetch_resource.py")
    parser.add_argument("--chrom-sizes", required=True, help="The genome's chrom.sizes from the vault")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    entry = resource_registry.load_entry(args.registry)
    name, kind = entry["name"], entry["kind"]
    label = f"{name}:{args.variant}"

    with open(os.path.join(args.fetched, f"{name}.fetch.json")) as fh:
        record = json.load(fh)
    source = os.path.join(args.fetched, record["file"])

    sizes = read_chrom_sizes(args.chrom_sizes)
    found = inspect(source, kind)
    reasons, summary = gate(found, sizes, float(entry.get("min_matched", 0.95)))

    # Four decimal places: 51 records out of 4.6 million on sequences the genome
    # lacks is 99.9989%, and "100.00%" beside "51 are not" reads as a contradiction.
    print(f"Resource {label} against {args.genome}: {summary['records']:,} records on "
          f"{summary['contigs_in_resource']} sequence(s); "
          f"{summary['share_matched']:.4%} on sequences the genome has.", flush=True)
    if reasons:
        print(f"ERROR: {label} was NOT registered for {args.genome}:", file=sys.stderr)
        for reason in reasons:
            print(f"  - {reason}.", file=sys.stderr)
        print(f"  {entry['variants'][args.variant].get('naming', '')}", file=sys.stderr)
        print("  Request it for a genome it matches, or choose another build of it. "
              "Nothing was added to the vault.", file=sys.stderr, flush=True)
        sys.exit(1)
    if summary["contigs_absent_from_genome"]:
        absent = summary["contigs_absent_from_genome"]
        print(f"       {sum(absent.values()):,} record(s) are on {len(absent)} sequence(s) "
              "this genome does not have; they are listed in the asset's provenance.",
              flush=True)

    os.makedirs(args.outdir, exist_ok=True)
    data = os.path.join(args.outdir, f"{name}.data.gz")
    shutil.copyfile(source, data)
    indexed = False
    if entry.get("index") == "tabix":
        tabix(data, kind)
        indexed = True
        print("       tabix index built.", flush=True)

    provenance = dict(record)
    provenance.update({
        "genome": args.genome,
        "kind": kind,
        "indexed": indexed,
        # Recompressing to BGZF changes the stored bytes but not the content; the
        # checksums above are of the file as published.
        "stored_as_published": is_bgzf(source) or not indexed,
        "naming": entry["variants"][args.variant].get("naming", ""),
        "file_header": found["header"],
        "compatibility": summary,
    })
    with open(os.path.join(args.outdir, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=4)
    with open(os.path.join(args.outdir, "TAG"), "w") as fh:
        fh.write(str(record["version"]) + "\n")


if __name__ == "__main__":
    main()
