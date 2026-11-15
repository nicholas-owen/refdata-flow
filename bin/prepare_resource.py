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

RENAMING FIRST, WHERE THE REGISTRY SAYS SO

One kind of mismatch is not a refusal but a known translation. dbSNP names its
sequences by RefSeq accession (NC_000001.11), as no genome in a vault does. Its
registry entry therefore names an NCBI assembly report, which lists every
sequence of the build under each of its names - 1, chr1, CM000663.2,
NC_000001.11 - with its length. A sequence is renamed to whichever of its names
this genome uses, and only when the genome's sequence of that name has the
length the report gives; a record on a sequence the genome does not have keeps
the name it was published with. The gate then runs on the result, so a renamed
file has to pass exactly what any other does. Nothing else in a record is touched.

A VEP CACHE

An Ensembl VEP cache is not a list of positions but a directory per sequence,
each holding one file per million bases. The same two questions are asked of it
in the form it allows: what share of the genome's bases lie on sequences the
cache has a directory for, and does any directory hold a region that starts
beyond the end of its sequence. The archive is read to its end to answer them,
which also proves it arrived whole. It is registered still packed; the recipe
unpacks it into the vault.

AFTERWARDS

A resource whose registry entry asks for a tabix index gets one, built here by
pysam's bundled htslib, after re-compressing to BGZF if it was published as plain
gzip. One that asks for none is stored as fetched.

Writes, to --outdir: the file, as <asset>.data.gz (<asset>.data.tar.gz for a VEP
cache); its .tbi if indexed; provenance.json, which is the fetch record plus what
the gate found; and TAG, the version to register it under. A file that is stored
exactly as fetched is hard-linked, not copied, where the filesystem allows: a
30 GB resource should not occupy the work directory twice.
"""
import argparse
import gzip
import json
import os
import re
import shutil
import sys
import tarfile

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


def read_names_report(path):
    """
    An NCBI assembly report as [([names of one sequence], its length)].

    The report has one line per sequence and four columns of names for it:
    Sequence-Name (1), GenBank-Accn (CM000663.2), RefSeq-Accn (NC_000001.11) and
    UCSC-style-name (chr1); 'na' where a scheme has no name for it.
    """
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                continue
            names = [name for name in (fields[0], fields[4], fields[6], fields[9])
                     if name and name != "na"]
            if names:
                rows.append((names, int(fields[8]) if fields[8].isdigit() else None))
    return rows


def name_mapping(rows, sizes):
    """
    {a name the genome does not use: the name it does use for that sequence}.

    A sequence is only ever mapped onto a genome sequence of the length the
    report gives for it: the right name on a different build is not a match, and
    is counted so the record can say so. A name the genome itself uses is never
    mapped away from, whatever the report says about it.

    Returns (mapping, summary).
    """
    mapping, matched, different_length = {}, 0, {}
    for names, length in rows:
        held = [name for name in names if name in sizes]
        target = next((name for name in held if length is None or sizes[name] == length), None)
        if target is None:
            if held:
                different_length[held[0]] = {"genome": sizes[held[0]], "report": length}
            continue
        matched += 1
        for name in names:
            if name != target and name not in sizes:
                mapping.setdefault(name, target)
    return mapping, {
        "sequences_in_report": len(rows),
        "sequences_this_genome_has": matched,
        "same_name_different_length": dict(sorted(different_length.items())),
    }


def rename_line(line, kind, mapping, renamed):
    """One line with its sequence name rewritten, if the mapping has it."""
    if line.startswith("#"):
        if kind == "vcf" and line.startswith("##contig=<ID="):
            start = len("##contig=<ID=")
            ends = [i for i in (line.find(",", start), line.find(">", start)) if i != -1]
            if ends:
                new = mapping.get(line[start:min(ends)])
                if new:
                    return line[:start] + new + line[min(ends):]
        return line
    name, tab, rest = line.partition("\t")
    new = mapping.get(name) if tab else None
    if new is None:
        return line
    renamed[name] = renamed.get(name, 0) + 1
    return new + tab + rest


class BgzfWriter:
    """
    Lines out to a BGZF file, a few megabytes at a time. BGZF, whatever the input
    was, because it is still gzip to anything that reads gzip and is what a tabix
    index needs.
    """

    def __init__(self, path):
        import pysam

        self.out = pysam.BGZFile(path, "wb")
        self.lines, self.held = [], 0

    def write(self, line):
        self.lines.append(line)
        self.held += len(line)
        if self.held >= 4 << 20:
            self.flush()

    def flush(self):
        if self.lines:
            self.out.write("".join(self.lines).encode("utf-8", "surrogateescape"))
            self.lines, self.held = [], 0

    def close(self):
        self.flush()
        self.out.close()


def inspect(path, kind, mapping=None, renamed_to=None):
    """
    One pass over a gzipped VCF or BED. Returns what the gate needs:

        records           number of data lines
        per_contig        {name: [records, furthest end position]}
        declared_lengths  {name: length} from ##contig lines (VCF only)
        header            the VCF's ## lines other than INFO/FILTER/FORMAT/contig,
                          which is where a file says what it is

    A BED interval's end is its third column. A VCF record is taken to end at its
    position: enough to tell a wrong build, without parsing REF or END.

    With a `mapping`, every line is first renamed through it and written to
    `renamed_to`, and it is the renamed file that is described; `renamed` in the
    result is {published name: records renamed}. The file is then read so that
    every byte that is not a sequence name comes out as it went in.
    """
    records, per_contig, declared, header, renamed = 0, {}, {}, [], {}
    writer = BgzfWriter(renamed_to) if mapping is not None else None
    reading = ({"encoding": "utf-8", "errors": "surrogateescape", "newline": "\n"}
               if writer else {"errors": "replace"})
    with gzip.open(path, "rt", **reading) as fh:
        for line in fh:
            if writer:
                line = rename_line(line, kind, mapping, renamed)
                writer.write(line)
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
    if writer:
        writer.close()
    return {"records": records, "per_contig": per_contig,
            "declared_lengths": declared, "header": header, "renamed": renamed}


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


REGION_RE = re.compile(r"^(\d+)-(\d+)")


def inspect_vep_cache(path):
    """
    One pass over a VEP cache tarball, read to its end. Returns:

        files     number of files in it
        root      the directory the cache sits in: <species>/<release>_<assembly>
        contigs   {sequence directory: [files, furthest region start]}
        info      the lines of the cache's own info.txt, where it says what it is
        outside   member names that would be written outside the directory it is
                  unpacked into (absolute, or climbing with ..)

    The layout is <species>/<release>_<assembly>/<sequence>/<start>-<end>.gz, one
    file per million bases, beside a few files such as all_vars.gz that have no
    region in their name.
    """
    files, root, contigs, info, outside = 0, None, {}, [], []
    with tarfile.open(path, "r|gz") as archive:
        for member in archive:
            parts = [part for part in member.name.split("/") if part not in ("", ".")]
            if member.name.startswith("/") or ".." in parts:
                outside.append(member.name)
            if not member.isfile():
                continue
            files += 1
            if len(parts) == 3 and parts[2] == "info.txt":
                info = archive.extractfile(member).read(1 << 16).decode(
                    "utf-8", "replace").splitlines()
            if len(parts) < 4:
                continue
            root = root or "/".join(parts[:2])
            entry = contigs.setdefault(parts[2], [0, 0])
            entry[0] += 1
            region = REGION_RE.match(parts[3])
            if region and int(region.group(1)) > entry[1]:
                entry[1] = int(region.group(1))
    return {"files": files, "root": root, "contigs": contigs, "info": info, "outside": outside}


def gate_vep_cache(found, sizes, min_matched):
    """
    Compare a VEP cache's sequence directories with the genome's sequences.

    Returns (reasons it cannot be registered, summary), as gate() does. The share
    here is of the GENOME: the bases on sequences the cache has a directory for.
    A cache is not expected to cover a sequence with nothing annotated on it, and
    it carries directories for patches and haplotypes a primary assembly leaves
    out, so counting directories either way round would measure the wrong thing.
    """
    reasons = []
    contigs = found["contigs"]
    total = sum(sizes.values())
    covered = sum(length for name, length in sizes.items() if name in contigs)
    share = (covered / total) if total else 0.0
    absent = {name: count for name, (count, _) in contigs.items() if name not in sizes}
    summary = {
        "files": found["files"],
        "sequence_directories": len(contigs),
        "genome_bases": total,
        "genome_bases_on_sequences_in_the_cache": covered,
        "share_matched": round(share, 6),
        "min_matched": min_matched,
        "directories_absent_from_genome": dict(sorted(absent.items())),
    }
    if found["outside"]:
        reasons.append(f"it holds {len(found['outside'])} member(s) that would be written "
                       f"outside the asset, for example '{found['outside'][0]}'")
    if not contigs:
        reasons.append("it contains no sequence directories, so it is not a VEP cache "
                       "as Ensembl publishes them")
        return reasons, summary
    if share < min_matched:
        ours = ", ".join(list(sizes)[:4])
        theirs = ", ".join(sorted(contigs, key=lambda n: -contigs[n][0])[:4])
        reasons.append(
            f"it has a directory for sequences holding only {share:.1%} of this genome's "
            f"bases (at least {min_matched:.0%} is required). The genome names its "
            f"sequences {ours}, ...; the cache names them {theirs}, .... They use "
            "different naming, or are for different assemblies"
        )
    beyond = {
        name: (start, sizes[name]) for name, (_, start) in contigs.items()
        if name in sizes and start > sizes[name]
    }
    if beyond:
        name, (start, length) = sorted(beyond.items())[0]
        reasons.append(
            f"it has regions that start beyond the end of {len(beyond)} sequence(s): for "
            f"example one on {name} starts at {start:,}, and in this genome {name} is "
            f"{length:,} long. The names match but the coordinates do not, so it is for "
            "a different build"
        )
    summary["sequences_with_positions_beyond_their_end"] = len(beyond)
    return reasons, summary


def store(source, destination, will_be_rewritten):
    """
    Put the fetched file where the registration step will find it: a hard link
    where that is possible, a copy where it is not. A file about to be
    recompressed in place is always copied, or the link would carry the rewrite
    back into the fetch step's own output.
    """
    if not will_be_rewritten:
        try:
            os.link(os.path.realpath(source), destination)
            return
        except OSError:
            pass
    shutil.copyfile(source, destination)


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
    variant = entry["variants"][args.variant]
    min_matched = float(entry.get("min_matched", 0.95))
    sizes = read_chrom_sizes(args.chrom_sizes)

    def refuse(reasons):
        shutil.rmtree(args.outdir, ignore_errors=True)
        print(f"ERROR: {label} was NOT registered for {args.genome}:", file=sys.stderr)
        for reason in reasons:
            print(f"  - {reason}.", file=sys.stderr)
        print(f"  {variant.get('naming', '')}", file=sys.stderr)
        print("  Request it for a genome it matches, or choose another build of it. "
              "Nothing was added to the vault.", file=sys.stderr, flush=True)
        sys.exit(1)

    provenance = dict(record)
    provenance.update({"genome": args.genome, "kind": kind,
                       "naming": variant.get("naming", "")})

    if kind == "vep_cache":
        data = os.path.join(args.outdir, f"{name}.data.tar.gz")
        try:
            found = inspect_vep_cache(source)
        except (tarfile.TarError, EOFError, OSError) as e:
            refuse([f"the archive could not be read to its end ({e}); the transfer was "
                    "probably cut short"])
        reasons, summary = gate_vep_cache(found, sizes, min_matched)
        print(f"Resource {label} against {args.genome}: {summary['files']:,} files in "
              f"{summary['sequence_directories']} sequence director(ies); "
              f"{summary['share_matched']:.4%} of the genome's bases are on sequences "
              "the cache covers.", flush=True)
        if reasons:
            refuse(reasons)
        if summary["directories_absent_from_genome"]:
            print(f"       {len(summary['directories_absent_from_genome'])} of its "
                  "directories are for sequences this genome does not have; they are "
                  "listed in the asset's provenance.", flush=True)
        os.makedirs(args.outdir, exist_ok=True)
        store(source, data, will_be_rewritten=False)
        provenance.update({"indexed": False, "stored_as_published": True,
                           "unpacked_in_the_vault": True, "cache_root": found["root"],
                           "file_header": found["info"], "compatibility": summary})
    else:
        data = os.path.join(args.outdir, f"{name}.data.gz")
        wants_index = entry.get("index") == "tabix"
        mapping = None
        if variant.get("rename"):
            report = os.path.join(args.fetched, record["names_report"]["file"])
            rows = read_names_report(report)
            if not rows:
                refuse([f"the table its sequence names are rewritten through "
                        f"({variant['rename']['report']}) lists no sequences"])
            mapping, names_summary = name_mapping(rows, sizes)
            print(f"Resource {label}: renaming sequences to {args.genome}'s own names; "
                  f"the report lists {names_summary['sequences_in_report']:,} sequences, "
                  f"of which this genome has {names_summary['sequences_this_genome_has']:,}.",
                  flush=True)
            # The renamed file is written as it is read, so its directory is needed
            # now; a refusal removes it again.
            os.makedirs(args.outdir, exist_ok=True)
        found = inspect(source, kind, mapping, data if mapping is not None else None)
        reasons, summary = gate(found, sizes, min_matched)

        # Four decimal places: 51 records out of 4.6 million on sequences the genome
        # lacks is 99.9989%, and "100.00%" beside "51 are not" reads as a contradiction.
        print(f"Resource {label} against {args.genome}: {summary['records']:,} records on "
              f"{summary['contigs_in_resource']} sequence(s); "
              f"{summary['share_matched']:.4%} on sequences the genome has.", flush=True)
        if reasons:
            refuse(reasons)
        if summary["contigs_absent_from_genome"]:
            absent = summary["contigs_absent_from_genome"]
            print(f"       {sum(absent.values()):,} record(s) are on {len(absent)} sequence(s) "
                  "this genome does not have; they are listed in the asset's provenance.",
                  flush=True)

        os.makedirs(args.outdir, exist_ok=True)
        if mapping is None:
            store(source, data, will_be_rewritten=wants_index and not is_bgzf(source))
        else:
            names_summary.update({
                "report": variant["rename"]["report"],
                "sequences_renamed": len(found["renamed"]),
                "records_renamed": sum(found["renamed"].values()),
                "examples": {old: mapping[old] for old in sorted(found["renamed"])[:5]},
            })
            provenance["renamed"] = names_summary
            print(f"       {names_summary['records_renamed']:,} record(s) on "
                  f"{names_summary['sequences_renamed']} sequence(s) were renamed.", flush=True)
        indexed = False
        if wants_index:
            tabix(data, kind)
            indexed = True
            print("       tabix index built.", flush=True)

        provenance.update({
            "indexed": indexed,
            # Recompressing to BGZF changes the stored bytes but not the content; the
            # checksums above are of the file as published. A renamed file differs
            # from the published one in content too, by exactly the names.
            "stored_as_published": mapping is None and (is_bgzf(source) or not indexed),
            "file_header": found["header"],
            "compatibility": summary,
        })

    with open(os.path.join(args.outdir, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=4)
    with open(os.path.join(args.outdir, "TAG"), "w") as fh:
        fh.write(str(record["version"]) + "\n")


if __name__ == "__main__":
    main()
