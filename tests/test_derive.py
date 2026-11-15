# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of bin/derive_annotation.sh and bin/derive_twobit.sh with the real tools.

No Nextflow and no vault: each script is run on a hand-written annotation or
genome, and what it writes is read back. tests/test_nextflow_stub.py covers the
wiring around them, with the tools stubbed out.

The annotation is written out of order on purpose, has a gene with no name, a
coding transcript on the minus strand with an intron, and a transcript that
appears only on an exon line - so a wrong answer cannot pass by accident.

Needs the tools on PATH, which the pipeline's Python environment does not have:
gtfToGenePred, genePredToBed, bgzip and tabix for the first class, faToTwoBit and
twoBitInfo for the second. Each class is skipped, with the reason, where they are
missing. The pipeline's conda environments have them; see the README.
"""
import gzip
import json
import os
import shutil
import subprocess
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
DERIVE_ANNOTATION = os.path.join(REPO, "bin", "derive_annotation.sh")
DERIVE_TWOBIT = os.path.join(REPO, "bin", "derive_twobit.sh")


def missing(*tools):
    if os.name != "posix":
        return "the derive scripts are bash; this test needs a POSIX system"
    absent = [tool for tool in tools if shutil.which(tool) is None]
    return f"{', '.join(absent)} not on PATH" if absent else None


def attributes(gene, transcript=None, name=None):
    text = f'gene_id "{gene}";'
    if transcript:
        text += f' transcript_id "{transcript}";'
    if name:
        text += f' gene_name "{name}";'
    return text


GTF = "\n".join([
    "#!genome-build toyAsm",
    f"2\tsrc\tgene\t1\t30\t.\t-\t.\t{attributes('g2', name='Beta')}",
    f"2\tsrc\ttranscript\t1\t30\t.\t-\t.\t{attributes('g2', 't2', 'Beta')}",
    f"2\tsrc\texon\t21\t30\t.\t-\t.\t{attributes('g2', 't2', 'Beta')}",
    f"2\tsrc\texon\t1\t10\t.\t-\t.\t{attributes('g2', 't2', 'Beta')}",
    f"2\tsrc\tCDS\t3\t10\t.\t-\t0\t{attributes('g2', 't2', 'Beta')}",
    f"1\tsrc\tgene\t5\t60\t.\t+\t.\t{attributes('g1')}",
    f"1\tsrc\ttranscript\t5\t40\t.\t+\t.\t{attributes('g1', 't1')}",
    f"1\tsrc\texon\t5\t40\t.\t+\t.\t{attributes('g1', 't1')}",
    f"1\tsrc\texon\t50\t60\t.\t+\t.\t{attributes('g1', 't1b')}",
]) + "\n"


class ScriptTestCase(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="refdata-flow-test-")
        self.addCleanup(directory.cleanup)
        self.work = directory.name

    def path(self, *parts):
        return os.path.join(self.work, *parts)

    def run_script(self, script, *arguments, env=None):
        return subprocess.run(["bash", script, *arguments], cwd=self.work, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def lines(self, *parts):
        with open(self.path(*parts)) as f:
            return f.read().splitlines()

    def record(self, *parts):
        with open(self.path(*parts)) as f:
            return json.load(f)


@unittest.skipIf(missing("gtfToGenePred", "genePredToBed", "bgzip", "tabix"),
                 missing("gtfToGenePred", "genePredToBed", "bgzip", "tabix"))
class Annotation(ScriptTestCase):
    def derive(self, gtf=GTF):
        with gzip.open(self.path("toy.gtf.gz"), "wt") as f:
            f.write(gtf)
        return self.run_script(DERIVE_ANNOTATION, "toyAsm_ensembl", "ensembl_gtf", "116",
                               self.path("toy.gtf.gz"), "derived", "1.2.3")

    def test_writes_four_assets_each_with_its_record(self):
        result = self.derive()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(sorted(os.listdir(self.path("derived"))),
                         ["bed12.provenance.json", "gtf_indexed.provenance.json",
                          "refFlat.txt", "refflat.provenance.json", "sorted.gtf.gz",
                          "sorted.gtf.gz.tbi", "transcripts.bed", "tx2gene.provenance.json",
                          "tx2gene.tsv"])
        for asset in ("tx2gene", "gtf_indexed", "bed12", "refflat"):
            record = self.record("derived", f"{asset}.provenance.json")
            self.assertEqual((record["asset"], record["genome"], record["annotation"]),
                             (asset, "toyAsm_ensembl", "ensembl_gtf:116"))
            self.assertTrue(record["tool"] and record["command"], asset)
            self.assertRegex(record["derived_utc"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\+00:00$")
        tx2gene = self.record("derived", "tx2gene.provenance.json")
        self.assertEqual((tx2gene["tool"], tx2gene["tool_version"]), ("refdata-flow", "1.2.3"))
        self.assertEqual((tx2gene["transcripts"], tx2gene["genes"]), (3, 2))

    def test_tx2gene_has_every_transcript_once_in_gtf_order_and_no_empty_name(self):
        self.assertEqual(self.derive().returncode, 0)
        # t1b is named only on an exon line. g1 has no gene_name, so its id stands in.
        self.assertEqual(self.lines("derived", "tx2gene.tsv"),
                         ["t2\tg2\tBeta", "t1\tg1\tg1", "t1b\tg1\tg1"])

    def test_sorted_gtf_keeps_every_record_and_the_providers_order_among_ties(self):
        self.assertEqual(self.derive().returncode, 0)
        with gzip.open(self.path("derived", "sorted.gtf.gz"), "rt") as f:
            written = f.read().splitlines()
        self.assertEqual(sorted(written), sorted(GTF.splitlines()))
        self.assertEqual(written[0], "#!genome-build toyAsm")
        keys = [(line.split("\t")[0], int(line.split("\t")[3]), line.split("\t")[2])
                for line in written[1:]]
        self.assertEqual(keys, [("1", 5, "gene"), ("1", 5, "transcript"), ("1", 5, "exon"),
                                ("1", 50, "exon"),
                                ("2", 1, "gene"), ("2", 1, "transcript"), ("2", 1, "exon"),
                                ("2", 3, "CDS"), ("2", 21, "exon")])

    def test_sorted_gtf_is_bgzipped_and_its_index_answers_a_region(self):
        self.assertEqual(self.derive().returncode, 0)
        with open(self.path("derived", "sorted.gtf.gz"), "rb") as f:
            self.assertEqual(f.read(4), b"\x1f\x8b\x08\x04")       # gzip with the BGZF extra field
        result = subprocess.run(["tabix", self.path("derived", "sorted.gtf.gz"), "2:21-30"],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([line.split("\t")[2] for line in result.stdout.splitlines()],
                         ["gene", "transcript", "exon"])

    def test_bed12_is_one_line_per_transcript_with_no_coding_region_where_there_is_none(self):
        self.assertEqual(self.derive().returncode, 0)
        self.assertEqual(self.lines("derived", "transcripts.bed"), [
            "1\t4\t40\tt1\t0\t+\t40\t40\t0\t1\t36,\t0,",
            "1\t49\t60\tt1b\t0\t+\t60\t60\t0\t1\t11,\t0,",
            "2\t0\t30\tt2\t0\t-\t2\t10\t0\t2\t10,10,\t0,20,",
        ])

    def test_refflat_leads_with_the_gene_id(self):
        self.assertEqual(self.derive().returncode, 0)
        self.assertEqual(self.lines("derived", "refFlat.txt"), [
            "g1\tt1\t1\t+\t4\t40\t40\t40\t1\t4,\t40,",
            "g1\tt1b\t1\t+\t49\t60\t60\t60\t1\t49,\t60,",
            "g2\tt2\t2\t-\t0\t30\t2\t10\t2\t0,20,\t10,30,",
        ])
        self.assertEqual(self.record("derived", "refflat.provenance.json")["gene_column"],
                         "gene_id")

    def test_attribute_that_only_ends_in_gene_id_is_not_taken_for_it(self):
        gtf = ('1\tsrc\texon\t1\t9\t.\t+\t.\tref_gene_id "WRONG"; gene_id "g1"; '
               'transcript_id "t1"; ref_gene_name "WRONG";\n')
        self.assertEqual(self.derive(gtf).returncode, 0)
        self.assertEqual(self.lines("derived", "tx2gene.tsv"), ["t1\tg1\tg1"])

    def test_annotation_with_no_transcripts_is_an_error(self):
        result = self.derive(f"1\tsrc\tgene\t1\t9\t.\t+\t.\t{attributes('g1')}\n")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("carries a transcript_id", result.stdout)

    def test_missing_tool_is_named_before_anything_is_written(self):
        with gzip.open(self.path("toy.gtf.gz"), "wt") as f:
            f.write(GTF)
        empty = self.path("no-tools")
        os.makedirs(empty)
        for tool in ("bash", "gzip"):
            os.symlink(shutil.which(tool), os.path.join(empty, tool))
        result = self.run_script(DERIVE_ANNOTATION, "toyAsm_ensembl", "ensembl_gtf", "116",
                                 self.path("toy.gtf.gz"), "derived", env={"PATH": empty})
        self.assertEqual(result.returncode, 3, result.stdout)
        self.assertIn("is not on PATH", result.stdout)
        self.assertFalse(os.path.exists(self.path("derived")))


@unittest.skipIf(missing("faToTwoBit", "twoBitInfo"), missing("faToTwoBit", "twoBitInfo"))
class TwoBit(ScriptTestCase):
    def derive(self, fai):
        with open(self.path("toy.fa"), "w") as f:
            f.write(">1 a description\nACGTacgtAC\n>2\nNNNNACGT\n")
        with open(self.path("toy.fa.fai"), "w") as f:
            f.write(fai)
        return self.run_script(DERIVE_TWOBIT, "toyAsm_ensembl", self.path("toy.fa"),
                               self.path("toy.fa.fai"), "toy.2bit", "toy.provenance.json")

    def test_twobit_holds_the_genomes_sequences_by_name_and_length(self):
        result = self.derive("1\t10\t18\t10\t11\n2\t8\t32\t8\t9\n")
        self.assertEqual(result.returncode, 0, result.stdout)
        info = subprocess.run(["twoBitInfo", self.path("toy.2bit"), "stdout"],
                              text=True, capture_output=True)
        self.assertEqual(sorted(info.stdout.splitlines()), ["1\t10", "2\t8"])
        record = self.record("toy.provenance.json")
        self.assertEqual((record["asset"], record["tool"], record["sequences"], record["bases"]),
                         ("twobit", "faToTwoBit", 2, 18))
        self.assertEqual(record["soft_masking"], "kept")

    def test_twobit_that_disagrees_with_the_genome_is_refused(self):
        result = self.derive("1\t10\t18\t10\t11\n2\t9\t32\t9\t10\n")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("does not hold the genome's sequences", result.stdout)
        self.assertFalse(os.path.exists(self.path("toy.provenance.json")))


if __name__ == "__main__":
    unittest.main()
