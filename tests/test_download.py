# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of the provenance record bin/download.py writes.

Nothing is downloaded. The helpers are given small files written here, and the
one test that runs main() replaces genomepy.install_genome with a function that
leaves behind what genomepy would: a bgzipped FASTA, its sizes file, a GTF and a
README.txt.
"""
import gzip
import hashlib
import io
import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import download  # noqa: E402

FASTA = ">1\n" + "ACGTacgt" * 5 + "\n>2\n" + "TTGGCCAA" * 3 + "\n"

GTF = (
    "#!genome-build toyAsm\n"
    '1\thavana\tgene\t1\t30\t.\t+\t.\tgene_id "G1"; gene_biotype "lncRNA";\n'
    '1\thavana\ttranscript\t1\t30\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_biotype "lncRNA";\n'
    '1\thavana\texon\t1\t30\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_biotype "lncRNA";\n'
    '2\tensembl\tgene\t1\t20\t.\t-\t.\tgene_id "G2"; gene_biotype "protein_coding";\n'
    '2\tensembl\ttranscript\t1\t20\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding";\n'
    '2\tensembl\texon\t1\t20\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_biotype "protein_coding";\n'
    '2\tensembl\ttranscript\t5\t20\t.\t-\t.\tgene_id "G2"; transcript_id "T3"; gene_biotype "protein_coding";\n'
)

# As genomepy 0.16.4 wrote it for a real Ensembl download.
README = """name: toyAsm
provider: Ensembl
original name: toyAsm
original filename: Mus_musculus.toyAsm.dna_sm.primary_assembly.fa.gz
assembly_accession: GCA_000001635.9
tax_id: 10090
mask: soft
genome url: http://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.toyAsm.dna_sm.primary_assembly.fa.gz
annotation url: http://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.toyAsm.116.gtf.gz
sanitized annotation: no changes made.
genomepy version: 0.16.4
date: 2026-11-12 18:01:58
regex: 'alt' (inverted match)

No contigs were removed.
"""


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


class FileTestCase(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = directory.name

    def write_gz(self, name, text):
        path = os.path.join(self.directory, name)
        with gzip.open(path, "wt") as f:
            f.write(text)
        return path


class DescribeFasta(FileTestCase):
    def test_content_hash_is_of_the_sequence_text_not_the_compressed_file(self):
        path = self.write_gz("toyAsm.fa.gz", FASTA)
        record = download.describe_fasta(path)
        self.assertEqual(record["content_sha256"], sha256(FASTA))
        self.assertEqual(record["content_bytes"], len(FASTA))
        with open(path, "rb") as f:
            stored = f.read()
        self.assertEqual(record["file_sha256"], hashlib.sha256(stored).hexdigest())
        self.assertEqual(record["file_bytes"], len(stored))
        self.assertNotEqual(record["file_sha256"], record["content_sha256"])

    def test_recompressing_changes_the_file_hash_and_not_the_content_hash(self):
        first = download.describe_fasta(self.write_gz("a.fa.gz", FASTA))
        path = os.path.join(self.directory, "b.fa.gz")
        with gzip.open(path, "wt", compresslevel=1) as f:
            f.write(FASTA)
        second = download.describe_fasta(path)
        self.assertEqual(first["content_sha256"], second["content_sha256"])
        self.assertNotEqual(first["file_sha256"], second["file_sha256"])

    def test_sequence_count_and_length_come_from_the_sizes_file(self):
        sizes = os.path.join(self.directory, "toyAsm.fa.sizes")
        with open(sizes, "w") as f:
            f.write("1\t40\n2\t24\n")
        record = download.describe_fasta(self.write_gz("toyAsm.fa.gz", FASTA), sizes)
        self.assertEqual((record["sequences"], record["total_length"]), (2, 64))


class DescribeGtf(FileTestCase):
    def setUp(self):
        super().setUp()
        self.record = download.describe_gtf(self.write_gz("toyAsm.annotation.gtf.gz", GTF))

    def test_content_hash_covers_every_line_including_headers(self):
        self.assertEqual(self.record["content_sha256"], sha256(GTF))

    def test_summary_says_what_is_in_the_annotation(self):
        self.assertEqual(self.record["records"], 7)
        self.assertEqual(self.record["feature_counts"],
                         {"exon": 2, "gene": 2, "transcript": 3})
        self.assertEqual((self.record["genes"], self.record["transcripts"]), (2, 3))
        self.assertEqual(self.record["contigs"], ["1", "2"])
        self.assertEqual(self.record["biotypes"], ["lncRNA", "protein_coding"])
        self.assertEqual(self.record["header_lines"], ["#!genome-build toyAsm"])

    def test_gencode_style_gene_type_is_read_as_a_biotype(self):
        record = download.describe_gtf(self.write_gz(
            "g.gtf.gz", 'chr1\tHAVANA\tgene\t1\t9\t.\t+\t.\tgene_id "G"; gene_type "miRNA";\n'))
        self.assertEqual(record["biotypes"], ["miRNA"])

    def test_gtf_with_no_gene_rows_still_counts_genes_and_transcripts(self):
        # UCSC's GTFs have exon and CDS rows only.
        record = download.describe_gtf(self.write_gz(
            "u.gtf.gz",
            'chrI\tncbiRefSeq\texon\t1\t9\t.\t+\t.\tgene_id "A"; transcript_id "NM_1";\n'
            'chrI\tncbiRefSeq\tCDS\t1\t9\t.\t+\t0\tgene_id "A"; transcript_id "NM_1";\n'))
        self.assertEqual((record["genes"], record["transcripts"]), (1, 1))
        self.assertEqual(record["biotypes"], [])


class ParseReadme(FileTestCase):
    def test_key_value_lines_become_fields_and_free_text_is_kept(self):
        path = os.path.join(self.directory, "README.txt")
        with open(path, "w") as f:
            f.write(README)
        fields = download.parse_genomepy_readme(path)
        self.assertEqual(fields["assembly_accession"], "GCA_000001635.9")
        self.assertEqual(fields["tax_id"], "10090")
        self.assertEqual(fields["mask"], "soft")
        self.assertTrue(fields["genome_url"].endswith("primary_assembly.fa.gz"))
        self.assertEqual(fields["sanitized_annotation"], "no changes made.")
        self.assertEqual(fields["regex"], "'alt' (inverted match)")
        self.assertEqual(fields["notes"], ["No contigs were removed."])


class WholeRecord(FileTestCase):
    """main(), with the download itself replaced."""

    def run_download(self, *extra, readme=README, annotation=True, provider="Ensembl"):
        calls = []

        def fake_install(**kwargs):
            calls.append(kwargs)
            # As genomepy does: the directory is named for localname when one is
            # given, and for the name otherwise.
            target = os.path.join(kwargs["genomes_dir"], kwargs.get("localname") or kwargs["name"])
            os.makedirs(target)
            with gzip.open(os.path.join(target, "toyAsm.fa.gz"), "wt") as f:
                f.write(FASTA)
            with open(os.path.join(target, "toyAsm.fa.sizes"), "w") as f:
                f.write("1\t40\n2\t24\n")
            if kwargs["annotation"]:
                with gzip.open(os.path.join(target, "toyAsm.annotation.gtf.gz"), "wt") as f:
                    f.write(GTF)
            if readme is not None:
                with open(os.path.join(target, "README.txt"), "w") as f:
                    f.write(readme)

        argv = ["download.py", "--species", "Mus_musculus", "--assembly", "toyAsm",
                "--provider", provider, "--pipeline-version", "9.9.9", *extra]
        if annotation:
            argv.append("--annotation")
        cwd = os.getcwd()
        os.chdir(self.directory)
        self.addCleanup(os.chdir, cwd)
        code = 0
        with mock.patch.object(download.genomepy, "install_genome", fake_install), \
                mock.patch.object(download, "check_annotation_tools", lambda tag: None), \
                mock.patch.object(download, "_http_last_modified", lambda url, tag: "20260101"), \
                mock.patch.object(sys, "argv", argv), \
                redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()) as err:
            try:
                download.main()
            except SystemExit as e:
                code = e.code or 0
        os.chdir(cwd)
        path = os.path.join(self.directory, "toyAsm", "provenance.json")
        record = None
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                record = json.load(f)
        return code, record, calls, err.getvalue()

    def test_record_identifies_what_was_asked_for_fetched_and_delivered(self):
        code, record, calls, err = self.run_download("--release", "116")
        self.assertEqual(code, 0, err)
        self.assertEqual(record["schema"], download.PROVENANCE_SCHEMA)
        self.assertEqual(record["pipeline_version"], "9.9.9")
        # What was asked for.
        self.assertEqual((record["assembly"], record["provider"], record["release_requested"]),
                         ("toyAsm", "Ensembl", "116"))
        # Where it came from.
        self.assertEqual(record["assembly_accession"], "GCA_000001635.9")
        self.assertIn("/release-116/fasta/", record["genome_url"])
        self.assertIn("/release-116/gtf/", record["annotation_url"])
        self.assertEqual((record["annotation_release"], record["annotation_release_kind"]),
                         ("116", "ensembl_release"))
        # What was delivered.
        self.assertEqual(record["genome"]["content_sha256"], sha256(FASTA))
        self.assertEqual((record["genome"]["sequences"], record["genome"]["total_length"]),
                         (2, 64))
        self.assertEqual(record["annotation"]["content_sha256"], sha256(GTF))
        self.assertEqual(record["annotation"]["transcripts"], 3)
        self.assertEqual(record["genomepy_readme"]["notes"], ["No contigs were removed."])
        self.assertIn("genomepy", record["tool_versions"])
        self.assertNotIn("provenance_errors", record)

    def test_timestamp_is_utc_and_says_so(self):
        _, record, _, _ = self.run_download()
        self.assertTrue(record["retrieved_utc"].endswith("+00:00"), record["retrieved_utc"])
        self.assertNotIn("timestamp", record)

    def test_recorded_install_arguments_are_the_ones_passed_to_genomepy(self):
        _, record, calls, _ = self.run_download()
        self.assertEqual(record["install_arguments"],
                         {"mask": "soft", "keep_alt": False, "invert_match": False, "bgzip": True})
        for key, value in record["install_arguments"].items():
            self.assertEqual(calls[0][key], value, key)

    def test_sequence_only_download_has_no_annotation_block(self):
        code, record, _, err = self.run_download(annotation=False)
        self.assertEqual(code, 0, err)
        self.assertNotIn("annotation", record)
        self.assertIn("genome", record)

    # -- a genome requested by address ---------------------------------------------

    FASTA_URL = "https://example.org/pub/toyAsm.analysis_set.fna.gz"
    GTF_URL = "https://example.org/pub/toyAsm.genes.gtf.gz"

    def test_address_is_what_genomepy_is_asked_for_and_the_name_is_ours(self):
        code, record, calls, err = self.run_download(
            "--url", self.FASTA_URL, "--annotation-url", self.GTF_URL, provider="URL")
        self.assertEqual(code, 0, err)
        call = calls[0]
        self.assertEqual((call["name"], call["localname"], call["provider"]),
                         (self.FASTA_URL, "toyAsm", "URL"))
        self.assertEqual(call["to_annotation"], self.GTF_URL)
        # Files and record are named for the genome, exactly as for a named assembly.
        self.assertEqual(record["assembly"], "toyAsm")
        self.assertEqual((record["url_requested"], record["annotation_url_requested"]),
                         (self.FASTA_URL, self.GTF_URL))
        self.assertIn("genome", record)

    def test_sequence_only_by_address_passes_no_annotation_address(self):
        code, record, calls, err = self.run_download("--url", self.FASTA_URL, provider="URL",
                                                     annotation=False)
        self.assertEqual(code, 0, err)
        self.assertNotIn("to_annotation", calls[0])
        self.assertEqual(record["annotation_url_requested"], "")

    def test_by_address_without_an_address_is_refused_before_anything_is_fetched(self):
        code, record, calls, err = self.run_download(provider="URL", annotation=False)
        self.assertNotEqual(code, 0)
        self.assertEqual(calls, [])
        self.assertIn("needs --url", err)

    def test_annotation_by_address_is_never_left_for_genomepy_to_guess(self):
        code, _, calls, err = self.run_download("--url", self.FASTA_URL, provider="URL")
        self.assertNotEqual(code, 0)
        self.assertEqual(calls, [])
        self.assertIn("--annotation-url", err)

    def test_address_given_for_a_named_provider_is_refused_not_ignored(self):
        code, _, calls, err = self.run_download("--url", self.FASTA_URL, annotation=False)
        self.assertNotEqual(code, 0)
        self.assertEqual(calls, [])
        self.assertIn("only for provider 'url'", err)

    def test_named_provider_record_has_no_address_fields(self):
        _, record, calls, _ = self.run_download(annotation=False)
        self.assertNotIn("url_requested", record)
        self.assertNotIn("localname", calls[0])

    def test_a_part_that_cannot_be_recorded_is_named_and_the_download_still_succeeds(self):
        # No README.txt: the download is good, the record is short, and it says so.
        code, record, _, err = self.run_download(readme=None, annotation=False)
        self.assertEqual(code, 0, err)
        self.assertIn("genome", record)
        self.assertNotIn("genomepy_readme", record)
        self.assertEqual(len(record["provenance_errors"]), 1)
        self.assertIn("README.txt", record["provenance_errors"][0])
        self.assertIn("could not record", err)


if __name__ == "__main__":
    unittest.main()
