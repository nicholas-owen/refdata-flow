# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
End-to-end test of bin/update_refgenie.py against real refgenie, on toy genomes.

Three downloads are faked on disk, in the layout main.nf publishes:

    Ensembl/Mus_musculus/toyAsm/               contigs 1, 2       annotation 116
    GENCODE/Mus_musculus/toyAsm/               contigs chr1, chr2 annotation M33
    Ensembl/Mus_musculus/release-110/toyAsm/   contigs 1, 2       annotation 110

That is one assembly name from two providers, plus an older pinned release of the
first, ingested after the newer one. Each sequence is 40 bases, so the whole run
takes a few seconds and needs no network.

Needs the pipeline's own environment: refgenie beside the interpreter, and
`samtools` on PATH (run_pipeline.sh writes a shim for it into .venv/bin). Skipped
with a reason when either is missing, and on platforms without symlinks.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
INGEST = os.path.join(REPO, "bin", "update_refgenie.py")
BIN_DIR = os.path.dirname(sys.executable)
REFGENIE = os.path.join(BIN_DIR, "refgenie")
ENV = dict(os.environ, PATH=BIN_DIR + os.pathsep + os.environ.get("PATH", ""))

GTF_ATTRIBUTES = ('gene_id "g1"; gene_version "1"; gene_name "G1"; gene_source "s"; '
                  'gene_biotype "protein_coding";')
EXON_ATTRIBUTES = ('gene_id "g1"; gene_version "1"; transcript_id "t1"; '
                   'transcript_version "1"; exon_number "1"; gene_name "G1"; '
                   'gene_source "s"; gene_biotype "protein_coding";')


def missing_requirement():
    if os.name != "posix":
        return "refgenie builds a symlink tree; this test needs a POSIX filesystem"
    if not os.path.isfile(REFGENIE):
        return f"refgenie is not installed beside {sys.executable}"
    if shutil.which("samtools", path=ENV["PATH"]) is None:
        return "samtools is not on PATH (run_pipeline.sh creates the shim in .venv/bin)"
    return None


def fake_download(directory, contig_prefix, release, kind):
    """Write what bin/download.py leaves behind for one assembly."""
    os.makedirs(directory)
    with gzip.open(os.path.join(directory, "toyAsm.fa.gz"), "wt") as f:
        f.write(f">{contig_prefix}1\n" + "ACGT" * 10 + "\n")
        f.write(f">{contig_prefix}2\n" + "TTGGCCAA" * 5 + "\n")
    with gzip.open(os.path.join(directory, "toyAsm.annotation.gtf.gz"), "wt") as f:
        f.write(f"{contig_prefix}1\tsrc\tgene\t1\t30\t.\t+\t.\t{GTF_ATTRIBUTES}\n")
        f.write(f"{contig_prefix}1\tsrc\texon\t1\t30\t.\t+\t.\t{EXON_ATTRIBUTES}\n")
    with open(os.path.join(directory, "provenance.json"), "w") as f:
        json.dump({"assembly": "toyAsm", "annotation_release": release,
                   "annotation_release_kind": kind}, f)


@unittest.skipIf(missing_requirement(), missing_requirement())
class ToyVault(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp(prefix="refdata-flow-test-")
        cls.addClassCleanup(shutil.rmtree, cls.work, ignore_errors=True)
        cls.outdir = os.path.join(cls.work, "out")
        cls.config_file = os.path.join(cls.outdir, "refgenie_export.yaml")

        fake_download(os.path.join(cls.outdir, "Ensembl", "Mus_musculus", "toyAsm"),
                      "", "116", "ensembl_release")
        fake_download(os.path.join(cls.outdir, "GENCODE", "Mus_musculus", "toyAsm"),
                      "chr", "M33", "gencode_release")
        fake_download(os.path.join(cls.outdir, "Ensembl", "Mus_musculus", "release-110", "toyAsm"),
                      "", "110", "ensembl_release")

        cls.resolved = os.path.join(cls.work, "requests_resolved.csv")
        with open(cls.resolved, "w", newline="") as f:
            f.write("species,assembly,provider,annotation,query,release,aliases\n"
                    "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens\n"
                    "Mus_musculus,toyAsm,GENCODE,TRUE,toyAsm,,mouse_gc\n"
                    "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens\n")

        cls.first_run = cls.ingest()

    @classmethod
    def ingest(cls, *extra):
        return subprocess.run(
            [sys.executable, INGEST, cls.resolved, cls.outdir, *extra],
            cwd=cls.work, env=ENV, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    # -- helpers ---------------------------------------------------------------

    def genomes(self):
        with open(self.config_file) as f:
            return yaml.safe_load(f)["genomes"]

    def genome(self, alias):
        matches = [g for g in self.genomes().values() if alias in g["aliases"]]
        self.assertEqual(len(matches), 1, f"'{alias}' should name exactly one genome")
        return matches[0]

    def seek(self, registry_path):
        result = subprocess.run(
            [sys.executable, REFGENIE, "seek", "-c", self.config_file, registry_path],
            cwd=self.work, env=ENV, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip().splitlines()[-1]

    def first_contig(self, registry_path):
        path = self.seek(registry_path)
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            return f.readline().split("\t")[0].lstrip(">").strip()

    # -- the first ingest ------------------------------------------------------

    def test_01_ingest_succeeds(self):
        self.assertEqual(self.first_run.returncode, 0, self.first_run.stdout[-3000:])

    def test_02_two_providers_make_two_genomes(self):
        self.assertEqual(len(self.genomes()), 2)
        self.assertIsNot(self.genome("toyAsm_ensembl"), self.genome("toyAsm_gencode"))

    def test_03_each_genome_keeps_its_own_aliases(self):
        self.assertEqual(set(self.genome("toyAsm_ensembl")["aliases"]),
                         {"toyAsm_ensembl", "toyAsm", "mouse_ens"})
        # The bare name went to the first genome ingested and stays there.
        self.assertEqual(set(self.genome("toyAsm_gencode")["aliases"]),
                         {"toyAsm_gencode", "mouse_gc"})

    def test_04_each_annotation_is_filed_under_its_own_genome(self):
        self.assertEqual(set(self.genome("toyAsm_ensembl")["assets"]), {"fasta", "ensembl_gtf"})
        self.assertEqual(set(self.genome("toyAsm_gencode")["assets"]), {"fasta", "gencode_gtf"})

    def test_05_annotation_contigs_match_the_genome_they_are_filed_under(self):
        # The regression this suite exists for: a chr-prefixed GENCODE GTF used to
        # be attached to the Ensembl genome, whose contigs are 1 and 2.
        self.assertEqual(self.first_contig("toyAsm_gencode/fasta"), "chr1")
        self.assertEqual(self.first_contig("toyAsm_gencode/gencode_gtf"), "chr1")
        self.assertEqual(self.first_contig("toyAsm_ensembl/fasta"), "1")
        self.assertEqual(self.first_contig("toyAsm_ensembl/ensembl_gtf"), "1")

    def test_06_custom_alias_reaches_the_right_genome(self):
        self.assertEqual(self.first_contig("mouse_gc/fasta"), "chr1")
        self.assertEqual(self.first_contig("mouse_ens/fasta"), "1")

    def test_07_both_releases_are_held_as_tags(self):
        asset = self.genome("toyAsm_ensembl")["assets"]["ensembl_gtf"]
        self.assertEqual(set(asset["tags"]), {"110", "116"})

    def test_08_older_release_ingested_later_does_not_take_the_default(self):
        asset = self.genome("toyAsm_ensembl")["assets"]["ensembl_gtf"]
        self.assertEqual(asset["default_tag"], "116")
        self.assertIn("/116/", self.seek("toyAsm_ensembl/ensembl_gtf"))
        self.assertIn("/110/", self.seek("toyAsm_ensembl/ensembl_gtf:110"))

    # -- a second ingest, with cleanup -----------------------------------------
    # Named to sort after everything above: it removes the download directories.

    def test_90_rerun_with_cleanup_changes_nothing_in_the_vault(self):
        before = {a: sorted(g["aliases"]) for g in self.genomes().values() for a in g["aliases"]}
        result = self.ingest("--cleanup")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        after = {a: sorted(g["aliases"]) for g in self.genomes().values() for a in g["aliases"]}
        self.assertEqual(before, after)
        self.assertEqual(
            self.genome("toyAsm_ensembl")["assets"]["ensembl_gtf"]["default_tag"], "116")

        downloads = {
            ("Ensembl", "Mus_musculus", "toyAsm"): "116",
            ("GENCODE", "Mus_musculus", "toyAsm"): "M33",
            ("Ensembl", "Mus_musculus", "release-110", "toyAsm"): "110",
        }
        for parts, release in downloads.items():
            self.assertFalse(os.path.exists(os.path.join(self.outdir, *parts)),
                             f"download not removed: {'/'.join(parts)}")
            # Each release keeps its own record; they must not overwrite each other.
            with open(os.path.join(self.outdir, "provenance", *parts, "provenance.json")) as f:
                self.assertEqual(json.load(f)["annotation_release"], release)

    def test_91_rerun_after_cleanup_still_succeeds(self):
        # The download directories are gone, so only the alias fallback runs.
        result = self.ingest()
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])


if __name__ == "__main__":
    unittest.main()
