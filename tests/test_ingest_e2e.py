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
import csv
import gzip
import hashlib
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


def fake_download(directory, contig_prefix, release, kind, identity=None):
    """
    Write what bin/download.py leaves behind for one assembly.

    `identity` adds the fields that say which request the download answers
    (provider, release_requested), as current downloads record them. Left out, the
    provenance looks like one written before those fields existed.
    """
    os.makedirs(directory)
    with gzip.open(os.path.join(directory, "toyAsm.fa.gz"), "wt") as f:
        f.write(f">{contig_prefix}1\n" + "ACGT" * 10 + "\n")
        f.write(f">{contig_prefix}2\n" + "TTGGCCAA" * 5 + "\n")
    with gzip.open(os.path.join(directory, "toyAsm.annotation.gtf.gz"), "wt") as f:
        f.write(f"{contig_prefix}1\tsrc\tgene\t1\t30\t.\t+\t.\t{GTF_ATTRIBUTES}\n")
        f.write(f"{contig_prefix}1\tsrc\texon\t1\t30\t.\t+\t.\t{EXON_ATTRIBUTES}\n")
    with open(os.path.join(directory, "provenance.json"), "w") as f:
        json.dump({"assembly": "toyAsm", "annotation_release": release,
                   "annotation_release_kind": kind, **(identity or {})}, f)


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
        self.assertEqual(set(self.genome("toyAsm_ensembl")["assets"]), {"fasta", "ensembl_gtf", "provenance"})
        self.assertEqual(set(self.genome("toyAsm_gencode")["assets"]), {"fasta", "gencode_gtf", "provenance"})

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

    # -- the sequence dictionary -------------------------------------------------

    def test_09_dictionary_sits_beside_the_fasta_under_every_alias(self):
        # GATK and Picard, given -R <name>.fa, look for <name>.dict in the same
        # directory. That has to hold whichever alias the user reached the FASTA by.
        for alias in ("toyAsm_ensembl", "toyAsm", "mouse_ens", "toyAsm_gencode", "mouse_gc"):
            fasta = self.seek(f"{alias}/fasta")
            dictionary = self.seek(f"{alias}/fasta.dict")
            self.assertEqual(os.path.dirname(dictionary), os.path.dirname(fasta), alias)
            self.assertEqual(os.path.basename(fasta), f"{alias}.fa")
            self.assertEqual(os.path.basename(dictionary), f"{alias}.dict")
            self.assertTrue(os.path.isfile(fasta + ".fai"), alias)

    def test_10_dictionary_describes_the_sequence_it_sits_beside(self):
        # Checked against values computed here, not against the tool's own output:
        # the samtools on PATH is a pysam-backed shim in the pip environment.
        sequences = {"1": "ACGT" * 10, "2": "TTGGCCAA" * 5}
        for alias, prefix in (("toyAsm_ensembl", ""), ("toyAsm_gencode", "chr")):
            with open(self.seek(f"{alias}/fasta.dict")) as f:
                lines = f.read().splitlines()
            self.assertTrue(lines[0].startswith("@HD"), lines[0])
            records = [dict(field.split(":", 1) for field in line.split("\t")[1:])
                       for line in lines[1:]]
            self.assertEqual([r["SN"] for r in records], [prefix + n for n in sequences])
            for record, sequence in zip(records, sequences.values()):
                self.assertEqual(record["LN"], str(len(sequence)))
                self.assertEqual(record["M5"], hashlib.md5(sequence.encode()).hexdigest())
                # Names the FASTA beside it; never a path into this machine.
                self.assertNotIn("/", record["UR"])

    # -- the download's record ---------------------------------------------------

    def test_11_each_download_record_is_in_the_vault_under_its_release(self):
        asset = self.genome("toyAsm_ensembl")["assets"]["provenance"]
        self.assertEqual(set(asset["tags"]), {"110", "116"})
        self.assertEqual(str(asset["default_tag"]), "116")
        # Reachable by any alias, and named for it.
        for alias, release in (("mouse_ens", "116"), ("toyAsm", "116"), ("mouse_gc", "M33")):
            path = self.seek(f"{alias}/provenance")
            self.assertEqual(os.path.basename(path), f"{alias}.provenance.json")
            with open(path) as f:
                self.assertEqual(json.load(f)["annotation_release"], release)
        with open(self.seek("mouse_ens/provenance:110")) as f:
            self.assertEqual(json.load(f)["annotation_release"], "110")

    def test_12_vault_fasta_is_byte_for_byte_the_downloads_decompressed_content(self):
        # This is what makes the record's genome.content_sha256 checkable: refgenie
        # stores the FASTA decompressed, so `sha256sum` of the vault's file must
        # equal the hash of the text inside the downloaded .fa.gz.
        for alias, prefix in (("toyAsm_ensembl", ""), ("toyAsm_gencode", "chr")):
            text = f">{prefix}1\n" + "ACGT" * 10 + f"\n>{prefix}2\n" + "TTGGCCAA" * 5 + "\n"
            with open(self.seek(f"{alias}/fasta"), "rb") as f:
                self.assertEqual(hashlib.sha256(f.read()).hexdigest(),
                                 hashlib.sha256(text.encode()).hexdigest(), alias)

    def test_13_index_has_a_line_for_every_download_record(self):
        with open(os.path.join(self.outdir, "provenance", "index.tsv"), newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(sorted(r["annotation_release"] for r in rows), ["110", "116", "M33"])
        self.assertEqual({r["digest"] for r in rows}, set(self.genomes()))
        for row in rows:
            aliases = row["aliases"].split(";")
            if row["annotation_release"] == "M33":
                self.assertEqual(sorted(aliases), ["mouse_gc", "toyAsm_gencode"])
                self.assertIn("gencode_gtf:M33", row["assets"].split(";"))
            else:
                self.assertIn("mouse_ens", aliases)
                # Every tag the genome holds, on each of its lines.
                self.assertIn("ensembl_gtf:110,116", row["assets"].split(";"))
            self.assertIn("fasta:default", row["assets"].split(";"))

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

    def test_92_download_records_outlive_the_downloads(self):
        # In the vault itself now, not only in the side directory --cleanup fills.
        for registry_path in ("mouse_ens/provenance", "mouse_ens/provenance:110",
                              "mouse_gc/provenance"):
            self.assertTrue(os.path.isfile(self.seek(registry_path)), registry_path)
        with open(os.path.join(self.outdir, "provenance", "index.tsv"), newline="") as f:
            self.assertEqual(len(list(csv.DictReader(f, delimiter="\t"))), 3)


@unittest.skipIf(missing_requirement(), missing_requirement())
class RegisterDerivedAsset(unittest.TestCase):
    """
    bin/register_asset.py: file something BUILT from a genome against that genome.

    The transcript FASTA here is written by hand; deriving one is DERIVE_TXOME's
    job and is tested where gffread is available (tests/test_nextflow_stub.py).
    """

    REGISTER = os.path.join(REPO, "bin", "register_asset.py")

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="refdata-flow-test-")
        self.addCleanup(shutil.rmtree, self.work, ignore_errors=True)
        self.outdir = os.path.join(self.work, "out")
        self.config_file = os.path.join(self.outdir, "refgenie_export.yaml")
        fake_download(os.path.join(self.outdir, "Ensembl", "Mus_musculus", "toyAsm"),
                      "", "116", "ensembl_release")
        fake_download(os.path.join(self.outdir, "Ensembl", "Mus_musculus", "release-110", "toyAsm"),
                      "", "110", "ensembl_release")
        resolved = os.path.join(self.work, "requests_resolved.csv")
        with open(resolved, "w", newline="") as f:
            f.write("species,assembly,provider,annotation,query,release,aliases\n"
                    "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse\n"
                    "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse\n")
        result = subprocess.run([sys.executable, INGEST, resolved, self.outdir], cwd=self.work,
                                env=ENV, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])

    def derived(self, tag):
        """A transcript FASTA and its provenance record, as DERIVE_TXOME leaves them."""
        transcripts = os.path.join(self.work, f"transcripts-{tag}.fa.gz")
        with gzip.open(transcripts, "wt") as f:
            f.write(f">t1\nACGTACGTAC\n>t2_release_{tag}\nTTGGCCAATT\n")
        provenance = os.path.join(self.work, f"provenance-{tag}.json")
        with open(provenance, "w") as f:
            json.dump({"asset": "fasta_txome", "annotation": f"ensembl_gtf:{tag}"}, f)
        return transcripts, provenance

    def register(self, tag, genome="toyAsm_ensembl", transcripts=None):
        made, provenance = self.derived(tag)
        return subprocess.run(
            [sys.executable, self.REGISTER, self.outdir, "--genome", genome,
             "--asset", "fasta_txome", "--tag", tag,
             "--file", f"transcripts={transcripts or made}", "--file", f"provenance={provenance}",
             "--parent", f"annotation=ensembl_gtf:{tag}"],
            cwd=self.work, env=ENV, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    def txome(self):
        with open(self.config_file) as f:
            genomes = yaml.safe_load(f)["genomes"]
        (digest, genome), = genomes.items()
        return digest, genome["assets"].get("fasta_txome")

    def seek(self, registry_path):
        result = subprocess.run(
            [sys.executable, REFGENIE, "seek", "-c", self.config_file, registry_path],
            cwd=self.work, env=ENV, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip().splitlines()[-1]

    def test_asset_is_built_tagged_and_reachable_under_every_alias(self):
        result = self.register("116")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        digest, asset = self.txome()
        self.assertEqual(set(asset["tags"]), {"116"})
        for alias in ("toyAsm_ensembl", "toyAsm", "mouse"):
            fasta = self.seek(f"{alias}/fasta_txome")
            self.assertTrue(os.path.isfile(fasta), fasta)
            self.assertTrue(os.path.isfile(fasta + ".fai"), alias)
            with open(fasta) as f:
                self.assertEqual(f.readline().strip(), ">t1")
        with open(self.seek("mouse/fasta_txome.provenance")) as f:
            self.assertEqual(json.load(f)["annotation"], "ensembl_gtf:116")

    def test_both_inputs_are_recorded_as_parents(self):
        self.assertEqual(self.register("116").returncode, 0)
        digest, asset = self.txome()
        self.assertEqual(set(asset["tags"]["116"]["asset_parents"]),
                         {f"{digest}/fasta:default", f"{digest}/ensembl_gtf:116"})

    def test_older_release_registered_later_does_not_take_the_default(self):
        self.assertEqual(self.register("116").returncode, 0)
        result = self.register("110")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        _, asset = self.txome()
        self.assertEqual(set(asset["tags"]), {"110", "116"})
        self.assertEqual(str(asset["default_tag"]), "116")
        self.assertIn("/116/", self.seek("mouse/fasta_txome"))
        with open(self.seek("mouse/fasta_txome:110")) as f:
            self.assertIn(">t2_release_110", f.read())

    def test_genome_that_is_not_in_the_vault_is_refused(self):
        result = self.register("116", genome="toyAsm_gencode")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is not in the vault", result.stdout)
        self.assertIsNone(self.txome()[1])

    def test_missing_input_file_is_refused(self):
        result = self.register("116", transcripts=os.path.join(self.work, "nothing.fa.gz"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--file expects", result.stdout)
        self.assertIsNone(self.txome()[1])


@unittest.skipIf(missing_requirement(), missing_requirement())
class AliasRetention(unittest.TestCase):
    """
    An alias, once granted, stays.

    `refgenie build <name>/fasta` resets a genome's alias list to <name> alone
    every time it is called, built or skipped. So each ingest used to leave a
    genome with only the aliases of the request that ingested it last.
    """

    HEADER = "species,assembly,provider,annotation,query,release,aliases\n"

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="refdata-flow-test-")
        self.addCleanup(shutil.rmtree, self.work, ignore_errors=True)
        self.outdir = os.path.join(self.work, "out")
        self.config_file = os.path.join(self.outdir, "refgenie_export.yaml")
        fake_download(os.path.join(self.outdir, "Ensembl", "Mus_musculus", "toyAsm"),
                      "", "116", "ensembl_release")
        fake_download(os.path.join(self.outdir, "Ensembl", "Mus_musculus", "release-110", "toyAsm"),
                      "", "110", "ensembl_release")

    def ingest(self, rows):
        resolved = os.path.join(self.work, "requests_resolved.csv")
        with open(resolved, "w", newline="") as f:
            f.write(self.HEADER + rows)
        result = subprocess.run(
            [sys.executable, INGEST, resolved, self.outdir],
            cwd=self.work, env=ENV, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        with open(self.config_file) as f:
            genomes = yaml.safe_load(f)["genomes"]
        self.assertEqual(len(genomes), 1)
        return set(next(iter(genomes.values()))["aliases"])

    def seek(self, registry_path):
        result = subprocess.run(
            [sys.executable, REFGENIE, "seek", "-c", self.config_file, registry_path],
            cwd=self.work, env=ENV, text=True, capture_output=True,
        )
        return result.returncode, result.stdout.strip().splitlines()[-1:]

    def test_alias_from_an_earlier_request_survives_a_later_one(self):
        self.ingest("Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,first_name\n")
        aliases = self.ingest("Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,second_name\n")
        self.assertTrue({"first_name", "second_name", "toyAsm", "toyAsm_ensembl"} <= aliases)
        # Still usable, not merely listed.
        code, path = self.seek("first_name/fasta")
        self.assertEqual(code, 0)
        self.assertTrue(os.path.isfile(path[0]))

    def test_request_with_no_aliases_removes_none(self):
        self.ingest("Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,first_name\n")
        aliases = self.ingest("Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,\n")
        self.assertIn("first_name", aliases)

    def test_two_rows_for_one_genome_keep_both_rows_aliases(self):
        # Two releases of one assembly are one genome. The second row's build
        # used to wipe the alias the first row had just been given.
        aliases = self.ingest("Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,current_name\n"
                              "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,older_name\n")
        self.assertTrue({"current_name", "older_name"} <= aliases)


@unittest.skipIf(missing_requirement(), missing_requirement())
class OlderVault(unittest.TestCase):
    """
    A vault as earlier versions built it, met by the current ingest.

    Built here with refgenie directly, the old way: the genome named for the
    assembly alone, and the fasta asset from refgenie's stock recipe, so with no
    sequence dictionary.
    """

    RESOLVED = ("species,assembly,provider,annotation,query,release,aliases\n"
                "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens\n")

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="refdata-flow-test-")
        self.addCleanup(shutil.rmtree, self.work, ignore_errors=True)
        self.outdir = os.path.join(self.work, "out")
        self.config_file = os.path.join(self.outdir, "refgenie_export.yaml")
        self.download = os.path.join(self.outdir, "Ensembl", "Mus_musculus", "toyAsm")
        fake_download(self.download, "", "116", "ensembl_release")
        self.resolved = os.path.join(self.work, "requests_resolved.csv")
        with open(self.resolved, "w", newline="") as f:
            f.write(self.RESOLVED)

        self.refgenie("init")
        self.refgenie("build", "toyAsm/fasta",
                      "--files", "fasta=" + os.path.join(self.download, "toyAsm.fa.gz"))
        self.refgenie("build", "toyAsm/ensembl_gtf:115",
                      "--files", "ensembl_gtf=" + os.path.join(self.download, "toyAsm.annotation.gtf.gz"),
                      "--recipe", os.path.join(REPO, "recipes", "ensembl_gtf.json"))
        self.assertEqual(self.fasta_seek_keys("toyAsm") & {"dict"}, set())

    def refgenie(self, *arguments):
        result = subprocess.run(
            [sys.executable, REFGENIE, *arguments, "-c", self.config_file],
            cwd=self.work, env=ENV, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        return result.stdout

    def ingest(self):
        return subprocess.run(
            [sys.executable, INGEST, self.resolved, self.outdir],
            cwd=self.work, env=ENV, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    def genome(self, alias):
        with open(self.config_file) as f:
            genomes = yaml.safe_load(f)["genomes"]
        matches = [g for g in genomes.values() if alias in g["aliases"]]
        self.assertEqual(len(matches), 1, f"'{alias}' should name exactly one genome")
        self.assertEqual(len(genomes), 1, "the ingest must not create a second genome")
        return matches[0]

    def fasta_seek_keys(self, alias):
        return set(self.genome(alias)["assets"]["fasta"]["tags"]["default"]["seek_keys"])

    def test_reingest_gives_it_the_vault_name_and_a_dictionary(self):
        result = self.ingest()
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])

        genome = self.genome("toyAsm_ensembl")
        # Same genome, now reachable by the new name as well as the old ones.
        self.assertTrue({"toyAsm_ensembl", "toyAsm", "mouse_ens"} <= set(genome["aliases"]))
        self.assertIn("dict", self.fasta_seek_keys("toyAsm_ensembl"))
        # Rebuilding the fasta asset must not cost the genome what it already held.
        self.assertTrue({"115", "116"} <= set(genome["assets"]["ensembl_gtf"]["tags"]))

        dictionary = self.refgenie("seek", "toyAsm/fasta.dict").strip().splitlines()[-1]
        self.assertTrue(os.path.isfile(dictionary))
        self.assertEqual(os.path.basename(dictionary), "toyAsm.dict")

    def test_second_reingest_does_not_rebuild_again(self):
        self.assertEqual(self.ingest().returncode, 0)
        result = self.ingest()
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("without a sequence dictionary", result.stdout)

    def test_without_its_download_the_row_is_reported_not_guessed_at(self):
        # No migration path was written for this case: with no source to ingest
        # from, the old entry is never given its vault name, and the ingest says so
        # rather than matching on the bare assembly name another provider may share.
        shutil.rmtree(self.download)
        result = self.ingest()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no genome named 'toyAsm_ensembl'", result.stdout)


@unittest.skipIf(missing_requirement(), missing_requirement())
class StagedSources(unittest.TestCase):
    """
    --source: ingest from directories handed over by name, as REGISTER does, with
    nothing in the published layout under outdir.
    """

    RESOLVED = ("species,assembly,provider,annotation,query,release,aliases\n"
                "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens\n"
                "Mus_musculus,toyAsm,GENCODE,TRUE,toyAsm,,mouse_gc\n")

    def setUp(self):
        self.work = tempfile.mkdtemp(prefix="refdata-flow-test-")
        self.addCleanup(shutil.rmtree, self.work, ignore_errors=True)
        self.outdir = os.path.join(self.work, "out")
        self.resolved = os.path.join(self.work, "requests_resolved.csv")
        with open(self.resolved, "w", newline="") as f:
            f.write(self.RESOLVED)

        # Two directories with the same name, told apart only by their parent -
        # the shape Nextflow stages them in.
        self.ensembl = os.path.join(self.work, "sources", "1", "toyAsm")
        self.gencode = os.path.join(self.work, "sources", "2", "toyAsm")
        fake_download(self.ensembl, "", "116", "ensembl_release",
                      {"provider": "Ensembl", "release_requested": ""})
        fake_download(self.gencode, "chr", "M33", "gencode_release",
                      {"provider": "GENCODE", "release_requested": ""})

    def ingest(self, *extra):
        return subprocess.run(
            [sys.executable, INGEST, self.resolved, self.outdir, *extra],
            cwd=self.work, env=ENV, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    def genomes(self):
        with open(os.path.join(self.outdir, "refgenie_export.yaml")) as f:
            return (yaml.safe_load(f) or {}).get("genomes") or {}

    def names(self):
        return {alias for genome in self.genomes().values() for alias in genome["aliases"]}

    def test_each_source_is_ingested_under_its_own_request(self):
        manifest = os.path.join(self.work, "registered.tsv")
        result = self.ingest("--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}",
                             "--manifest", manifest)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])

        by_name = {alias: genome for genome in self.genomes().values()
                   for alias in genome["aliases"]}
        self.assertEqual(set(by_name["toyAsm_ensembl"]["assets"]), {"fasta", "ensembl_gtf", "provenance"})
        self.assertEqual(set(by_name["toyAsm_gencode"]["assets"]), {"fasta", "gencode_gtf", "provenance"})

        with open(manifest, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual([r["vault_name"] for r in rows], ["toyAsm_ensembl", "toyAsm_gencode"])
        digests = {r["vault_name"]: r["digest"] for r in rows}
        self.assertEqual(set(digests.values()), set(self.genomes()))
        self.assertNotEqual(digests["toyAsm_ensembl"], digests["toyAsm_gencode"])

        # What a derive process is handed: the genome and the annotation release
        # the row means, as files in the vault's data/ tree, named by digest.
        for row, asset, tag, contig in zip(rows, ("ensembl_gtf", "gencode_gtf"),
                                           ("116", "M33"), ("1", "chr1")):
            self.assertEqual((row["annotation_asset"], row["annotation_tag"]), (asset, tag))
            self.assertEqual(row["fasta_txome_tags"], "")
            for path in (row["fasta"], row["fasta"] + ".fai", row["gtf"]):
                self.assertTrue(os.path.isfile(path), path)
                self.assertIn(os.path.join("data", row["digest"]), path)
            with open(row["fasta"]) as f:
                self.assertEqual(f.readline().strip(), ">" + contig)
            with gzip.open(row["gtf"], "rt") as f:
                self.assertEqual(f.readline().split("\t")[0], contig)

    def test_assembly_report_is_kept_as_an_asset_where_the_download_has_one(self):
        report = "# Assembly name:  toyAsm\n1\tassembled-molecule\t1\tChromosome\tCM1.1\t=\tNC_1.1\n"
        with open(os.path.join(self.ensembl, "assembly_report.txt"), "w") as f:
            f.write(report)
        result = self.ingest("--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        by_name = {alias: genome for genome in self.genomes().values()
                   for alias in genome["aliases"]}
        self.assertIn("assembly_report", by_name["toyAsm_ensembl"]["assets"])
        # Not invented for a download that came without one.
        self.assertNotIn("assembly_report", by_name["toyAsm_gencode"]["assets"])
        stored = os.path.join(self.outdir, "alias", "mouse_ens", "assembly_report", "default",
                              "mouse_ens.assembly_report.txt")
        with open(stored) as f:
            self.assertEqual(f.read(), report)

    def test_manifest_leaves_annotation_blank_when_none_was_requested(self):
        with open(self.resolved, "w", newline="") as f:
            f.write(self.RESOLVED.replace("TRUE", "false"))
        manifest = os.path.join(self.work, "registered.tsv")
        result = self.ingest("--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}",
                             "--manifest", manifest)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        with open(manifest, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        for row in rows:
            self.assertTrue(os.path.isfile(row["fasta"]))
            self.assertEqual((row["annotation_asset"], row["annotation_tag"], row["gtf"]),
                             ("", "", ""))

    def test_staged_source_is_never_deleted_by_cleanup(self):
        result = self.ingest("--cleanup",
                             "--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        # They belong to whoever staged them. Provenance is still copied out.
        self.assertTrue(os.path.isfile(os.path.join(self.ensembl, "toyAsm.fa.gz")))
        self.assertTrue(os.path.isfile(os.path.join(self.gencode, "toyAsm.fa.gz")))
        self.assertTrue(os.path.isfile(os.path.join(
            self.outdir, "provenance", "GENCODE", "Mus_musculus", "toyAsm", "provenance.json")))

    def test_swapped_sources_are_refused_and_nothing_is_ingested_from_them(self):
        # The mistake that matters: each key pointing at the other's directory.
        result = self.ingest("--source", f"ensembl/toyAsm@={self.gencode}",
                             "--source", f"gencode/toyAsm@={self.ensembl}")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("is not the download for this request", result.stdout)
        self.assertEqual(self.genomes(), {})

    def test_source_for_the_wrong_release_is_refused(self):
        with open(self.resolved, "w", newline="") as f:
            f.write("species,assembly,provider,annotation,query,release,aliases\n"
                    "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens\n")
        # self.ensembl records release_requested '' - it is the unpinned download.
        result = self.ingest("--source", f"ensembl/toyAsm@110={self.ensembl}")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("release_requested", result.stdout)
        self.assertEqual(self.genomes(), {})

    def test_source_that_matches_no_request_fails_the_run(self):
        result = self.ingest("--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}",
                             "--source", f"ncbi/other@={self.gencode}")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("matches no row", result.stdout)
        # The two that did match are still ingested.
        self.assertTrue({"toyAsm_ensembl", "toyAsm_gencode"} <= self.names())

    def test_named_source_that_is_missing_is_a_failure_not_a_skip(self):
        result = self.ingest("--source", f"ensembl/toyAsm@={self.ensembl}",
                             "--source", f"gencode/toyAsm@={self.gencode}-gone")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not exist", result.stdout)
        self.assertNotIn("toyAsm_gencode", self.names())

    def test_malformed_source_argument_is_rejected(self):
        result = self.ingest("--source", "ensembl/toyAsm@")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--source expects", result.stdout)


if __name__ == "__main__":
    unittest.main()
