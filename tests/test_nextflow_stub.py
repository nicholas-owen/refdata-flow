# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests that execute main.nf. No network, no conda, no genomepy.

Most of them use `nextflow run -stub`, which swaps FETCH_GENOME's script for its
stub block: a toy assembly directory in the layout bin/download.py leaves behind.
Everything around the download is real - request validation, the duplicate guard,
the publish path - and the last stub test hands the output to
bin/update_refgenie.py, so the two stages are checked against each other rather
than against a layout each one assumes.

The retry tests run the real script block instead, with a fake `python3` first on
PATH that does nothing but exit with a chosen code. bin/download.py is launched
through `#!/usr/bin/env python3`, so that stands in for a download failing.

Needs `nextflow` on PATH. Skipped with a reason when it is not; each run costs a
JVM start, so this file takes a minute or so.
"""
import gzip
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
MAIN_NF = os.path.join(REPO, "main.nf")
INGEST = os.path.join(REPO, "bin", "update_refgenie.py")
BIN_DIR = os.path.dirname(sys.executable)
REFGENIE = os.path.join(BIN_DIR, "refgenie")

HEADER = "species,assembly,provider,annotation,query,release,aliases\n"
REQUESTS = (
    HEADER
    + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens\n"
    + "Mus_musculus,toyAsm,GENCODE,TRUE,toyAsm,,mouse_gc\n"
    + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens\n"
    + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n"
)


def missing_requirement():
    if os.name != "posix":
        return "Nextflow tasks are bash scripts; this test needs a POSIX system"
    if shutil.which("nextflow") is None:
        return "nextflow is not on PATH"
    return None


def ingest_unavailable():
    if not os.path.isfile(REFGENIE):
        return f"refgenie is not installed beside {sys.executable}"
    if shutil.which("samtools", path=BIN_DIR + os.pathsep + os.environ.get("PATH", "")) is None:
        return "samtools is not on PATH (run_pipeline.sh creates the shim in .venv/bin)"
    return None


@unittest.skipIf(missing_requirement(), missing_requirement())
class NextflowTestCase(unittest.TestCase):
    """One throwaway launch directory per class; work/ and .nextflow/ stay inside it."""

    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.mkdtemp(prefix="refdata-flow-nf-")
        cls.addClassCleanup(shutil.rmtree, cls.work, ignore_errors=True)

    @classmethod
    def write_requests(cls, name, content):
        path = os.path.join(cls.work, name)
        with open(path, "w", newline="") as f:
            f.write(content)
        return path

    @classmethod
    def nextflow(cls, requests, outdir, *extra, stub=True, env=None,
                 register=False, venv_bin=BIN_DIR, default_txome=False, default_derived=False):
        # Most tests here are about the download side and pass --register false, so
        # they run wherever nextflow does. The Registration tests turn it on and
        # point REGISTER at the interpreter running this suite, which is the
        # pipeline's own environment when the suite is run as the README says.
        command = ["nextflow", "run", MAIN_NF, "-profile", "standard", "-ansi-log", "false",
                   "--input", requests, "--outdir", outdir]
        # --check_releases false: the preflight would otherwise ask Ensembl what its
        # current release is, for an assembly that exists only in these tests.
        command += (["--venv_bin", venv_bin, "--check_releases", "false"] if register
                    else ["--register", "false"])
        # Transcriptomes are on by default. Most tests here are about something
        # else and count tasks or list a genome's assets exactly, so they run with
        # it off unless they say otherwise; the Transcriptome tests cover it, and
        # one of them covers the default itself.
        if register and "--txome" not in extra and not default_txome:
            command += ["--txome", "false"]
        # The same goes for the small derived assets (tx2gene and the rest): on by
        # default, off here unless a test asks, and covered by DerivedAssets.
        if register and "--derived" not in extra and not default_derived:
            command += ["--derived", "false"]
        command += extra
        if stub:
            command.append("-stub")
        return subprocess.run(
            command, cwd=cls.work, env=env, text=True, timeout=600,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )


class StubRun(NextflowTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", REQUESTS)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir)

    def published(self, *parts):
        return os.path.join(self.outdir, *parts)

    def test_01_run_succeeds(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])

    def test_02_unpinned_request_keeps_the_original_layout(self):
        for name in ("toyAsm.fa.gz", "toyAsm.annotation.gtf.gz", "provenance.json"):
            self.assertTrue(os.path.isfile(
                self.published("Ensembl", "Mus_musculus", "toyAsm", name)), name)

    def test_03_pinned_release_is_published_one_level_down(self):
        directory = self.published("Ensembl", "Mus_musculus", "release-110", "toyAsm")
        self.assertTrue(os.path.isfile(os.path.join(directory, "toyAsm.fa.gz")))
        # Its own record, not the unpinned row's: the two used to share a directory.
        with open(os.path.join(directory, "provenance.json")) as f:
            self.assertEqual(json.load(f)["annotation_release"], "110")
        with open(self.published("Ensembl", "Mus_musculus", "toyAsm", "provenance.json")) as f:
            self.assertNotEqual(json.load(f)["annotation_release"], "110")

    def test_04_sequence_only_request_publishes_no_annotation(self):
        directory = self.published("NCBI", "Foo_bar", "seqOnly")
        self.assertTrue(os.path.isfile(os.path.join(directory, "seqOnly.fa.gz")))
        self.assertFalse(os.path.exists(os.path.join(directory, "seqOnly.annotation.gtf.gz")))

    def test_05_stub_output_cannot_pass_for_a_real_download(self):
        with open(self.published("GENCODE", "Mus_musculus", "toyAsm", "provenance.json")) as f:
            self.assertIs(json.load(f)["stub"], True)

    def test_06_no_registration_was_attempted(self):
        self.assertNotIn("REGISTER", self.run_result.stdout)
        self.assertFalse(os.path.exists(self.published("refgenie_export.yaml")))

    @unittest.skipIf(ingest_unavailable(), ingest_unavailable())
    def test_07_ingest_run_by_hand_finds_everything_the_workflow_published(self):
        # bin/update_refgenie.py can still be run on its own against the published
        # layout, with no --source. main.nf and the script each build that path
        # independently; here one's output is the other's input.
        env = dict(os.environ, PATH=BIN_DIR + os.pathsep + os.environ.get("PATH", ""))
        result = subprocess.run(
            [sys.executable, INGEST, self.requests, self.outdir],
            cwd=self.work, env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("not found - nothing to ingest", result.stdout)

        with open(os.path.join(self.outdir, "refgenie_export.yaml")) as f:
            genomes = yaml.safe_load(f)["genomes"]
        by_name = {alias: genome for genome in genomes.values() for alias in genome["aliases"]}

        self.assertEqual(len(genomes), 3)
        self.assertIsNot(by_name["toyAsm_ensembl"], by_name["toyAsm_gencode"])
        self.assertEqual(set(by_name["toyAsm_gencode"]["assets"]), {"fasta", "gencode_gtf", "provenance"})
        self.assertEqual(set(by_name["seqOnly_ncbi"]["assets"]), {"fasta", "provenance"})
        # Both Ensembl rows reached the vault, each under its own tag.
        self.assertEqual(set(by_name["toyAsm_ensembl"]["assets"]["ensembl_gtf"]["tags"]),
                         {"999", "110"})


def load_vault(outdir):
    """{alias: genome entry} for every alias in the vault's config."""
    with open(os.path.join(outdir, "refgenie_export.yaml")) as f:
        genomes = yaml.safe_load(f)["genomes"]
    return {alias: genome for genome in genomes.values() for alias in genome["aliases"]}


@unittest.skipIf(ingest_unavailable(), ingest_unavailable())
class Registration(NextflowTestCase):
    """Download and catalogue in one Nextflow run, as run_pipeline.sh now does it."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", REQUESTS)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir, register=True)

    def published(self, *parts):
        return os.path.join(self.outdir, *parts)

    def assert_vault_complete(self, outdir):
        vault = load_vault(outdir)
        self.assertEqual(len({id(genome) for genome in vault.values()}), 3)
        self.assertIsNot(vault["toyAsm_ensembl"], vault["toyAsm_gencode"])
        self.assertEqual(set(vault["toyAsm_gencode"]["assets"]), {"fasta", "gencode_gtf", "provenance"})
        self.assertEqual(set(vault["seqOnly_ncbi"]["assets"]), {"fasta", "provenance"})
        # The pinned and the unpinned Ensembl download share an assembly name and a
        # directory name. Both tags present means each was matched to its own row.
        self.assertEqual(set(vault["toyAsm_ensembl"]["assets"]["ensembl_gtf"]["tags"]),
                         {"999", "110"})
        # Built by recipes/fasta.json, so each genome has its sequence dictionary,
        # beside its FASTA under the alias name. The file is checked, not the
        # config: refgenie can list a seek key for a file it never wrote.
        for name in ("toyAsm_ensembl", "toyAsm_gencode", "seqOnly_ncbi"):
            directory = os.path.join(outdir, "alias", name, "fasta", "default")
            for filename in (f"{name}.fa", f"{name}.fa.fai", f"{name}.dict"):
                self.assertTrue(os.path.isfile(os.path.join(directory, filename)), filename)
        return vault

    def test_01_run_succeeds_with_one_register_task(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])
        self.assertEqual(self.run_result.stdout.count("Submitted process > REGISTER"), 1)
        # This class runs with --txome false, which must leave them out entirely.
        self.assertNotIn("DERIVE_TXOME", self.run_result.stdout)

    def test_02_every_download_is_in_the_vault(self):
        self.assert_vault_complete(self.outdir)

    def test_03_bare_assembly_name_goes_to_the_first_row_of_the_request(self):
        # Ensembl is listed first. Which download finished first must not matter.
        vault = load_vault(self.outdir)
        self.assertIs(vault["toyAsm"], vault["toyAsm_ensembl"])
        self.assertEqual(set(vault["toyAsm_gencode"]["aliases"]), {"toyAsm_gencode", "mouse_gc"})

    def test_04_sources_are_kept_by_default(self):
        self.assertTrue(os.path.isfile(
            self.published("Ensembl", "Mus_musculus", "release-110", "toyAsm", "toyAsm.fa.gz")))
        self.assertTrue(os.path.isfile(
            self.published("GENCODE", "Mus_musculus", "toyAsm", "toyAsm.fa.gz")))

    def test_05_second_run_downloads_nothing_and_changes_nothing(self):
        # No -resume: Nextflow's cache is not what prevents the downloads. The
        # preflight finds all four requests in the vault, so FETCH_GENOME is never
        # even created - which is the behaviour an HPC run without work/ needs.
        before = {alias: sorted(genome["assets"]) for alias, genome in load_vault(self.outdir).items()}
        result = self.nextflow(self.requests, self.outdir, register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("FETCH_GENOME", result.stdout)
        self.assertEqual(result.stdout.count("  SKIP "), 4)
        self.assertIn("0 to download, 4 already in the vault", result.stdout)
        # REGISTER still runs, with nothing staged: aliases may have been added.
        self.assertIn("Submitted process > REGISTER (0 downloaded)", result.stdout)
        after = {alias: sorted(genome["assets"]) for alias, genome in load_vault(self.outdir).items()}
        self.assertEqual(before, after)

    def test_05b_preflight_off_downloads_everything_again(self):
        # All four tasks run, where the test above had none - although every one
        # of these genomes is in the vault.
        result = self.nextflow(self.requests, self.outdir, "--preflight", "false", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("PREFLIGHT", result.stdout)
        self.assertEqual(result.stdout.count("Submitted process > FETCH_GENOME"), 4)
        self.assertIn("Submitted process > REGISTER (4 downloaded)", result.stdout)

    def test_05c_only_the_new_request_is_downloaded_and_new_aliases_still_apply(self):
        outdir = os.path.join(self.work, "grow.out")
        first = self.write_requests("grow1.csv", HEADER
                                    + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens\n"
                                    + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n")
        self.assertEqual(self.nextflow(first, outdir, register=True).returncode, 0)

        # Same two genomes, one of them with an extra alias, plus two new requests:
        # an older release of a genome that is there, and a genome that is not.
        second = self.write_requests("grow2.csv", HEADER
                                     + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens;added_later\n"
                                     + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n"
                                     + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens\n"
                                     + "Mus_musculus,toyAsm,GENCODE,TRUE,toyAsm,,mouse_gc\n")
        result = self.nextflow(second, outdir, register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertEqual(result.stdout.count("Submitted process > FETCH_GENOME"), 2)
        self.assertIn("2 to download, 2 already in the vault", result.stdout)

        vault = self.assert_vault_complete(outdir)
        # Applied although nothing was downloaded for that row - and not lost when
        # the pinned row, which shares its genome, was ingested after it.
        self.assertIn("added_later", vault["toyAsm_ensembl"]["aliases"])

    def test_05d_genome_requested_by_address_is_fetched_and_registered(self):
        # Under -stub nothing is fetched from the address; what is tested is that
        # such a row travels the whole workflow like any other: validated,
        # published under its provider, and filed under its own name.
        outdir = os.path.join(self.work, "byurl.out")
        requests = self.write_requests(
            "byurl.csv", HEADER.rstrip("\n") + ",url,annotation_url\n"
            + "Homo_sapiens,noalt,URL,TRUE,noalt,,hs_noalt,"
              "https://example.org/pub/noalt.fna.gz,https://example.org/pub/noalt.gtf.gz\n"
            + "Foo_bar,seqOnly,NCBI,false,seqOnly,,,,\n")
        result = self.nextflow(requests, outdir, register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertTrue(os.path.isfile(
            os.path.join(outdir, "URL", "Homo_sapiens", "noalt", "noalt.fa.gz")))
        vault = load_vault(outdir)
        self.assertIs(vault["noalt_url"], vault["hs_noalt"])
        # An annotation given by address is filed as exactly that.
        self.assertEqual(set(vault["noalt_url"]["assets"]), {"fasta", "url_gtf", "provenance"})

    def test_06_cleanup_false_on_the_command_line_still_keeps_sources(self):
        # A command-line value arrives as text, and the text 'false' is true to
        # Groovy. Read naively, this switched cleanup ON.
        outdir = os.path.join(self.work, "cleanup-false.out")
        result = self.nextflow(self.requests, outdir, "--cleanup", "false", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertTrue(os.path.isfile(
            os.path.join(outdir, "Ensembl", "Mus_musculus", "toyAsm", "toyAsm.fa.gz")))

    def test_07_cleanup_keeps_no_sources_but_keeps_their_provenance(self):
        outdir = os.path.join(self.work, "cleanup-true.out")
        result = self.nextflow(self.requests, outdir, "--cleanup", "true", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assert_vault_complete(outdir)
        for provider in ("Ensembl", "GENCODE", "NCBI"):
            self.assertFalse(os.path.exists(os.path.join(outdir, provider)), provider)
        with open(os.path.join(outdir, "provenance", "Ensembl", "Mus_musculus",
                               "release-110", "toyAsm", "provenance.json")) as f:
            self.assertEqual(json.load(f)["release_requested"], "110")
        with open(os.path.join(outdir, "provenance", "Ensembl", "Mus_musculus",
                               "toyAsm", "provenance.json")) as f:
            self.assertEqual(json.load(f)["release_requested"], "")

    def test_90_cleanup_removes_sources_an_earlier_run_kept(self):
        # Sorted last: it empties the published tree the tests above look at.
        self.assertTrue(os.path.isdir(self.published("Ensembl")))
        result = self.nextflow(self.requests, self.outdir, "-resume", "--cleanup", "true",
                               register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        for parts in (("Ensembl", "Mus_musculus", "toyAsm"),
                      ("Ensembl", "Mus_musculus", "release-110", "toyAsm"),
                      ("GENCODE", "Mus_musculus", "toyAsm"),
                      ("NCBI", "Foo_bar", "seqOnly")):
            self.assertFalse(os.path.exists(self.published(*parts)), "/".join(parts))
            self.assertTrue(os.path.isfile(
                self.published("provenance", *parts, "provenance.json")), "/".join(parts))
        self.assert_vault_complete(self.outdir)


def read_fasta(path):
    """{record name: sequence} from a FASTA file, plain or gzipped."""
    import gzip
    records, name = {}, None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = ""
            elif name is not None:
                records[name] += line
    return records


@unittest.skipIf(ingest_unavailable(), ingest_unavailable())
class Transcriptome(NextflowTestCase):
    """
    The path from REGISTER's manifest, through DERIVE_TXOME, into the vault.

    The class's main run passes no --txome flag at all: a transcriptome is part of
    every run by default, and this is the test that it is.

    Under -stub the derivation itself is a stand-in (one short record per
    transcript in the GTF). Everything around it is real: which rows get a
    transcriptome, the staging of vault files into the task, and the registration.
    RealDerivation below runs gffread.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", REQUESTS)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir, register=True,
                                      default_txome=True)

    def refgenie(self, *arguments):
        env = dict(os.environ, PATH=BIN_DIR + os.pathsep + os.environ.get("PATH", ""))
        result = subprocess.run(
            [sys.executable, REFGENIE, *arguments, "-c",
             os.path.join(self.outdir, "refgenie_export.yaml")],
            cwd=self.work, env=env, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.strip().splitlines()[-1]

    def test_01_run_succeeds(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])

    def test_02_one_transcriptome_per_genome_and_release(self):
        derived = sorted(line.split("DERIVE_TXOME (")[1].rstrip(")")
                         for line in self.run_result.stdout.splitlines()
                         if "Submitted process > DERIVE_TXOME" in line)
        # Both Ensembl releases and the GENCODE genome. Not the sequence-only
        # genome: with no annotation there is nothing to derive from.
        self.assertEqual(derived, ["toyAsm_ensembl:110", "toyAsm_ensembl:999",
                                   "toyAsm_gencode:999"])

    def test_03_tagged_with_the_annotation_release_newest_as_default(self):
        vault = load_vault(self.outdir)
        txome = vault["toyAsm_ensembl"]["assets"]["fasta_txome"]
        self.assertEqual(set(txome["tags"]), {"110", "999"})
        self.assertEqual(str(txome["default_tag"]), "999")
        self.assertNotIn("fasta_txome", vault["seqOnly_ncbi"]["assets"])

    def test_04_each_release_records_the_gtf_it_came_from(self):
        vault = load_vault(self.outdir)
        with open(os.path.join(self.outdir, "refgenie_export.yaml")) as f:
            genomes = yaml.safe_load(f)["genomes"]
        digest = next(d for d, g in genomes.items() if "toyAsm_ensembl" in g["aliases"])
        tags = vault["toyAsm_ensembl"]["assets"]["fasta_txome"]["tags"]
        for release in ("110", "999"):
            self.assertEqual(set(tags[release]["asset_parents"]),
                             {f"{digest}/fasta:default", f"{digest}/ensembl_gtf:{release}"})

    def test_05_reachable_by_alias_with_its_index_and_provenance(self):
        fasta = self.refgenie("seek", "mouse_gc/fasta_txome")
        self.assertEqual(os.path.basename(fasta), "mouse_gc.fa")
        self.assertEqual(list(read_fasta(fasta)), ["t1"])
        self.assertTrue(os.path.isfile(fasta + ".fai"))
        with open(self.refgenie("seek", "mouse_gc/fasta_txome.provenance")) as f:
            self.assertEqual(json.load(f)["annotation"], "gencode_gtf:999")
        # Not to be confused with the genome, which has the same file name one
        # directory over.
        self.assertNotEqual(fasta, self.refgenie("seek", "mouse_gc/fasta"))

    def test_06_second_run_derives_nothing(self):
        result = self.nextflow(self.requests, self.outdir, "--txome", "true", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("DERIVE_TXOME", result.stdout)
        self.assertNotIn("REGISTER_ASSET", result.stdout)

    def test_07_transcriptome_can_be_added_to_genomes_already_in_the_vault(self):
        # A vault built with --txome false, or before transcriptomes existed, and
        # then run normally: nothing is downloaded, and every annotated genome
        # gets one.
        outdir = os.path.join(self.work, "later.out")
        without = self.nextflow(self.requests, outdir, "--txome", "false", register=True)
        self.assertEqual(without.returncode, 0, without.stdout[-3000:])
        self.assertNotIn("DERIVE_TXOME", without.stdout)
        self.assertNotIn("fasta_txome", load_vault(outdir)["toyAsm_gencode"]["assets"])
        result = self.nextflow(self.requests, outdir, "--txome", "true", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("FETCH_GENOME", result.stdout)
        self.assertEqual(result.stdout.count("Submitted process > DERIVE_TXOME"), 3)
        self.assertIn("fasta_txome", load_vault(outdir)["toyAsm_gencode"]["assets"])


@unittest.skipIf(ingest_unavailable(), ingest_unavailable())
class Resources(NextflowTestCase):
    """
    The `resources` column: fetch, check against the genome, register.

    Under -stub the fetch writes a two-line file on the genome's first sequence.
    The compatibility check and the registration that follow are the real ones.
    """

    HEADER = HEADER.rstrip("\n") + ",url,annotation_url,resources\n"
    REQUESTS = (HEADER
                + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens,,,clinvar:GRCh38;blacklist:hg38\n"
                + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens,,,clinvar:GRCh38\n"
                + "Foo_bar,seqOnly,NCBI,false,seqOnly,,,,,\n")

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", cls.REQUESTS)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir, register=True)

    def submitted(self, process, output=None):
        return sorted(line.split(f"{process} (")[1].rstrip(")")
                      for line in (output or self.run_result.stdout).splitlines()
                      if f"Submitted process > {process} (" in line)

    def test_01_run_succeeds(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])

    def test_02_each_resource_is_fetched_once_per_genome(self):
        # Two rows are the same genome and both ask for clinvar: one fetch, not two.
        # The row that asks for nothing gets nothing.
        expected = ["blacklist:hg38 for toyAsm_ensembl", "clinvar:GRCh38 for toyAsm_ensembl"]
        self.assertEqual(self.submitted("FETCH_RESOURCE"), expected)
        self.assertEqual(self.submitted("PREPARE_RESOURCE"), expected)
        self.assertEqual(self.submitted("REGISTER_ASSET"),
                         ["toyAsm_ensembl/blacklist:stub1", "toyAsm_ensembl/clinvar:stub1"])

    def test_03_registered_with_index_provenance_and_the_genome_as_parent(self):
        vault = load_vault(self.outdir)
        genome = vault["toyAsm_ensembl"]
        self.assertNotIn("clinvar", vault["seqOnly_ncbi"]["assets"])
        clinvar = genome["assets"]["clinvar"]["tags"]["stub1"]
        self.assertEqual(set(clinvar["seek_keys"]) - {"dir"}, {"clinvar", "tbi", "provenance"})
        self.assertEqual([p.split("/")[1] for p in clinvar["asset_parents"]], ["fasta:default"])
        directory = os.path.join(self.outdir, "alias", "mouse_ens", "clinvar", "stub1")
        for name in ("mouse_ens.vcf.gz", "mouse_ens.vcf.gz.tbi", "mouse_ens.provenance.json"):
            self.assertTrue(os.path.isfile(os.path.join(directory, name)), name)
        # The blacklist asks for no index and gets none.
        self.assertEqual(sorted(os.listdir(
            os.path.join(self.outdir, "alias", "mouse_ens", "blacklist", "stub1"))),
            ["_refgenie_build", "mouse_ens.bed.gz", "mouse_ens.provenance.json"])

    def test_04_provenance_says_what_the_compatibility_check_found(self):
        with open(os.path.join(self.outdir, "alias", "mouse_ens", "clinvar", "stub1",
                               "mouse_ens.provenance.json")) as f:
            record = json.load(f)
        self.assertEqual((record["resource"], record["variant"], record["genome"]),
                         ("clinvar", "GRCh38", "toyAsm_ensembl"))
        self.assertEqual(record["compatibility"]["share_matched"], 1.0)
        self.assertTrue(record["indexed"])

    def test_05_second_run_fetches_and_registers_nothing(self):
        result = self.nextflow(self.requests, self.outdir, register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        # The fetch step is asked, and answers that the vault is up to date.
        self.assertEqual(len(self.submitted("FETCH_RESOURCE", result.stdout)), 2)
        self.assertEqual(result.stdout.count("  SKIP   clinvar:GRCh38"), 1)
        self.assertNotIn("PREPARE_RESOURCE", result.stdout)
        self.assertNotIn("REGISTER_ASSET", result.stdout)

    def test_06_resources_and_transcriptomes_register_in_the_same_run(self):
        outdir = os.path.join(self.work, "both.out")
        result = self.nextflow(self.requests, outdir, register=True, default_txome=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertEqual(self.submitted("REGISTER_ASSET", result.stdout),
                         ["toyAsm_ensembl/blacklist:stub1", "toyAsm_ensembl/clinvar:stub1",
                          "toyAsm_ensembl/fasta_txome:110", "toyAsm_ensembl/fasta_txome:999"])

        # The index and the RO-Crate describe the vault as the run LEFT it: with
        # what REGISTER_ASSET added, not only what REGISTER wrote before that.
        with open(os.path.join(outdir, "provenance", "index.tsv")) as f:
            index = f.read()
        for held in ("blacklist:stub1", "clinvar:stub1", "fasta_txome:110,999"):
            self.assertIn(held, index)
        with open(os.path.join(outdir, "ro-crate-metadata.json")) as f:
            graph = json.load(f)["@graph"]
        actions = [e["@id"] for e in graph if e["@type"] == "CreateAction"]
        self.assertEqual(sorted(a.split("-")[1] for a in actions if a.startswith("#fetch-")),
                         ["blacklist", "clinvar"])
        self.assertEqual(len([a for a in actions if a.startswith("#derive-")]), 2)
        self.assertTrue([a for a in actions if a.startswith("#download-")])
        files = [e for e in graph if e["@type"] == "File"]
        self.assertTrue(any(e["@id"].endswith(".vcf.gz.tbi") for e in files))
        for entity in files:
            path = os.path.join(outdir, entity["@id"])
            self.assertTrue(os.path.isfile(path), entity["@id"])
            self.assertEqual(entity["contentSize"], str(os.path.getsize(path)), entity["@id"])
        assets = load_vault(outdir)["toyAsm_ensembl"]["assets"]
        self.assertTrue({"clinvar", "blacklist", "fasta_txome"} <= set(assets))


@unittest.skipIf(ingest_unavailable(), ingest_unavailable())
class ReshapedResourcesAndDatabases(NextflowTestCase):
    """
    The resources that are not one file stored as fetched, and the databases tree.

    dbSNP has to be renamed before it fits, and a VEP cache is a directory that
    follows the annotation release. Under -stub the fetch writes a file named the
    RefSeq way with the table that maps it back, and a cache with a directory per
    sequence; the renaming, both checks, the unpacking and the filing are real.
    """

    HEADER = HEADER.rstrip("\n") + ",url,annotation_url,resources\n"
    REQUESTS = (HEADER
                + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse_ens,,,vep_cache:GRCm39;dbsnp:GRCh38\n"
                + "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,110,mouse_ens,,,vep_cache:GRCm39;dbsnp:GRCh38\n"
                + "Foo_bar,seqOnly,NCBI,false,seqOnly,,,,,\n")
    DATABASES = "# what this run should hold\ndatabase\nswissprot\nkraken2_viral\n"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", cls.REQUESTS)
        cls.databases = cls.write_requests("databases.csv", cls.DATABASES)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir, "--databases", cls.databases,
                                      register=True)

    def submitted(self, process, output=None):
        return sorted(line.split(f"{process} (")[1].rstrip(")")
                      for line in (output or self.run_result.stdout).splitlines()
                      if f"Submitted process > {process} (" in line)

    def test_01_run_succeeds(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])

    def test_02_a_cache_is_fetched_per_annotation_release_and_dbsnp_once(self):
        self.assertEqual(self.submitted("FETCH_RESOURCE"),
                         ["dbsnp:GRCh38 for toyAsm_ensembl", "vep_cache:GRCm39 for toyAsm_ensembl",
                          "vep_cache:GRCm39 for toyAsm_ensembl"])
        self.assertEqual(self.submitted("REGISTER_ASSET"),
                         ["toyAsm_ensembl/dbsnp:stub1", "toyAsm_ensembl/vep_cache:110",
                          "toyAsm_ensembl/vep_cache:999"])

    def test_03_cache_is_unpacked_in_the_vault_and_reachable_as_a_directory(self):
        cache = load_vault(self.outdir)["toyAsm_ensembl"]["assets"]["vep_cache"]
        self.assertEqual(str(cache["default_tag"]), "999")
        self.assertEqual(cache["tags"]["999"]["seek_keys"]["vep_cache"], ".")
        self.assertEqual([p.split("/")[1] for p in cache["tags"]["999"]["asset_parents"]],
                         ["fasta:default"])
        for release in ("110", "999"):
            directory = os.path.join(self.outdir, "alias", "mouse_ens", "vep_cache", release)
            root = os.path.join(directory, "stub_species", f"{release}_stub")
            self.assertTrue(os.path.isfile(os.path.join(root, "info.txt")), root)
            self.assertEqual(sorted(name for name in os.listdir(root) if name != "info.txt"),
                             ["1", "2", "Un_toyAsm"])
            self.assertTrue(os.path.isfile(os.path.join(root, "1", "1-1000000.gz")))
            # The tarball itself is not kept beside what was unpacked from it.
            self.assertEqual(sorted(os.listdir(directory)),
                             ["_refgenie_build", "mouse_ens.provenance.json", "stub_species"])
        seek = subprocess.run(
            [sys.executable, REFGENIE, "seek", "mouse_ens/vep_cache", "-c",
             os.path.join(self.outdir, "refgenie_export.yaml")],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.strip()
        self.assertTrue(os.path.isdir(os.path.join(seek, "stub_species", "999_stub")), seek)

    def test_04_dbsnp_is_stored_under_the_genomes_own_sequence_names(self):
        directory = os.path.join(self.outdir, "alias", "mouse_ens", "dbsnp", "stub1")
        with gzip.open(os.path.join(directory, "mouse_ens.vcf.gz"), "rt") as f:
            records = [line.split("\t")[0] for line in f if not line.startswith("#")]
        self.assertEqual(records, ["1", "1"])
        self.assertTrue(os.path.isfile(os.path.join(directory, "mouse_ens.vcf.gz.tbi")))
        with open(os.path.join(directory, "mouse_ens.provenance.json")) as f:
            record = json.load(f)
        self.assertEqual(record["renamed"]["records_renamed"], 2)
        self.assertEqual(record["renamed"]["examples"], {"NC_STUB.1": "1"})
        self.assertFalse(record["stored_as_published"])

    def test_05_databases_are_filed_beside_the_vault_with_a_catalogue(self):
        self.assertEqual(self.submitted("FETCH_DATABASE"), ["kraken2_viral", "swissprot"])
        self.assertEqual(self.submitted("FILE_DATABASE"), ["kraken2_viral", "swissprot"])
        tree = os.path.join(self.outdir, "databases")
        self.assertEqual(sorted(os.listdir(os.path.join(tree, "swissprot", "stub1"))),
                         ["provenance.json", "swissprot.fasta.gz"])
        # Published as an archive, so unpacked: the directory is what a tool is given.
        self.assertEqual(sorted(os.listdir(os.path.join(tree, "kraken2_viral", "stub1"))),
                         ["hash.k2d", "provenance.json", "taxo.k2d"])
        with open(os.path.join(tree, "databases.tsv")) as f:
            lines = f.read().splitlines()
        self.assertEqual([line.split("\t")[:3] for line in lines[1:]],
                         [["kraken2_viral", "stub1", "databases/kraken2_viral/stub1/"],
                          ["swissprot", "stub1", "databases/swissprot/stub1/"]])

    def test_06_crate_holds_the_databases_and_the_registrations_that_ran_beside_them(self):
        with open(os.path.join(self.outdir, "ro-crate-metadata.json")) as f:
            graph = {entity["@id"]: entity for entity in json.load(f)["@graph"]}
        actions = sorted(key for key, entity in graph.items() if entity["@type"] == "CreateAction")
        self.assertEqual([a for a in actions if a.startswith("#fetch-database-")],
                         ["#fetch-database-kraken2_viral-stub1", "#fetch-database-swissprot-stub1"])
        self.assertEqual(sorted(a.split("-")[1] for a in actions
                                if a.startswith("#fetch-") and "database" not in a),
                         ["dbsnp", "vep_cache", "vep_cache"])
        # A cache is described as its directory, and not file by file.
        fetch = next(graph[a] for a in actions if a.startswith("#fetch-vep_cache-") and a.endswith("-999"))
        directory = fetch["result"][0]["@id"]
        self.assertTrue(directory.endswith("/vep_cache/999/"), directory)
        self.assertFalse(any("stub_species" in key for key in graph))
        for key, entity in graph.items():
            if entity["@type"] == "File":
                self.assertTrue(os.path.isfile(os.path.join(self.outdir, key)), key)

    def test_07_second_run_fetches_registers_and_files_nothing(self):
        result = self.nextflow(self.requests, self.outdir, "--databases", self.databases,
                               register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertEqual(result.stdout.count("  SKIP   vep_cache:GRCm39"), 2)
        self.assertEqual(result.stdout.count("  SKIP   dbsnp:GRCh38"), 1)
        self.assertEqual(result.stdout.count("  SKIP   swissprot"), 1)
        self.assertEqual(result.stdout.count("  SKIP   kraken2_viral"), 1)
        for process in ("PREPARE_RESOURCE", "REGISTER_ASSET", "FILE_DATABASE"):
            self.assertNotIn(f"Submitted process > {process}", result.stdout)

    def test_07b_resumed_run_asks_the_publishers_again(self):
        # run_pipeline.sh always passes -resume. This run's fetch tasks have the
        # same genome, resource and versions held as the run before it, and would
        # be replayed from the cache - "still current", unasked - if nothing told
        # them apart. A publisher's new release would then never be seen.
        result = self.nextflow(self.requests, self.outdir, "--databases", self.databases,
                               "-resume", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertEqual(len(self.submitted("FETCH_RESOURCE", result.stdout)), 3)
        self.assertEqual(len(self.submitted("FETCH_DATABASE", result.stdout)), 2)
        self.assertNotIn("Cached process > FETCH_", result.stdout)

    def test_08_cache_for_a_genome_with_no_annotation_fails_once_and_says_why(self):
        requests = self.write_requests(
            "no_annotation.csv", self.HEADER + "Foo_bar,seqOnly,NCBI,false,seqOnly,,,,,vep_cache:GRCm39\n")
        result = self.nextflow(requests, os.path.join(self.work, "none.out"), register=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("follows an annotation release", result.stdout)
        self.assertEqual(result.stdout.count("Submitted process > FETCH_RESOURCE"), 1)
        self.assertNotIn("Re-submitted process > FETCH_RESOURCE", result.stdout)

    def test_09_database_with_no_definition_stops_the_run_before_any_task(self):
        listed = self.write_requests("unknown.csv", "database\nswissprot\nuniref\n")
        result = self.nextflow(self.requests, os.path.join(self.work, "unknown.out"),
                               "--databases", listed, register=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unknown database 'uniref'", result.stdout)
        self.assertIn("kraken2_standard, kraken2_viral, pfam, silva, swissprot", result.stdout)
        self.assertNotIn("Submitted process", result.stdout)


def gffread_unavailable():
    return ingest_unavailable() or (
        None if shutil.which("gffread") else "gffread is not on PATH")


@unittest.skipIf(ingest_unavailable(), ingest_unavailable())
class DerivedAssets(NextflowTestCase):
    """
    The small derived assets: tx2gene, gtf_indexed, bed12 and refflat from each
    annotation, twobit from each genome.

    The class's main run passes no --derived flag: they are part of every run by
    default, and this is the test that they are. Under -stub the files themselves
    are stand-ins; which genomes get what, the staging and the registration are
    real. tests/test_derive.py runs the real tools.
    """

    FROM_ANNOTATION = ["bed12", "gtf_indexed", "refflat", "tx2gene"]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.outdir = os.path.join(cls.work, "out")
        cls.requests = cls.write_requests("requests_resolved.csv", REQUESTS)
        cls.run_result = cls.nextflow(cls.requests, cls.outdir, register=True,
                                      default_derived=True)

    def submitted(self, process, output=None):
        return sorted(line.split(f"{process} (")[1].rstrip(")")
                      for line in (output or self.run_result.stdout).splitlines()
                      if f"Submitted process > {process} (" in line)

    def refgenie(self, *arguments):
        env = dict(os.environ, PATH=BIN_DIR + os.pathsep + os.environ.get("PATH", ""))
        result = subprocess.run(
            [sys.executable, REFGENIE, *arguments, "-c",
             os.path.join(self.outdir, "refgenie_export.yaml")],
            cwd=self.work, env=env, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        # `seek` prints a path; `remove` prints nothing on stdout.
        lines = result.stdout.strip().splitlines()
        return lines[-1] if lines else ""

    def test_01_run_succeeds(self):
        self.assertEqual(self.run_result.returncode, 0, self.run_result.stdout[-3000:])

    def test_02_one_task_per_annotation_and_one_per_genome(self):
        # Both Ensembl releases and the GENCODE genome have an annotation to derive
        # from. The sequence-only genome does not, and still gets its .2bit.
        self.assertEqual(self.submitted("DERIVE_ANNOTATION"),
                         ["toyAsm_ensembl:110", "toyAsm_ensembl:999", "toyAsm_gencode:999"])
        self.assertEqual(self.submitted("DERIVE_TWOBIT"),
                         ["seqOnly_ncbi", "toyAsm_ensembl", "toyAsm_gencode"])
        self.assertEqual(len(self.submitted("REGISTER_ASSET")), 3 * 4 + 3)

    def test_03_annotation_assets_are_tagged_by_release_with_the_annotation_as_parent(self):
        vault = load_vault(self.outdir)
        with open(os.path.join(self.outdir, "refgenie_export.yaml")) as f:
            genomes = yaml.safe_load(f)["genomes"]
        digest = next(d for d, g in genomes.items() if "toyAsm_ensembl" in g["aliases"])
        for asset in self.FROM_ANNOTATION:
            held = vault["toyAsm_ensembl"]["assets"][asset]
            self.assertEqual(set(held["tags"]), {"110", "999"}, asset)
            self.assertEqual(str(held["default_tag"]), "999", asset)
            for release in ("110", "999"):
                self.assertEqual(held["tags"][release]["asset_parents"],
                                 [f"{digest}/ensembl_gtf:{release}"], asset)
            self.assertNotIn(asset, vault["seqOnly_ncbi"]["assets"])

    def test_04_twobit_is_per_genome_with_the_fasta_as_parent(self):
        vault = load_vault(self.outdir)
        for genome in ("toyAsm_ensembl", "toyAsm_gencode", "seqOnly_ncbi"):
            held = vault[genome]["assets"]["twobit"]
            self.assertEqual(set(held["tags"]), {"default"}, genome)
            self.assertEqual([p.split("/")[1] for p in held["tags"]["default"]["asset_parents"]],
                             ["fasta:default"], genome)

    def test_05_each_is_reachable_by_alias_under_a_name_that_says_what_it_is(self):
        expected = {
            "tx2gene": "mouse_gc.tx2gene.tsv",
            "gtf_indexed": "mouse_gc.sorted.gtf.gz",
            "gtf_indexed.tbi": "mouse_gc.sorted.gtf.gz.tbi",
            "bed12": "mouse_gc.bed",
            "refflat": "mouse_gc.refFlat.txt",
            "twobit": "mouse_gc.2bit",
            "tx2gene.provenance": "mouse_gc.provenance.json",
        }
        for seek, name in expected.items():
            path = self.refgenie("seek", f"mouse_gc/{seek}")
            self.assertEqual(os.path.basename(path), name)
            self.assertTrue(os.path.isfile(path), path)
        with open(self.refgenie("seek", "mouse_gc/tx2gene")) as f:
            self.assertEqual([line.split("\t")[0] for line in f], ["t1"])
        with open(self.refgenie("seek", "mouse_gc/bed12.provenance")) as f:
            record = json.load(f)
        self.assertEqual((record["asset"], record["annotation"]), ("bed12", "gencode_gtf:999"))

    def test_06_index_and_crate_list_them(self):
        with open(os.path.join(self.outdir, "provenance", "index.tsv")) as f:
            index = f.read()
        for held in ("tx2gene:110,999", "gtf_indexed:110,999", "twobit:default"):
            self.assertIn(held, index)
        with open(os.path.join(self.outdir, "ro-crate-metadata.json")) as f:
            graph = json.load(f)["@graph"]
        derived = sorted(e["@id"].split("-")[1] for e in graph
                         if e["@type"] == "CreateAction" and e["@id"].startswith("#derive-"))
        self.assertEqual(derived, sorted(self.FROM_ANNOTATION * 3 + ["twobit"] * 3))
        tx2gene = next(e for e in graph if e["@id"].startswith("#derive-tx2gene-")
                       and e["@id"].endswith("-110"))
        self.assertTrue(tx2gene["object"]["@id"].endswith(".gtf.gz"), tx2gene["object"])
        self.assertIn("/ensembl_gtf/110/", tx2gene["object"]["@id"])

    def test_07_second_run_derives_and_registers_nothing(self):
        result = self.nextflow(self.requests, self.outdir, register=True, default_derived=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        for process in ("DERIVE_ANNOTATION", "DERIVE_TWOBIT", "REGISTER_ASSET"):
            self.assertNotIn(process, result.stdout)

    def test_08_only_what_the_vault_lacks_is_registered(self):
        # One asset of one genome removed by hand. The task that makes all four
        # runs again for that annotation, and only the missing one is registered.
        self.refgenie("remove", "toyAsm_gencode/bed12:999", "-f")
        self.assertNotIn("bed12", load_vault(self.outdir)["toyAsm_gencode"]["assets"])
        result = self.nextflow(self.requests, self.outdir, register=True, default_derived=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertEqual(self.submitted("DERIVE_ANNOTATION", result.stdout), ["toyAsm_gencode:999"])
        self.assertEqual(self.submitted("DERIVE_TWOBIT", result.stdout), [])
        self.assertEqual(self.submitted("REGISTER_ASSET", result.stdout),
                         ["toyAsm_gencode/bed12:999"])
        self.assertIn("bed12", load_vault(self.outdir)["toyAsm_gencode"]["assets"])

    def test_09_derived_false_makes_none_and_a_later_run_adds_them(self):
        outdir = os.path.join(self.work, "later.out")
        without = self.nextflow(self.requests, outdir, register=True)
        self.assertEqual(without.returncode, 0, without.stdout[-3000:])
        self.assertNotIn("DERIVE_", without.stdout)
        self.assertNotIn("twobit", load_vault(outdir)["seqOnly_ncbi"]["assets"])
        result = self.nextflow(self.requests, outdir, register=True, default_derived=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("FETCH_GENOME", result.stdout)
        self.assertEqual(len(self.submitted("DERIVE_ANNOTATION", result.stdout)), 3)
        self.assertIn("twobit", load_vault(outdir)["seqOnly_ncbi"]["assets"])


@unittest.skipIf(gffread_unavailable(), gffread_unavailable())
class RealDerivation(NextflowTestCase):
    """
    DERIVE_TXOME with the real gffread, on a genome built for the purpose.

    No -stub here. The genome is put in the vault first by the ingest, so the
    preflight skips the download and the workflow goes straight from REGISTER to
    the real derivation. Needs gffread on PATH, which the pipeline's Python
    environment does not have: run the suite inside an environment that does
    (`conda create -n gffread -c bioconda gffread=0.12.7`) to include these.

    Contig 1 is half lower case, as a soft-masked genome is. Transcript t2 is on
    the minus strand with an intron, so a wrong answer cannot pass by accident.
    """

    CONTIG_1 = "ACGTacgt" * 5                    # 40 bases
    CONTIG_2 = "AAAACCCCGGGGTTTTACGTAAAACCCCGGGGTTTTACGT"
    HEADER_ROW = "Mus_musculus,toyAsm,Ensembl,TRUE,toyAsm,,mouse\n"

    @staticmethod
    def attributes(transcript):
        return (f'gene_id "g_{transcript}"; gene_version "1"; transcript_id "{transcript}"; '
                f'transcript_version "1"; exon_number "1"; gene_name "G"; gene_source "s"; '
                f'gene_biotype "lncRNA";')

    def build_vault(self, name, gtf_contig_prefix=""):
        """Ingest a hand-written download; returns (outdir, resolved CSV)."""
        import gzip
        outdir = os.path.join(self.work, name)
        download = os.path.join(outdir, "Ensembl", "Mus_musculus", "toyAsm")
        os.makedirs(download)
        with gzip.open(os.path.join(download, "toyAsm.fa.gz"), "wt") as f:
            f.write(f">1\n{self.CONTIG_1}\n>2\n{self.CONTIG_2}\n")
        p = gtf_contig_prefix
        with gzip.open(os.path.join(download, "toyAsm.annotation.gtf.gz"), "wt") as f:
            f.write(f"{p}1\tsrc\texon\t1\t30\t.\t+\t.\t{self.attributes('t1')}\n")
            f.write(f"{p}2\tsrc\texon\t1\t10\t.\t-\t.\t{self.attributes('t2')}\n")
            f.write(f"{p}2\tsrc\texon\t21\t30\t.\t-\t.\t{self.attributes('t2')}\n")
        with open(os.path.join(download, "provenance.json"), "w") as f:
            json.dump({"assembly": "toyAsm", "annotation_release": "116",
                       "annotation_release_kind": "ensembl_release"}, f)
        resolved = self.write_requests(name + ".csv", HEADER + self.HEADER_ROW)
        env = dict(os.environ, PATH=BIN_DIR + os.pathsep + os.environ.get("PATH", ""))
        result = subprocess.run([sys.executable, INGEST, resolved, outdir], cwd=self.work,
                                env=env, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        return outdir, resolved

    def test_transcripts_are_the_spliced_strand_aware_upper_case_sequence(self):
        outdir, resolved = self.build_vault("real.out")
        result = self.nextflow(resolved, outdir, "--txome", "true", register=True, stub=False)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertNotIn("FETCH_GENOME", result.stdout)

        complement = str.maketrans("ACGT", "TGCA")
        spliced = self.CONTIG_2[0:10] + self.CONTIG_2[20:30]
        expected = {
            "t1": self.CONTIG_1[0:30].upper(),
            "t2": spliced.translate(complement)[::-1],
        }
        directory = os.path.join(outdir, "alias", "mouse", "fasta_txome", "116")
        self.assertEqual(read_fasta(os.path.join(directory, "mouse.fa")), expected)

        with open(os.path.join(directory, "mouse.provenance.json")) as f:
            record = json.load(f)
        self.assertEqual(record["tool"], "gffread")
        self.assertRegex(record["tool_version"], r"\d+\.\d+")
        self.assertEqual((record["transcripts"], record["sequence_case"]), (2, "upper"))
        self.assertEqual(record["annotation"], "ensembl_gtf:116")
        # When it was made, in the form the other derived assets use; it is what
        # gives the derivation an end time in the vault's RO-Crate.
        self.assertRegex(record["derived_utc"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\+00:00$")
        with open(os.path.join(outdir, "ro-crate-metadata.json")) as f:
            derivation, = [e for e in json.load(f)["@graph"]
                           if e["@id"].startswith("#derive-fasta_txome-")]
        self.assertEqual(derivation["endTime"], record["derived_utc"])

    def test_annotation_whose_contigs_are_not_in_the_genome_is_an_error(self):
        # A chr-named GTF beside a genome whose contigs are 1 and 2. gffread
        # exits 0 having written nothing; filing that as a transcriptome would
        # put an empty asset in the vault under a confident name.
        outdir, resolved = self.build_vault("mismatch.out", gtf_contig_prefix="chr")
        result = self.nextflow(resolved, outdir, "--txome", "true", register=True, stub=False)
        self.assertNotEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertIn("wrote no transcripts", result.stdout)
        self.assertNotIn("fasta_txome", load_vault(outdir)["toyAsm_ensembl"]["assets"])


class RequestValidation(NextflowTestCase):
    def rejected(self, name, content, *extra, **options):
        result = self.nextflow(self.write_requests(name, content),
                               os.path.join(self.work, name + ".out"), *extra, **options)
        self.assertNotEqual(result.returncode, 0, result.stdout[-2000:])
        self.assertNotIn("Submitted process", result.stdout, "a task started before the check")
        return result.stdout

    def test_missing_python_environment_stops_the_run_before_any_download(self):
        # Found after the downloads, this would cost hours.
        output = self.rejected("novenv.csv", HEADER + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n",
                               register=True, venv_bin=os.path.join(self.work, "no-such-venv"))
        self.assertIn("No Python environment", output)

    def test_number_parameter_that_is_not_a_number_is_refused(self):
        output = self.rejected("badnumber.csv", HEADER + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n",
                               "--max_retries", "lots")
        self.assertIn("--max_retries must be a whole number", output)

    def test_duplicate_download_is_refused_before_any_task_starts(self):
        output = self.rejected("duplicate.csv", HEADER
                               + "Mus_musculus,toyAsm,Ensembl,TRUE,a,,\n"
                               + "Mus_musculus,toyAsm,ensembl,false,b,,\n")
        self.assertIn("Duplicate request for toyAsm", output)

    def test_same_assembly_at_two_releases_is_not_a_duplicate(self):
        result = self.nextflow(
            self.write_requests("releases.csv", HEADER
                                + "Mus_musculus,toyAsm,Ensembl,TRUE,a,115,\n"
                                + "Mus_musculus,toyAsm,Ensembl,TRUE,a,116,\n"),
            os.path.join(self.work, "releases.out"))
        self.assertEqual(result.returncode, 0, result.stdout[-2000:])

    URL_HEADER = HEADER.rstrip("\n") + ",url,annotation_url\n"

    def test_address_with_shell_characters_is_refused(self):
        output = self.rejected("badurl.csv", self.URL_HEADER
                               + "Homo_sapiens,noalt,URL,false,noalt,,,https://example.org/a.fa.gz;id,\n")
        self.assertIn("Unusable url", output)

    def test_genome_by_address_without_an_address_is_refused(self):
        output = self.rejected("nourl.csv", self.URL_HEADER
                               + "Homo_sapiens,noalt,URL,false,noalt,,,,\n")
        self.assertIn("needs an address in the url column", output)

    def test_address_on_a_named_provider_is_refused(self):
        output = self.rejected("strayurl.csv", self.URL_HEADER
                               + "Mus_musculus,toyAsm,Ensembl,false,a,,,https://example.org/a.fa.gz,\n")
        self.assertIn("only used with provider 'url'", output)

    def test_non_numeric_release_is_refused(self):
        output = self.rejected("release.csv", HEADER
                               + "Mus_musculus,toyAsm,Ensembl,TRUE,a,../110,\n")
        self.assertIn("non-numeric release", output)

    def test_unsafe_species_is_refused(self):
        output = self.rejected("species.csv", HEADER
                               + "Mus musculus; rm,toyAsm,Ensembl,TRUE,a,,\n")
        self.assertIn("Unsafe characters detected in species", output)

    def test_older_request_file_without_a_release_column_runs(self):
        result = self.nextflow(
            self.write_requests("old.csv", "species,assembly,provider,annotation\n"
                                           "Foo_bar,seqOnly,NCBI,false\n"),
            os.path.join(self.work, "old.out"))
        self.assertEqual(result.returncode, 0, result.stdout[-2000:])
        self.assertTrue(os.path.isfile(
            os.path.join(self.work, "old.out", "NCBI", "Foo_bar", "seqOnly", "seqOnly.fa.gz")))


class Retries(NextflowTestCase):
    """The real script block, with `python3` replaced by something that only fails."""

    def failing_download(self, exit_code, *extra):
        fake_bin = os.path.join(self.work, f"fake-python-exit-{exit_code}")
        os.makedirs(fake_bin, exist_ok=True)
        fake = os.path.join(fake_bin, "python3")
        with open(fake, "w") as f:
            f.write(f"#!/bin/sh\nexit {exit_code}\n")
        os.chmod(fake, os.stat(fake).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        name = f"fail-{exit_code}-{len(extra)}"
        requests = self.write_requests(
            name + ".csv", HEADER + "Foo_bar,seqOnly,NCBI,false,seqOnly,,\n")
        env = dict(os.environ, PATH=fake_bin + os.pathsep + os.environ.get("PATH", ""))
        return self.nextflow(requests, os.path.join(self.work, name + ".out"),
                             "--retry_wait", "0", *extra, stub=False, env=env)

    def test_failed_download_is_retried_then_fails_the_run(self):
        result = self.failing_download(1)
        self.assertNotEqual(result.returncode, 0, result.stdout[-2000:])
        # max_retries defaults to 2: three attempts, so two retries and no third.
        self.assertIn("Execution is retried (1)", result.stdout)
        self.assertIn("Execution is retried (2)", result.stdout)
        self.assertNotIn("Execution is retried (3)", result.stdout)

    def test_missing_tool_exit_code_is_not_retried(self):
        result = self.failing_download(3)
        self.assertNotEqual(result.returncode, 0, result.stdout[-2000:])
        self.assertNotIn("Execution is retried", result.stdout)

    def test_retries_can_be_turned_off(self):
        result = self.failing_download(1, "--max_retries", "0")
        self.assertNotEqual(result.returncode, 0, result.stdout[-2000:])
        self.assertNotIn("Execution is retried", result.stdout)


if __name__ == "__main__":
    unittest.main()
