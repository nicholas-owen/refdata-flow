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
                 register=False, venv_bin=BIN_DIR, default_txome=False):
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


def gffread_unavailable():
    return ingest_unavailable() or (
        None if shutil.which("gffread") else "gffread is not on PATH")


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
