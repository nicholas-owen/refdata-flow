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
                 register=False, venv_bin=BIN_DIR):
        # Most tests here are about the download side and pass --register false, so
        # they run wherever nextflow does. The Registration tests turn it on and
        # point REGISTER at the interpreter running this suite, which is the
        # pipeline's own environment when the suite is run as the README says.
        command = ["nextflow", "run", MAIN_NF, "-profile", "standard", "-ansi-log", "false",
                   "--input", requests, "--outdir", outdir]
        command += ["--venv_bin", venv_bin] if register else ["--register", "false"]
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
        self.assertEqual(set(by_name["toyAsm_gencode"]["assets"]), {"fasta", "gencode_gtf"})
        self.assertEqual(set(by_name["seqOnly_ncbi"]["assets"]), {"fasta"})
        # Both Ensembl rows reached the vault, each under its own tag.
        self.assertEqual(set(by_name["toyAsm_ensembl"]["assets"]["ensembl_gtf"]["tags"]),
                         {"stub", "110"})


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
        self.assertEqual(set(vault["toyAsm_gencode"]["assets"]), {"fasta", "gencode_gtf"})
        self.assertEqual(set(vault["seqOnly_ncbi"]["assets"]), {"fasta"})
        # The pinned and the unpinned Ensembl download share an assembly name and a
        # directory name. Both tags present means each was matched to its own row.
        self.assertEqual(set(vault["toyAsm_ensembl"]["assets"]["ensembl_gtf"]["tags"]),
                         {"stub", "110"})
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

    def test_05_resumed_run_registers_again_and_changes_nothing(self):
        before = {alias: sorted(genome["assets"]) for alias, genome in load_vault(self.outdir).items()}
        result = self.nextflow(self.requests, self.outdir, "-resume", register=True)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertIn("Cached process > FETCH_GENOME", result.stdout)
        # REGISTER is never cached: the vault may have changed outside the workflow.
        self.assertIn("Submitted process > REGISTER", result.stdout)
        after = {alias: sorted(genome["assets"]) for alias, genome in load_vault(self.outdir).items()}
        self.assertEqual(before, after)

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
