# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of bin/preflight.py: which requests the vault already satisfies.

The vault here is a hand-written config plus the few directories and files the
preflight looks for, so no genome is built and no provider is contacted. Where a
test needs "the release the provider publishes now", it supplies the answer.
tests/test_nextflow_stub.py runs the preflight for real, inside the workflow.
"""
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import preflight  # noqa: E402

ENSEMBL, GENCODE, NCBI = "e" * 48, "c" * 48, "n" * 48


def fasta_asset(filename):
    return {"default_tag": "default",
            "tags": {"default": {"seek_keys": {"fasta": "x.fa", "dict": filename}}}}


class VaultTestCase(unittest.TestCase):
    """
    Three genomes, each wrong or right in a different way:

      GRCm39_ensembl  complete, with annotation releases 115 and 116
      GRCm39_gencode  complete genome; the config lists annotation M33, but that
                      release has no directory in the vault
      ASM1_ncbi       in the config, but its sequence dictionary was never written
    """

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.outdir = directory.name
        self.config_file = os.path.join(self.outdir, "refgenie_export.yaml")
        config = {
            "genomes": {
                ENSEMBL: {"aliases": ["GRCm39_ensembl", "GRCm39"],
                          "assets": {"fasta": fasta_asset("e.dict"),
                                     "ensembl_gtf": {"default_tag": "116",
                                                     "tags": {"115": {}, "116": {}}}}},
                GENCODE: {"aliases": ["GRCm39_gencode"],
                          "assets": {"fasta": fasta_asset("c.dict"),
                                     "gencode_gtf": {"default_tag": "M33", "tags": {"M33": {}}}}},
                NCBI: {"aliases": ["ASM1_ncbi"], "assets": {"fasta": fasta_asset("n.dict")}},
            }
        }
        with open(self.config_file, "w") as f:
            yaml.dump(config, f)
        self.config = config

        for digest, filename in ((ENSEMBL, "e.dict"), (GENCODE, "c.dict")):
            asset_dir = os.path.join(self.outdir, "data", digest, "fasta", "default")
            os.makedirs(asset_dir)
            with open(os.path.join(asset_dir, filename), "w") as f:
                f.write("@HD\tVN:1.0\n")
        for tag in ("115", "116"):
            os.makedirs(os.path.join(self.outdir, "data", ENSEMBL, "ensembl_gtf", tag))

    def assess(self, assembly, provider, annotation="TRUE", release="", current="unasked",
               check_releases=True, config="vault"):
        """Returns (fetch, reason, providers the release question was put to)."""
        asked = []

        def current_release(asked_provider, asked_assembly):
            asked.append(asked_provider)
            return current

        row = {"assembly": assembly, "provider": provider,
               "annotation": annotation, "release": release}
        fetch, reason = preflight.assess(
            row, self.config if config == "vault" else config, self.config_file,
            current_release=current_release, check_releases=check_releases)
        return fetch, reason, asked


class GenomePresence(VaultTestCase):
    def test_everything_is_fetched_when_there_is_no_vault(self):
        fetch, reason, _ = self.assess("GRCm39", "Ensembl", config=None)
        self.assertTrue(fetch)
        self.assertIn("no vault", reason)

    def test_genome_not_in_the_vault_is_fetched(self):
        fetch, reason, _ = self.assess("sacCer3", "UCSC")
        self.assertTrue(fetch)
        self.assertIn("sacCer3_ucsc", reason)

    def test_same_assembly_from_another_provider_does_not_count(self):
        # 'GRCm39' is an alias in the vault - of the Ensembl genome.
        fetch, _, _ = self.assess("GRCm39", "NCBI", annotation="false")
        self.assertTrue(fetch)

    def test_sequence_only_request_is_satisfied_by_the_genome(self):
        fetch, _, asked = self.assess("GRCm39", "Ensembl", annotation="false")
        self.assertFalse(fetch)
        self.assertEqual(asked, [])

    def test_genome_without_its_dictionary_on_disk_is_fetched(self):
        # The config lists the seek key; the file was never written.
        fetch, reason, _ = self.assess("ASM1", "NCBI", annotation="false")
        self.assertTrue(fetch)
        self.assertIn("dictionary", reason)


class PinnedRelease(VaultTestCase):
    def test_pinned_release_in_the_vault_is_skipped_without_asking_anyone(self):
        fetch, _, asked = self.assess("GRCm39", "Ensembl", release="115")
        self.assertFalse(fetch)
        self.assertEqual(asked, [])

    def test_pinned_release_not_in_the_vault_is_fetched(self):
        fetch, reason, _ = self.assess("GRCm39", "Ensembl", release="110")
        self.assertTrue(fetch)
        self.assertIn("110", reason)


class UnpinnedRelease(VaultTestCase):
    def test_skipped_when_the_vault_has_the_current_release(self):
        fetch, _, asked = self.assess("GRCm39", "Ensembl", current="116")
        self.assertFalse(fetch)
        self.assertEqual(asked, ["Ensembl"])

    def test_fetched_when_the_provider_has_moved_on(self):
        fetch, reason, _ = self.assess("GRCm39", "Ensembl", current="117")
        self.assertTrue(fetch)
        self.assertIn("117", reason)

    def test_skipped_when_the_current_release_cannot_be_determined(self):
        # Offline, or the provider is down. The genome is there; do not refetch it.
        fetch, reason, _ = self.assess("GRCm39", "Ensembl", current=None)
        self.assertFalse(fetch)
        self.assertIn("could not be checked", reason)

    def test_release_check_can_be_turned_off(self):
        fetch, _, asked = self.assess("GRCm39", "Ensembl", current="117", check_releases=False)
        self.assertFalse(fetch)
        self.assertEqual(asked, [])

    def test_release_listed_in_the_config_but_absent_from_disk_is_fetched(self):
        fetch, reason, asked = self.assess("GRCm39", "GENCODE", current="M33")
        self.assertTrue(fetch)
        self.assertIn("gencode_gtf", reason)
        self.assertEqual(asked, [], "nothing to compare a current release with")


class CurrentRelease(unittest.TestCase):
    """The tag comes out of the same rules that tag an asset at ingest."""

    def release_for(self, provider, source):
        with mock.patch.object(preflight, "annotation_source", return_value=source), \
                redirect_stderr(io.StringIO()):
            return preflight.current_annotation_release(provider, "GRCm39")

    def test_ensembl_release_number_from_the_download_url(self):
        self.assertEqual(self.release_for(
            "Ensembl",
            "http://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"),
            "116")

    def test_gencode_release_from_the_download_path(self):
        self.assertEqual(self.release_for(
            "GENCODE",
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/"
            "gencode.vM33.annotation.gtf.gz"),
            "M33")

    def test_any_failure_is_reported_as_unknown_not_raised(self):
        with mock.patch.object(preflight, "annotation_source",
                               side_effect=ConnectionError("Ensembl appears to be offline")), \
                redirect_stderr(io.StringIO()) as err:
            self.assertIsNone(preflight.current_annotation_release("Ensembl", "GRCm39"))
        self.assertIn("Could not determine", err.getvalue())


class CommandLine(VaultTestCase):
    def run_preflight(self, rows, *extra):
        resolved = os.path.join(self.outdir, "requests_resolved.csv")
        with open(resolved, "w", newline="") as f:
            f.write("species,assembly,provider,annotation,query,release,aliases\n"
                    + "\n".join(rows) + "\n")
        keys_file = os.path.join(self.outdir, "fetch_keys.txt")
        out = io.StringIO()
        argv = ["preflight.py", resolved, self.outdir, "--fetch-keys", keys_file, *extra]
        with mock.patch.object(sys, "argv", argv), redirect_stdout(out), \
                redirect_stderr(io.StringIO()):
            preflight.main()
        with open(keys_file) as f:
            return f.read().splitlines(), out.getvalue()

    def test_writes_the_key_of_each_request_to_fetch_and_reports_every_row(self):
        keys, report = self.run_preflight([
            "Mus_musculus,GRCm39,Ensembl,TRUE,GRCm39,,mouse",
            "Mus_musculus,GRCm39,Ensembl,TRUE,GRCm39,110,mouse",
            "Saccharomyces_cerevisiae,sacCer3,UCSC,false,sacCer3,,yeast",
        ], "--no-release-check")
        # The format main.nf builds and REGISTER's --source uses: provider in lower
        # case, and '@' followed by the pinned release or by nothing.
        self.assertEqual(keys, ["ensembl/GRCm39@110", "ucsc/sacCer3@"])
        self.assertEqual(report.count("  SKIP "), 1)
        self.assertEqual(report.count("  FETCH"), 2)
        self.assertIn("2 to download, 1 already in the vault", report)

    def test_a_failure_while_checking_a_row_fetches_it(self):
        with mock.patch.object(preflight, "assess", side_effect=RuntimeError("boom")):
            keys, report = self.run_preflight(
                ["Mus_musculus,GRCm39,Ensembl,TRUE,GRCm39,,mouse"], "--no-release-check")
        self.assertEqual(keys, ["ensembl/GRCm39@"])
        self.assertIn("the check itself failed", report)

    def test_no_requests_writes_an_empty_file(self):
        keys, report = self.run_preflight([], "--no-release-check")
        self.assertEqual(keys, [])
        self.assertIn("0 to download", report)


if __name__ == "__main__":
    unittest.main()
