# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Unit tests of the naming, alias and default-tag rules in bin/update_refgenie.py.

These need only PyYAML. The refgenie config is a hand-written file, so no vault
is built; tests/test_ingest_e2e.py covers the same rules against real refgenie.
"""
import logging
import os
import sys
import tempfile
import unittest

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import update_refgenie as ingest  # noqa: E402

LOGGER = logging.getLogger("test_update_refgenie")
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False

ENSEMBL, GENCODE = "e" * 48, "c" * 48


class VaultGenomeName(unittest.TestCase):
    def test_names_the_provider_and_replaces_dots(self):
        self.assertEqual(ingest.vault_genome_name("GRCh38.p14", "Ensembl"), "GRCh38_p14_ensembl")

    def test_provider_case_does_not_make_a_second_genome(self):
        self.assertEqual(ingest.vault_genome_name("GRCm39", "Ensembl"),
                         ingest.vault_genome_name("GRCm39", "ensembl"))

    def test_same_assembly_from_two_providers_gets_two_names(self):
        self.assertNotEqual(ingest.vault_genome_name("GRCm39", "Ensembl"),
                            ingest.vault_genome_name("GRCm39", "GENCODE"))

    def test_missing_provider_is_refused(self):
        with self.assertRaises(ValueError):
            ingest.vault_genome_name("GRCm39", "  ")


class DefaultTagOrdering(unittest.TestCase):
    def test_sort_key_reads_each_tag_vocabulary(self):
        self.assertEqual(ingest.release_sort_key("116"), (116,))
        self.assertEqual(ingest.release_sort_key("M33"), (33,))
        self.assertEqual(ingest.release_sort_key("20240312"), (20240312,))
        self.assertIsNone(ingest.release_sort_key("default"))
        self.assertIsNone(ingest.release_sort_key(None))

    def test_newer_release_takes_the_default(self):
        self.assertTrue(ingest.should_become_default("110", "116"))
        self.assertTrue(ingest.should_become_default("M32", "M33"))

    def test_older_release_does_not(self):
        self.assertFalse(ingest.should_become_default("116", "110"))
        self.assertFalse(ingest.should_become_default("20240312", "20230101"))

    def test_comparison_is_numeric_not_alphabetical(self):
        # As strings, '99' sorts above '116'.
        self.assertTrue(ingest.should_become_default("99", "116"))
        self.assertFalse(ingest.should_become_default("116", "99"))

    def test_same_release_again_keeps_the_default(self):
        self.assertTrue(ingest.should_become_default("116", "116"))

    def test_identified_release_replaces_an_unidentified_one(self):
        self.assertTrue(ingest.should_become_default("default", "116"))
        self.assertTrue(ingest.should_become_default(None, "116"))


class ConfigTestCase(unittest.TestCase):
    """A config holding one assembly from two providers, as two digests."""

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.config_file = os.path.join(directory.name, "refgenie_export.yaml")
        self.write({
            "genomes": {
                ENSEMBL: {
                    "aliases": ["GRCm39_ensembl", "GRCm39", "mouse"],
                    "assets": {"ensembl_gtf": {"default_tag": "116",
                                               "tags": {"110": {}, "116": {}}}},
                },
                GENCODE: {"aliases": ["GRCm39_gencode"], "assets": {}},
            }
        })

    def write(self, config):
        with open(self.config_file, "w") as f:
            yaml.dump(config, f)

    def aliases(self, digest):
        with open(self.config_file) as f:
            return set(yaml.safe_load(f)["genomes"][digest]["aliases"])

    def merge(self, vault_name, auto_names=(), requested=()):
        return ingest.merge_aliases(vault_name, "GRCm39", set(auto_names), list(requested),
                                    self.config_file, LOGGER)


class MergeAliases(ConfigTestCase):
    def test_free_aliases_are_added_and_existing_ones_kept(self):
        problems = self.merge("GRCm39_gencode", requested=["mm_gc"])
        self.assertEqual(problems, [])
        self.assertEqual(self.aliases(GENCODE), {"GRCm39_gencode", "mm_gc"})

    def test_bare_assembly_name_stays_with_the_first_genome_without_failing(self):
        # Requesting one assembly from two providers is legitimate. The second
        # must not take 'GRCm39' from the first, and must not fail the run either.
        problems = self.merge("GRCm39_gencode", auto_names=["GRCm39"])
        self.assertEqual(problems, [])
        self.assertNotIn("GRCm39", self.aliases(GENCODE))
        self.assertIn("GRCm39", self.aliases(ENSEMBL))

    def test_requested_alias_owned_by_another_genome_is_a_problem(self):
        problems = self.merge("GRCm39_gencode", requested=["mouse", "mm_gc"])
        self.assertEqual(len(problems), 1)
        self.assertIn("mouse", problems[0])
        self.assertNotIn("mouse", self.aliases(GENCODE))
        self.assertIn("mouse", self.aliases(ENSEMBL))
        # The alias that was free is still applied.
        self.assertIn("mm_gc", self.aliases(GENCODE))

    def test_alias_asked_for_by_name_is_not_downgraded_to_a_derived_one(self):
        # 'GRCm39' is both the assembly name and, here, an explicit request.
        problems = self.merge("GRCm39_gencode", auto_names=["GRCm39"], requested=["GRCm39"])
        self.assertEqual(len(problems), 1)

    def test_genome_is_found_by_vault_name_only(self):
        # 'GRCm39' is in the config, but on the other provider's digest. Matching
        # on it would merge this row's aliases onto the wrong genome.
        problems = self.merge("GRCm39_ncbi", auto_names=["GRCm39"], requested=["mm_ncbi"])
        self.assertEqual(len(problems), 1)
        self.assertIn("GRCm39_ncbi", problems[0])
        self.assertNotIn("mm_ncbi", self.aliases(ENSEMBL))

    def test_merging_twice_changes_nothing(self):
        self.merge("GRCm39_gencode", requested=["mm_gc"])
        before = self.aliases(GENCODE)
        self.assertEqual(self.merge("GRCm39_gencode", requested=["mm_gc"]), [])
        self.assertEqual(self.aliases(GENCODE), before)


class CurrentDefaultTag(ConfigTestCase):
    def test_reads_the_default_through_any_alias_of_the_genome(self):
        self.assertEqual(
            ingest.current_default_tag(self.config_file, "GRCm39_ensembl", "ensembl_gtf"), "116")

    def test_none_when_the_asset_has_not_been_built(self):
        self.assertIsNone(
            ingest.current_default_tag(self.config_file, "GRCm39_gencode", "gencode_gtf"))

    def test_none_when_the_genome_or_config_is_missing(self):
        self.assertIsNone(
            ingest.current_default_tag(self.config_file, "nothing_here", "ensembl_gtf"))
        self.assertIsNone(
            ingest.current_default_tag(self.config_file + ".absent", "GRCm39_ensembl", "x"))


class ParseAliases(unittest.TestCase):
    def test_splits_on_semicolons_and_drops_blanks(self):
        self.assertEqual(ingest.parse_aliases(" yeast ; sacCer;;"), ["yeast", "sacCer"])
        self.assertEqual(ingest.parse_aliases(None), [])


if __name__ == "__main__":
    unittest.main()
