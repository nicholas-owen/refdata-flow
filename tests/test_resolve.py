# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Offline tests of bin/resolve.py.

The provider search is replaced with a fixed table, so nothing here touches the
network; everything else is the real resolver, driven through main() exactly as
run_pipeline.sh drives it. Run with the pipeline's own interpreter, which is the
one that has genomepy:

    .venv/bin/python -m unittest discover -s tests -v
"""
import csv
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import resolve  # noqa: E402

# Rows in the layout genomepy.search() yields:
#   assembly, provider, accession, taxid, annotations, species, other
# `annotations` is a bool for most providers and a list of bools for UCSC; both
# forms are here on purpose.
SEARCH_RESULTS = {
    "GRCm39": [
        ["GRCm39", "Ensembl", "GCA_000001635.9", 10090, True, "Mus musculus", ""],
        ["GRCm39", "GENCODE", "GCA_000001635.9", 10090, True, "Mus musculus", ""],
    ],
    "GRCh38": [["GRCh38.p14", "Ensembl", "GCA_000001405.29", 9606, True, "Homo sapiens", ""]],
    "human": [["GRCh38.p14", "Ensembl", "GCA_000001405.29", 9606, True, "Homo sapiens", ""]],
    "unannotated": [["ASM1v1", "NCBI", "GCA_1.1", 1, False, "Foo bar", ""]],
    "unannotated_ucsc": [
        ["fooBar1", "UCSC", "GCA_2.1", 2, [False, False, False, False], "Foo baz", ""]
    ],
    "sacCer3": [
        ["sacCer3", "UCSC", "GCA_000146045.2", 559292, [True, False, True, False],
         "Saccharomyces cerevisiae", ""]
    ],
}

RESOLVED_HEADER = "species,assembly,provider,annotation,query,release,aliases\n"


def rows_for(query, provider=None):
    """What genomepy.search() yields for a query: one provider's rows, or everyone's."""
    rows = SEARCH_RESULTS.get(query, [])
    if provider:
        rows = [row for row in rows if row[1].lower() == provider.lower()]
    return rows


class ResolverTestCase(unittest.TestCase):
    def setUp(self):
        self.searched = []
        self.searched_providers = []

        def fake_search(query, provider=None, exact=False, size=False):
            self.searched.append(query)
            self.searched_providers.append(provider)
            return iter(rows_for(query, provider))

        for target, attribute, replacement in (
            (resolve.genomepy, "search", fake_search),
            (resolve, "provider_reachable", lambda name: True),
            # The real check asks Ensembl whether the release holds the assembly.
            (resolve, "validate_release",
             lambda provider, release, assembly:
                 None if provider.lower() == "ensembl" else "not pinnable"),
        ):
            patcher = mock.patch.object(target, attribute, replacement)
            patcher.start()
            self.addCleanup(patcher.stop)

    def resolve(self, rows, cache=None,
                header="species,query,provider,annotation,aliases,release"):
        """Run the resolver on a request. Returns (exit code, resolved rows, stderr)."""
        with tempfile.TemporaryDirectory() as directory:
            raw = os.path.join(directory, "requests_raw.csv")
            out = os.path.join(directory, "requests_resolved.csv")
            with open(raw, "w", newline="") as f:
                f.write(header + "\n" + "\n".join(rows) + "\n")
            if cache is not None:
                with open(out, "w", newline="") as f:
                    f.write(cache)

            stderr = io.StringIO()
            code = 0
            with mock.patch.object(sys, "argv", ["resolve.py", raw, out, "--non-interactive"]), \
                    redirect_stdout(io.StringIO()), redirect_stderr(stderr):
                try:
                    resolve.main()
                except SystemExit as e:
                    code = e.code or 0
            with open(out, newline="") as f:
                return code, list(csv.DictReader(f)), stderr.getvalue()


class SameAssemblyTwoProviders(ResolverTestCase):
    REQUEST = [
        "Mus_musculus,GRCm39,ensembl,TRUE,mm_ens,",
        "Mus_musculus,GRCm39,gencode,TRUE,mm_gc,",
    ]

    def test_both_rows_resolve_with_their_own_aliases(self):
        code, rows, err = self.resolve(self.REQUEST)
        self.assertEqual(code, 0, err)
        self.assertEqual([(r["provider"], r["aliases"]) for r in rows],
                         [("Ensembl", "mm_ens"), ("GENCODE", "mm_gc")])

    def test_both_rows_are_answered_from_the_cache(self):
        # Both rows share (species, query, release). The cache has to tell them
        # apart by provider, or each run evicts one and searches for it again.
        cache = (RESOLVED_HEADER
                 + "Mus_musculus,GRCm39,Ensembl,TRUE,GRCm39,,mm_ens\n"
                 + "Mus_musculus,GRCm39,GENCODE,TRUE,GRCm39,,mm_gc\n")
        code, rows, err = self.resolve(self.REQUEST, cache=cache)
        self.assertEqual(code, 0, err)
        self.assertEqual(self.searched, [])
        self.assertEqual([r["provider"] for r in rows], ["Ensembl", "GENCODE"])


class TwoReleasesOfOneAssembly(ResolverTestCase):
    def test_resolve_side_by_side_and_may_share_an_alias(self):
        code, rows, err = self.resolve([
            "Homo_sapiens,GRCh38,ensembl,TRUE,human,115",
            "Homo_sapiens,GRCh38,ensembl,TRUE,human,116",
        ])
        self.assertEqual(code, 0, err)
        self.assertEqual([(r["assembly"], r["release"]) for r in rows],
                         [("GRCh38.p14", "115"), ("GRCh38.p14", "116")])

    def test_unpinned_cache_entry_does_not_answer_a_pinned_request(self):
        cache = RESOLVED_HEADER + "Homo_sapiens,STALE,Ensembl,TRUE,GRCh38,,\n"
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,ensembl,TRUE,,116"], cache=cache)
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["assembly"], "GRCh38.p14")
        self.assertEqual(self.searched, ["GRCh38"])


class DuplicateRequests(ResolverTestCase):
    def test_two_queries_for_one_download_are_rejected(self):
        code, rows, err = self.resolve([
            "Homo_sapiens,GRCh38,ensembl,TRUE,first,",
            "Homo_sapiens,human,ensembl,FALSE,second,",
        ])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(len(rows), 1)
        self.assertIn("line 3", err)                   # the row at fault
        self.assertIn("line 2 already requests", err)  # the row it duplicates

    def test_one_alias_on_two_genomes_is_rejected(self):
        code, rows, err = self.resolve([
            "Mus_musculus,GRCm39,ensembl,TRUE,mouse,",
            "Mus_musculus,GRCm39,gencode,TRUE,mouse,",
        ])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual([r["provider"] for r in rows], ["Ensembl"])
        self.assertIn("alias 'mouse' is already given", err)


class AnnotationGate(ResolverTestCase):
    def test_requested_annotation_that_does_not_exist_fails(self):
        code, rows, err = self.resolve(["Foo_bar,unannotated,ncbi,TRUE,,"])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(rows, [])
        self.assertIn("annotation was requested", err)

    def test_ucsc_list_of_flags_is_read_as_unannotated(self):
        code, rows, err = self.resolve(["Foo_baz,unannotated_ucsc,ucsc,true,,"])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("annotation was requested", err)

    def test_sequence_only_request_for_the_same_assembly_resolves(self):
        code, rows, err = self.resolve(["Foo_bar,unannotated,ncbi,FALSE,,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(len(rows), 1)

    def test_annotated_assembly_passes_and_keeps_its_aliases(self):
        code, rows, err = self.resolve(
            ["Saccharomyces_cerevisiae,sacCer3,ucsc,TRUE,yeast; sacCer,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["aliases"], "yeast;sacCer")


class ProviderOffline(ResolverTestCase):
    """
    genomepy does not fail when one provider is down: it logs a warning and
    searches the rest, so the outage shows up only as rows that are not there.
    """

    def setUp(self):
        super().setUp()
        # main() installs the listener; one test here calls search_assemblies()
        # directly, so install it explicitly instead of relying on test order.
        resolve.configure_genomepy_logging(io.StringIO(), False)

        def search_with_ensembl_down(query, provider=None, exact=False, size=False):
            # As genomepy behaves: a warning while some provider still answers, and
            # ConnectionError when none of the requested ones does - which, asked
            # for Ensembl alone, is straight away.
            if provider and provider.lower() == "ensembl":
                raise ConnectionError("ensembl appears to be offline")
            if not provider:
                resolve.genomepy_logger.warning("Ensembl appears to be offline.")
            return iter(row for row in rows_for(query, provider) if row[1] != "Ensembl")

        patcher = mock.patch.object(resolve.genomepy, "search", search_with_ensembl_down)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_request_for_the_offline_provider_is_reported_as_unreachable(self):
        # provider_reachable() answers True throughout this file: the provider is
        # back by the time anyone asks. A search takes minutes, so that is the
        # normal case, and it used to produce "no matches found" for a good request.
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,ensembl,TRUE,human,"])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(rows, [])
        self.assertIn("provider 'ensembl' is unreachable", err)
        self.assertNotIn("no matches found", err)

    def test_request_for_another_provider_is_unaffected(self):
        code, rows, err = self.resolve(["Mus_musculus,GRCm39,gencode,TRUE,mm_gc,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["provider"], "GENCODE")

    def test_unfiltered_search_says_its_candidates_may_be_incomplete(self):
        # One match is left, and it would otherwise be announced as the exact one.
        code, rows, err = self.resolve(["Mus_musculus,GRCm39,,TRUE,mm,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["provider"], "GENCODE")
        self.assertIn("ensembl could not be searched", err)

    def test_unfiltered_search_with_no_match_names_the_provider_it_missed(self):
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,,TRUE,human,"])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("among the providers that answered", err)

    def test_outage_on_one_row_does_not_stop_the_rows_after_it(self):
        # genomepy raises for the Ensembl row. Uncaught, that ended the whole run
        # with exit 1 and the GENCODE row was never looked at.
        code, rows, err = self.resolve([
            "Homo_sapiens,GRCh38,ensembl,TRUE,human,",
            "Mus_musculus,GRCm39,gencode,TRUE,mm_gc,",
        ])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual([r["provider"] for r in rows], ["GENCODE"])
        self.assertIn("provider 'ensembl' is unreachable", err)

    def test_every_provider_offline_is_reported_against_the_row(self):
        def nothing_answers(query, provider=None, exact=False, size=False):
            raise ConnectionError("Ensembl, UCSC, NCBI, GENCODE appears to be offline")

        with mock.patch.object(resolve.genomepy, "search", nothing_answers):
            code, rows, err = self.resolve(["Mus_musculus,GRCm39,,TRUE,mm,"])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(rows, [])
        self.assertIn("could not be searched", err)

    def test_offline_state_does_not_leak_into_the_next_search(self):
        _, offline = resolve.search_assemblies("GRCm39")
        self.assertEqual(offline, {"ensembl"})
        with mock.patch.object(resolve.genomepy, "search",
                               lambda query, provider=None, size=False: iter(rows_for(query))):
            _, offline = resolve.search_assemblies("GRCm39")
        self.assertEqual(offline, set())


class SearchScope(ResolverTestCase):
    def test_row_naming_a_provider_searches_only_that_provider(self):
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,Ensembl,TRUE,human,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(self.searched_providers, ["ensembl"])

    def test_row_naming_no_provider_searches_all_of_them(self):
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,,TRUE,human,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(self.searched_providers, [None])

    def test_unknown_provider_is_refused_before_any_search(self):
        # 'ensemble' used to search everything, match nothing, and be reported as
        # "no matches found" - a typo in one column blamed on another.
        code, rows, err = self.resolve([
            "Homo_sapiens,GRCh38,ensemble,TRUE,human,",
            "Mus_musculus,GRCm39,gencode,TRUE,mm_gc,",
        ])
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("unknown provider 'ensemble'", err)
        self.assertIn("ensembl", err)                      # the known names are listed
        self.assertEqual(self.searched_providers, ["gencode"])
        self.assertEqual([r["provider"] for r in rows], ["GENCODE"])


class GenomeByAddress(ResolverTestCase):
    """Provider 'url': the row is the request, and nothing is searched."""

    HEADER = "species,query,provider,annotation,aliases,release,url,annotation_url"
    FASTA = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/GRCh38_no_alt_analysis_set.fna.gz"
    GTF = "https://ftp.ebi.ac.uk/pub/databases/gencode/release_49/gencode.v49.annotation.gtf.gz"

    def by_address(self, *rows):
        return self.resolve(list(rows), header=self.HEADER)

    def test_row_is_taken_as_written_without_a_search(self):
        code, rows, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,false,hs_noalt,,{self.FASTA},")
        self.assertEqual(code, 0, err)
        self.assertEqual(self.searched, [])
        row = rows[0]
        # The query is the genome's name; the provider is spelled as genomepy does.
        self.assertEqual((row["assembly"], row["provider"]), ("GRCh38_no_alt", "URL"))
        self.assertEqual((row["url"], row["annotation_url"]), (self.FASTA, ""))
        self.assertEqual(row["aliases"], "hs_noalt")

    def test_annotation_is_given_by_address_too(self):
        code, rows, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,TRUE,,,{self.FASTA},{self.GTF}")
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["annotation_url"], self.GTF)

    def test_annotation_without_an_address_is_refused(self):
        code, rows, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,TRUE,,,{self.FASTA},")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("annotation_url", err)

    def test_annotation_address_must_be_a_format_genomepy_takes(self):
        code, _, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,TRUE,,,{self.FASTA},https://example.org/genes.txt")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("not a GTF, GFF or BED", err)

    def test_missing_address_is_refused(self):
        code, _, err = self.by_address("Homo_sapiens,GRCh38_no_alt,url,false,,,,")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("url column is empty", err)

    def test_address_with_shell_characters_is_refused(self):
        # The address reaches a shell command. Each of these would do something
        # there other than name a file.
        for address in ("https://example.org/a.fa.gz;rm", "https://example.org/$(id).fa.gz",
                        "https://example.org/a.fa.gz?x=1&y=2", "file:///etc/passwd",
                        "https://example.org/a b.fa.gz"):
            code, rows, err = self.by_address(
                f"Homo_sapiens,GRCh38_no_alt,url,false,,,{address},")
            self.assertEqual(code, resolve.EXIT_UNRESOLVED, address)
            self.assertEqual(rows, [], address)
            self.assertIn("cannot be used", err, address)

    def test_name_that_could_not_be_a_directory_is_refused(self):
        code, _, err = self.by_address(f"Homo_sapiens,GRCh38/no alt,url,false,,,{self.FASTA},")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("the name the genome will be known by", err)

    def test_release_cannot_be_pinned(self):
        code, _, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,false,,116,{self.FASTA},")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("cannot be pinned", err)

    def test_address_on_a_row_for_another_provider_is_refused_not_ignored(self):
        code, rows, err = self.by_address(
            f"Homo_sapiens,GRCh38,ensembl,TRUE,human,,{self.FASTA},")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(rows, [])
        self.assertEqual(self.searched, [])
        self.assertIn("only used with provider 'url'", err)

    def test_same_name_twice_is_a_duplicate(self):
        code, rows, err = self.by_address(
            f"Homo_sapiens,GRCh38_no_alt,url,false,a,,{self.FASTA},",
            f"Homo_sapiens,GRCh38_no_alt,url,false,b,,{self.FASTA}2,")
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(len(rows), 1)
        self.assertIn("already requests", err)

    def test_rows_for_other_providers_carry_empty_address_columns(self):
        code, rows, err = self.by_address("Homo_sapiens,GRCh38,ensembl,TRUE,human,,,")
        self.assertEqual(code, 0, err)
        self.assertEqual((rows[0]["url"], rows[0]["annotation_url"]), ("", ""))


class ResourcesColumn(ResolverTestCase):
    HEADER = "species,query,provider,annotation,aliases,release,url,annotation_url,resources"

    def test_request_is_checked_against_the_registry_and_carried_through(self):
        code, rows, err = self.resolve(
            ["Homo_sapiens,GRCh38,ensembl,TRUE,human,,,, clinvar:GRCh38 "], header=self.HEADER)
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["resources"], "clinvar:GRCh38")

    def test_unknown_resource_fails_the_row_before_any_search(self):
        code, rows, err = self.resolve(
            ["Homo_sapiens,GRCh38,ensembl,TRUE,human,,,,gnomad:GRCh38"], header=self.HEADER)
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertEqual(rows, [])
        self.assertEqual(self.searched, [])
        self.assertIn("unknown resource 'gnomad'", err)

    def test_build_the_resource_is_not_published_for_fails_the_row(self):
        code, _, err = self.resolve(
            ["Homo_sapiens,GRCh38,ensembl,TRUE,human,,,,blacklist:GRCh38"], header=self.HEADER)
        self.assertEqual(code, resolve.EXIT_UNRESOLVED)
        self.assertIn("not published for 'GRCh38'", err)
        self.assertIn("hg38", err)

    def test_rows_without_the_column_resolve_with_it_empty(self):
        code, rows, err = self.resolve(["Homo_sapiens,GRCh38,ensembl,TRUE,human,"])
        self.assertEqual(code, 0, err)
        self.assertEqual(rows[0]["resources"], "")


class LogStyling(unittest.TestCase):
    RED = "\x1b[31m"

    def lines(self, colorize):
        stream = io.StringIO()
        resolve.configure_genomepy_logging(stream=stream, colorize=colorize)
        self.addCleanup(resolve.configure_genomepy_logging, io.StringIO(), False)
        resolve.genomepy_logger.info("Downloading assembly summaries from UCSC")
        resolve.genomepy_logger.warning("Ensembl appears to be offline.")
        return stream.getvalue().splitlines()

    def test_warning_is_red_and_progress_lines_are_not(self):
        info, warning = self.lines(colorize=True)
        self.assertIn(self.RED, warning)
        self.assertNotIn(self.RED, info)
        # Level and message both: the message is the part that says what happened.
        self.assertGreaterEqual(warning.count(self.RED), 2)

    def test_no_escape_codes_when_output_is_not_a_terminal(self):
        info, warning = self.lines(colorize=False)
        self.assertNotIn("\x1b", info + warning)
        self.assertTrue(warning.endswith("| WARNING | Ensembl appears to be offline."), warning)
        self.assertTrue(info.endswith("| INFO | Downloading assembly summaries from UCSC"), info)


class OlderRequestFiles(ResolverTestCase):
    def test_request_without_aliases_or_release_columns_resolves(self):
        code, rows, err = self.resolve(["Saccharomyces_cerevisiae,sacCer3,ucsc,false"],
                                       header="species,query,provider,annotation")
        self.assertEqual(code, 0, err)
        self.assertEqual((rows[0]["aliases"], rows[0]["release"]), ("", ""))


if __name__ == "__main__":
    unittest.main()
