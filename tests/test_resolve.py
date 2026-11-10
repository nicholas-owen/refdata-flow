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


class ResolverTestCase(unittest.TestCase):
    def setUp(self):
        self.searched = []

        def fake_search(query, size=False):
            self.searched.append(query)
            return iter(SEARCH_RESULTS.get(query, []))

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


class OlderRequestFiles(ResolverTestCase):
    def test_request_without_aliases_or_release_columns_resolves(self):
        code, rows, err = self.resolve(["Saccharomyces_cerevisiae,sacCer3,ucsc,false"],
                                       header="species,query,provider,annotation")
        self.assertEqual(code, 0, err)
        self.assertEqual((rows[0]["aliases"], rows[0]["release"]), ("", ""))


if __name__ == "__main__":
    unittest.main()
