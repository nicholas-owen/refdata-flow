# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of third-party resources: the registry, the fetch step and the gate.

Nothing is fetched. The fetch tests replace fetch_resource.open_url, the one
function that touches the network, with canned responses. The gate tests are
given small VCF and BED files written here and a chrom.sizes to read them against.
tests/test_nextflow_stub.py runs the three steps together inside the workflow.
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

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.insert(0, os.path.join(REPO, "bin"))
import fetch_resource  # noqa: E402
import prepare_resource  # noqa: E402
import resource_registry  # noqa: E402

# A genome named the Ensembl way.
SIZES = {"1": 1000, "2": 800, "MT": 100}

VCF_HEAD = "##fileformat=VCFv4.1\n##reference=toy\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def vcf(*records, contig_lines=""):
    return VCF_HEAD.replace("#CHROM", contig_lines + "#CHROM") + "".join(
        f"{contig}\t{position}\t.\tA\tC\t.\t.\t.\n" for contig, position in records)


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = directory.name

    def path(self, *parts):
        return os.path.join(self.directory, *parts)

    def write_gz(self, name, text):
        with gzip.open(self.path(name), "wt") as f:
            f.write(text)
        return self.path(name)

    def write_json(self, name, value):
        with open(self.path(name), "w") as f:
            json.dump(value, f)
        return self.path(name)


class Registry(TempDirTestCase):
    def test_every_shipped_resource_loads_and_has_a_recipe_that_fits_it(self):
        registry = resource_registry.load_registry()
        self.assertTrue({"clinvar", "blacklist"} <= set(registry))
        for name, entry in registry.items():
            # Registered under its own name, by a recipe of that name, whose
            # files are what bin/prepare_resource.py hands over.
            with open(os.path.join(REPO, "recipes", f"{name}.json")) as f:
                recipe = json.load(f)
            self.assertEqual(recipe["name"], name)
            keys = {item["key"] for item in recipe["required_files"]}
            expected = {name, "provenance"} | ({"tbi"} if entry["index"] == "tabix" else set())
            self.assertEqual(keys, expected, name)
            self.assertIn(name, recipe["assets"])

    def test_request_is_parsed_into_name_and_build(self):
        registry = resource_registry.load_registry()
        pairs, problem = resource_registry.parse_request(" clinvar:GRCh38 ; blacklist:hg38;", registry)
        self.assertIsNone(problem)
        self.assertEqual(pairs, [("clinvar", "GRCh38"), ("blacklist", "hg38")])
        self.assertEqual(resource_registry.format_request(pairs), "clinvar:GRCh38;blacklist:hg38")
        self.assertEqual(resource_registry.parse_request("", registry), ([], None))

    def test_unusable_requests_say_what_is_available(self):
        registry = resource_registry.load_registry()
        for request, expected in (
            ("dbsnp:GRCh38", "unknown resource 'dbsnp' (known: blacklist, clinvar"),
            ("clinvar:hg38", "not published for 'hg38' (it has: GRCh37, GRCh38)"),
            ("clinvar", "must be written as name:build"),
            ("clinvar:GRCh38;clinvar:GRCh38", "listed twice"),
            ("clinvar:GRCh38;clinvar:GRCh37", "two builds on one row"),
        ):
            pairs, problem = resource_registry.parse_request(request, registry)
            self.assertEqual(pairs, [], request)
            self.assertIn(expected, problem, request)

    def test_a_faulty_definition_is_named_with_its_fault(self):
        good = {"name": "thing", "kind": "bed", "index": None,
                "variants": {"hg38": {"url": "https://example.org/x.bed.gz", "version": "v1"}}}
        self.assertEqual(resource_registry.load_entry(self.write_json("thing.json", good))["name"],
                         "thing")
        for change, expected in (
            ({"kind": "bam"}, "'kind' must be one of"),
            ({"index": "csi"}, "'index' must be"),
            ({"variants": {}}, "at least one entry under 'variants'"),
            ({"variants": {"hg38": {"url": "ftp://x/y", "version": "v1"}}}, "http(s) 'url'"),
            ({"variants": {"hg38": {"url": "https://x/y"}}}, "needs a 'version'"),
            ({"name": "other"}, "must be called 'other.json'"),
        ):
            with self.assertRaises(ValueError) as caught:
                resource_registry.load_entry(self.write_json("thing.json", dict(good, **change)))
            self.assertIn(expected, str(caught.exception))


class Gate(TempDirTestCase):
    def check(self, text, kind="vcf", min_matched=0.95, name="r.gz"):
        found = prepare_resource.inspect(self.write_gz(name, text), kind)
        return prepare_resource.gate(found, SIZES, min_matched)

    def test_resource_on_the_genomes_own_sequences_passes(self):
        reasons, summary = self.check(vcf(("1", 10), ("2", 700), ("MT", 5)))
        self.assertEqual(reasons, [])
        self.assertEqual((summary["records"], summary["share_matched"]), (3, 1.0))
        self.assertEqual(summary["contigs_absent_from_genome"], {})

    def test_differently_named_resource_is_refused(self):
        # The case that passes silently everywhere else: chr1 beside 1.
        reasons, summary = self.check(vcf(("chr1", 10), ("chr2", 700)))
        self.assertEqual(summary["share_matched"], 0.0)
        self.assertEqual(len(reasons), 1)
        self.assertIn("only 0.0% of its 2 records", reasons[0])
        self.assertIn("chr1", reasons[0])

    def test_right_names_on_the_wrong_build_are_refused_by_their_coordinates(self):
        reasons, summary = self.check(vcf(("1", 10), ("1", 1001)))
        self.assertEqual(summary["share_matched"], 1.0)
        self.assertEqual(len(reasons), 1)
        self.assertIn("beyond the end", reasons[0])
        self.assertIn("1 reaches 1,001", reasons[0])

    def test_last_base_of_a_sequence_is_not_beyond_its_end(self):
        reasons, _ = self.check(vcf(("1", 1000)))
        self.assertEqual(reasons, [])

    def test_lengths_the_file_declares_must_be_the_genomes(self):
        reasons, _ = self.check(vcf(("1", 10), contig_lines="##contig=<ID=1,length=999>\n"))
        self.assertEqual(len(reasons), 1)
        self.assertIn("declares different lengths", reasons[0])
        reasons, _ = self.check(vcf(("1", 10), contig_lines="##contig=<ID=1,length=1000>\n"))
        self.assertEqual(reasons, [])

    def test_a_few_records_on_sequences_the_genome_lacks_pass_and_are_reported(self):
        records = [("1", n) for n in range(1, 100)] + [("NT_187633.1", 5)]
        reasons, summary = self.check(vcf(*records))
        self.assertEqual(reasons, [])
        self.assertEqual(summary["contigs_absent_from_genome"], {"NT_187633.1": 1})
        self.assertEqual(summary["records_on_genome_contigs"], 99)

    def test_threshold_is_the_resources_own(self):
        records = [("1", 5)] * 9 + [("other", 5)]
        self.assertEqual(self.check(vcf(*records), min_matched=0.9)[0], [])
        self.assertEqual(len(self.check(vcf(*records), min_matched=1.0)[0]), 1)

    def test_empty_resource_is_refused(self):
        reasons, _ = self.check(VCF_HEAD)
        self.assertEqual(reasons, ["it contains no records"])

    def test_bed_intervals_are_read_by_their_end_and_free_text_columns_survive(self):
        bed = "track name=x\n1\t0\t500\tLow Mappability\n2\t10\t800\tHigh Signal Region\n"
        reasons, summary = self.check(bed, kind="bed", min_matched=1.0)
        self.assertEqual(reasons, [])
        self.assertEqual(summary["records"], 2)
        reasons, _ = self.check("2\t10\t801\tx\n", kind="bed", min_matched=1.0)
        self.assertIn("beyond the end", reasons[0])


class Prepare(TempDirTestCase):
    """main(): from fetch_resource.py's output directory to something registrable."""

    def setUp(self):
        super().setUp()
        with open(self.path("chrom.sizes"), "w") as f:
            f.write("".join(f"{name}\t{length}\n" for name, length in SIZES.items()))

    def prepare(self, name, kind, text, index):
        self.write_json(f"{name}.json", {
            "name": name, "kind": kind, "index": index, "min_matched": 0.95,
            "variants": {"toy": {"url": "https://example.org/x", "version": "v1",
                                 "naming": "Ensembl-style"}}})
        os.makedirs(self.path("fetched"), exist_ok=True)
        filename = f"published.{kind}.gz"
        with gzip.open(self.path("fetched", filename), "wt") as f:
            f.write(text)
        self.write_json(os.path.join("fetched", f"{name}.fetch.json"),
                        {"resource": name, "variant": "toy", "version": "v1", "file": filename,
                         "file_sha256": "as-published"})
        argv = ["prepare_resource.py", "--registry", self.path(f"{name}.json"), "--variant", "toy",
                "--genome", "toy_ensembl", "--fetched", self.path("fetched"),
                "--chrom-sizes", self.path("chrom.sizes"), "--outdir", self.path("prepared")]
        out, err, code = io.StringIO(), io.StringIO(), 0
        with mock.patch.object(sys, "argv", argv), redirect_stdout(out), redirect_stderr(err):
            try:
                prepare_resource.main()
            except SystemExit as e:
                code = e.code or 0
        return code, out.getvalue(), err.getvalue()

    def test_indexed_resource_comes_out_bgzipped_with_a_tabix_index_that_works(self):
        import pysam
        code, out, err = self.prepare("thing", "vcf", vcf(("1", 10), ("1", 500), ("2", 7)), "tabix")
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(os.listdir(self.path("prepared"))),
                         ["TAG", "provenance.json", "thing.data.gz", "thing.data.gz.tbi"])
        data = self.path("prepared", "thing.data.gz")
        self.assertTrue(prepare_resource.is_bgzf(data))
        with pysam.TabixFile(data) as index:
            self.assertEqual(len(list(index.fetch("1", 0, 100))), 1)
            self.assertEqual(len(list(index.fetch("1"))), 2)
        with open(self.path("prepared", "TAG")) as f:
            self.assertEqual(f.read().strip(), "v1")
        with open(self.path("prepared", "provenance.json")) as f:
            record = json.load(f)
        self.assertEqual(record["file_sha256"], "as-published")   # the fetch record is kept
        self.assertEqual(record["genome"], "toy_ensembl")
        self.assertTrue(record["indexed"])
        # It arrived as plain gzip and had to be recompressed, and says so.
        self.assertFalse(record["stored_as_published"])
        self.assertEqual(record["compatibility"]["records"], 3)
        self.assertIn("##reference=toy", record["file_header"])

    def test_unindexed_resource_is_stored_byte_for_byte_as_published(self):
        code, _, err = self.prepare("thing", "bed", "1\t0\t500\tLow Mappability\n", None)
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(os.listdir(self.path("prepared"))),
                         ["TAG", "provenance.json", "thing.data.gz"])
        with open(self.path("fetched", "published.bed.gz"), "rb") as a, \
                open(self.path("prepared", "thing.data.gz"), "rb") as b:
            self.assertEqual(a.read(), b.read())
        with open(self.path("prepared", "provenance.json")) as f:
            self.assertTrue(json.load(f)["stored_as_published"])

    def test_resource_that_does_not_fit_is_refused_and_nothing_is_written(self):
        code, _, err = self.prepare("thing", "bed", "chr1\t0\t500\tx\n", None)
        self.assertEqual(code, 1)
        self.assertIn("was NOT registered for toy_ensembl", err)
        self.assertIn("Ensembl-style", err)              # what the registry says it fits
        self.assertFalse(os.path.exists(self.path("prepared")))


class FakeResponse:
    def __init__(self, body=b"", headers=None):
        self._body, self.headers = io.BytesIO(body), headers or {}

    def read(self, size=-1):
        return self._body.read(size)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class Fetch(TempDirTestCase):
    URL = "https://example.org/pub/thing.vcf.gz"
    BODY = b"not really a vcf, but bytes are bytes" * 100
    MODIFIED = {"Last-Modified": "Mon, 09 Nov 2026 15:02:07 GMT"}

    def fetch(self, version, held="", checksum=None, responses=None, extra=()):
        variant = {"url": self.URL, "version": version}
        if checksum:
            variant["checksum"] = {"type": "md5", "url": self.URL + ".md5"}
        registry = self.write_json("thing.json", {
            "name": "thing", "kind": "vcf", "index": None, "terms": "https://example.org/terms",
            "variants": {"toy": variant}})
        responses = responses if responses is not None else {
            ("HEAD", self.URL): FakeResponse(headers=self.MODIFIED),
            ("GET", self.URL): FakeResponse(self.BODY, self.MODIFIED),
            ("GET", self.URL + ".md5"): FakeResponse(
                (checksum or "").encode() + b"  /some/path/thing_20261108.vcf.gz\n"),
        }
        self.calls = []

        def fake_open(url, method="GET"):
            self.calls.append((method, url))
            response = responses.get((method, url))
            if response is None:
                raise OSError(f"no route to {url}")
            return response

        argv = ["fetch_resource.py", "--registry", registry, "--variant", "toy",
                "--genome", "toy_ensembl", "--held", held, "--outdir", self.path("fetched"), *extra]
        out, err, code = io.StringIO(), io.StringIO(), 0
        with mock.patch.object(fetch_resource, "open_url", fake_open), \
                mock.patch.object(sys, "argv", argv), redirect_stdout(out), redirect_stderr(err):
            try:
                fetch_resource.main()
            except SystemExit as e:
                code = e.code or 0
        self.stdout, self.stderr = out.getvalue(), err.getvalue()
        return code, sorted(os.listdir(self.path("fetched")))

    def record(self):
        with open(self.path("fetched", "thing.fetch.json")) as f:
            return json.load(f)

    def test_fixed_version_already_held_is_skipped_without_touching_the_network(self):
        code, files = self.fetch("v2", held="v2")
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertEqual(self.calls, [])

    def test_dated_version_already_held_costs_one_head_request(self):
        code, files = self.fetch("last_modified", held="20261102,20261109")
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertEqual(self.calls, [("HEAD", self.URL)])

    def test_newer_dated_version_is_fetched_and_tagged_with_its_date(self):
        code, files = self.fetch("last_modified", held="20261102")
        self.assertEqual((code, files), (0, ["thing.fetch.json", "thing.vcf.gz"]))
        record = self.record()
        self.assertEqual(record["version"], "20261109")
        self.assertEqual(record["file_sha256"], hashlib.sha256(self.BODY).hexdigest())
        self.assertEqual(record["file_bytes"], len(self.BODY))
        self.assertEqual(record["url"], self.URL)
        self.assertTrue(record["retrieved_utc"].endswith("+00:00"))
        with open(self.path("fetched", "thing.vcf.gz"), "rb") as f:
            self.assertEqual(f.read(), self.BODY)

    def test_published_checksum_is_verified(self):
        code, _ = self.fetch("last_modified", checksum=hashlib.md5(self.BODY).hexdigest())
        self.assertEqual(code, 0, self.stderr)
        self.assertIn("md5 matches", self.record()["checksum_verified"])

    def test_checksum_mismatch_fails_and_discards_the_file(self):
        code, files = self.fetch("last_modified", checksum="0" * 32)
        self.assertEqual(code, 1)
        self.assertEqual(files, [])
        self.assertIn("does not match the publisher's md5", self.stderr)

    def test_no_checksum_published_is_recorded_as_such(self):
        self.fetch("v2")
        self.assertIn("no checksum is published", self.record()["checksum_verified"])
        self.assertEqual(self.record()["version"], "v2")

    def test_unreachable_publisher_never_causes_a_refetch_of_something_held(self):
        code, files = self.fetch("last_modified", held="20261102", responses={})
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertIn("could not be checked", self.stdout)

    def test_undated_head_falls_back_to_the_downloads_own_date(self):
        code, _ = self.fetch("last_modified", responses={
            ("HEAD", self.URL): FakeResponse(),
            ("GET", self.URL): FakeResponse(self.BODY, self.MODIFIED)})
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual(self.record()["version"], "20261109")

    def test_release_check_off_keeps_whatever_is_held(self):
        code, files = self.fetch("last_modified", held="20200101", extra=["--no-release-check"])
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertEqual(self.calls, [])

    def test_empty_download_is_an_error(self):
        code, files = self.fetch("v2", responses={("GET", self.URL): FakeResponse(b"")})
        self.assertEqual(code, 1)
        self.assertEqual(files, [])


if __name__ == "__main__":
    unittest.main()
