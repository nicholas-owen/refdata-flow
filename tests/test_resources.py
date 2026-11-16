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
import shutil
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
            ("gnomad:GRCh38", "unknown resource 'gnomad' (known: blacklist, clinvar"),
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


class Sources(TempDirTestCase):
    """What a variant, or a database, may say about where its file comes from."""

    def test_every_shipped_database_loads(self):
        databases = resource_registry.load_databases()
        self.assertTrue({"swissprot", "pfam", "silva", "kraken2_standard"} <= set(databases))
        self.assertEqual(databases["kraken2_standard"]["unpack"], "tar")
        self.assertIsNone(databases["swissprot"]["unpack"])

    def test_a_resource_that_follows_the_annotation_says_so_in_its_address(self):
        cache = resource_registry.load_registry()["vep_cache"]
        for label, variant in cache["variants"].items():
            self.assertEqual(variant["version"], "annotation_release", label)
            self.assertIn("{release}", variant["url"], label)

    def test_a_faulty_source_is_named_with_its_fault(self):
        def resource(**variant):
            return {"name": "thing", "kind": "vcf", "index": None, "variants": {"hg38": dict(
                {"url": "https://example.org/x.vcf.gz", "version": "v1"}, **variant)}}

        for entry, expected in (
            (resource(version={"url": "https://example.org/rel"}), "needs a 'url' and a 'pattern'"),
            (resource(version={"url": "https://example.org/rel", "pattern": "Release"}),
             "needs one group"),
            (resource(version="annotation_release"), "go together"),
            (resource(url="https://example.org/{release}/x.vcf.gz"), "go together"),
            (resource(checksum={"type": "crc"}), "a 'checksum' needs a 'type'"),
            (resource(checksum={"type": "md5"}), "either a 'url' or a 'header'"),
            (resource(checksum={"type": "md5", "url": "https://x/y", "header": "h"}),
             "either a 'url' or a 'header'"),
            (resource(rename={"report": "report.txt"}), "'rename' needs a 'report'"),
            (dict(resource(rename={"report": "https://example.org/r.txt"}), kind="vep_cache"),
             "cannot be renamed"),
            (dict(resource(), kind="vep_cache", index="tabix"), "cannot have an 'index'"),
        ):
            with self.assertRaises(ValueError) as caught:
                resource_registry.load_entry(self.write_json("thing.json", entry))
            self.assertIn(expected, str(caught.exception))

        good = {"name": "thing", "url": "https://example.org/x.fa.gz", "version": "1.0"}
        self.assertEqual(
            resource_registry.load_database(self.write_json("thing.json", good))["name"], "thing")
        for change, expected in (
            ({"unpack": "zip"}, "'unpack' must be 'tar' or null"),
            ({"version": "annotation_release", "url": "https://example.org/{release}/x"},
             "only for a resource"),
            ({"rename": {"report": "https://example.org/r.txt"}}, "'rename' is only for a resource"),
            ({"url": "x.fa.gz"}, "http(s) 'url'"),
        ):
            with self.assertRaises(ValueError) as caught:
                resource_registry.load_database(self.write_json("thing.json", dict(good, **change)))
            self.assertIn(expected, str(caught.exception))


class FetchSources(TempDirTestCase):
    """The ways of versioning and verifying a file that came after the first two."""

    URL = "https://example.org/pub/thing.vcf.gz"
    BODY = b"bytes are bytes " * 200

    def fetch(self, variant, responses, held="", extra=(), kind="vcf"):
        registry = self.write_json("thing.json", {
            "name": "thing", "kind": kind, "index": None,
            "variants": {"toy": dict({"url": self.URL, "version": "v1"}, **variant)}})
        self.calls = []

        def fake_open(url, method="GET"):
            self.calls.append((method, url))
            if (method, url) not in responses:
                raise OSError(f"no route to {url}")
            return responses[(method, url)]

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
        files = sorted(os.listdir(self.path("fetched"))) if os.path.isdir(self.path("fetched")) else []
        return code, files

    def record(self):
        with open(self.path("fetched", "thing.fetch.json")) as f:
            return json.load(f)

    # -- a release the publisher states in a file of its own --------------------------

    RELEASE = {"version": {"url": "https://example.org/pub/reldate.txt",
                           "pattern": r"Thing Release (\S+)"}}

    def test_release_stated_in_a_file_is_the_version(self):
        code, files = self.fetch(self.RELEASE, {
            ("GET", "https://example.org/pub/reldate.txt"):
                FakeResponse(b"Knowledgebase Release 2026_04 consists of:\nThing Release 2026_04 of 04-Nov\n"),
            ("GET", self.URL): FakeResponse(self.BODY)})
        self.assertEqual((code, files), (0, ["thing.fetch.json", "thing.vcf.gz"]), self.stderr)
        self.assertEqual(self.record()["version"], "2026_04")
        self.assertIn("reldate.txt", self.record()["version_from"])

    def test_release_file_may_be_gzipped_and_a_held_release_is_not_fetched(self):
        code, files = self.fetch(self.RELEASE, {
            ("GET", "https://example.org/pub/reldate.txt"):
                FakeResponse(gzip.compress(b"Thing Release 39.0\n"))}, held="38.1,39.0")
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertEqual(self.calls, [("GET", "https://example.org/pub/reldate.txt")])

    def test_release_that_cannot_be_read_never_refetches_something_held(self):
        code, files = self.fetch(self.RELEASE, {}, held="38.1")
        self.assertEqual((code, files), (0, ["SKIP"]))

    # -- a resource published per release of the genome's annotation -----------------

    PER_RELEASE = {"url": "https://example.org/pub/release-{release}/thing_{release}_toy.tar.gz",
                   "version": "annotation_release", "annotation": "ensembl_gtf"}
    AT_116 = "https://example.org/pub/release-116/thing_116_toy.tar.gz"

    def test_version_and_address_follow_the_genomes_annotation_release(self):
        code, files = self.fetch(self.PER_RELEASE, {("GET", self.AT_116): FakeResponse(self.BODY)},
                                 extra=["--annotation", "ensembl_gtf:116"], kind="vep_cache")
        self.assertEqual((code, files), (0, ["thing.fetch.json", "thing_116_toy.tar.gz"]),
                         self.stderr)
        record = self.record()
        self.assertEqual((record["version"], record["url"]), ("116", self.AT_116))
        # Nothing was asked of the publisher to find the version out.
        self.assertEqual(self.calls, [("GET", self.AT_116)])

    def test_release_already_held_is_skipped_and_another_release_is_not(self):
        code, files = self.fetch(self.PER_RELEASE, {}, held="110,116", kind="vep_cache",
                                 extra=["--annotation", "ensembl_gtf:116"])
        self.assertEqual((code, files, self.calls), (0, ["SKIP"], []))
        code, files = self.fetch(self.PER_RELEASE, {("GET", self.AT_116): FakeResponse(self.BODY)},
                                 held="110", kind="vep_cache",
                                 extra=["--annotation", "ensembl_gtf:116"])
        self.assertEqual(code, 0, self.stderr)
        self.assertIn("thing_116_toy.tar.gz", files)

    def test_release_the_publisher_has_no_file_for_is_refused_and_not_retried(self):
        # Seen for real: yeast's annotation is at Ensembl Genomes release 63, and
        # the address built from it did not exist. A 404 is an answer, not a fault.
        import urllib.error
        missing = urllib.error.HTTPError(self.AT_116, 404, "Not Found", {}, None)

        class Missing:
            def __enter__(self):
                raise missing

            def __exit__(self, *exc):
                return False

        code, files = self.fetch(self.PER_RELEASE, {("GET", self.AT_116): Missing()},
                                 kind="vep_cache", extra=["--annotation", "ensembl_gtf:116"])
        self.assertEqual((code, files), (fetch_resource.EXIT_NO_RETRY, []))
        self.assertIn("has no file at", self.stderr)
        self.assertIn("release 116 of this genome's annotation", self.stderr)

        # A server that is merely unwell is a different matter, and is retried.
        unwell = urllib.error.HTTPError(self.AT_116, 503, "Unavailable", {}, None)
        missing = unwell
        code, _ = self.fetch(self.PER_RELEASE, {("GET", self.AT_116): Missing()},
                             kind="vep_cache", extra=["--annotation", "ensembl_gtf:116"])
        self.assertEqual(code, 1)

    def test_genome_without_that_annotation_is_refused_and_not_retried(self):
        for annotation, said in (("", "without an annotation"),
                                 ("ucsc_gtf:20210210", "with 'ucsc_gtf:20210210'")):
            code, files = self.fetch(self.PER_RELEASE, {}, kind="vep_cache",
                                     extra=["--annotation", annotation])
            self.assertEqual((code, files), (fetch_resource.EXIT_NO_RETRY, []))
            self.assertIn(said, self.stderr)
            self.assertIn("'ensembl_gtf'", self.stderr)
            self.assertEqual(self.calls, [])

    # -- a checksum sent with the download, and one in a list of many ---------------

    def goog(self, body):
        import base64
        return {"x-goog-hash": "crc32c=AAAAAA==,md5=" + base64.b64encode(
            hashlib.md5(body).digest()).decode()}

    HEADER_CHECKSUM = {"checksum": {"type": "md5", "header": "x-goog-hash"}}

    def test_checksum_sent_with_the_download_is_verified(self):
        code, _ = self.fetch(self.HEADER_CHECKSUM,
                             {("GET", self.URL): FakeResponse(self.BODY, self.goog(self.BODY))})
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual(self.record()["checksum_verified"],
                         "md5 matches the publisher's x-goog-hash header")

    def test_download_that_differs_from_its_own_header_is_discarded(self):
        code, files = self.fetch(self.HEADER_CHECKSUM,
                                 {("GET", self.URL): FakeResponse(self.BODY, self.goog(b"other"))})
        self.assertEqual((code, files), (1, []))
        self.assertIn("does not match the publisher's md5", self.stderr)

    def test_header_without_that_checksum_is_recorded_not_failed(self):
        code, _ = self.fetch(self.HEADER_CHECKSUM, {
            ("GET", self.URL): FakeResponse(self.BODY, {"x-goog-hash": "crc32c=AAAAAA=="})})
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual(self.record()["checksum_verified"],
                         "the publisher sent no md5 with this file")

    def test_checksum_is_read_from_the_line_that_names_the_file(self):
        listing = ("0" * 32 + "  thing.vcf.dat.gz\n" + hashlib.md5(self.BODY).hexdigest()
                   + "  thing.vcf.gz\n" + "1" * 32 + "  other.gz\n").encode()
        code, _ = self.fetch({"checksum": {"type": "md5", "url": "https://example.org/pub/md5s"}}, {
            ("GET", self.URL): FakeResponse(self.BODY),
            ("GET", "https://example.org/pub/md5s"): FakeResponse(listing)})
        self.assertEqual(code, 0, self.stderr)
        self.assertIn("md5 matches", self.record()["checksum_verified"])

    # -- the size the publisher announces --------------------------------------------

    def test_file_bigger_than_the_space_free_is_refused_before_a_byte_is_written(self):
        response = FakeResponse(self.BODY, {"Content-Length": str(len(self.BODY))})
        with mock.patch.object(fetch_resource.shutil, "disk_usage",
                               return_value=mock.Mock(free=len(self.BODY) - 1)):
            code, files = self.fetch({}, {("GET", self.URL): response})
        self.assertEqual((code, files), (fetch_resource.EXIT_NO_RETRY, []))
        self.assertIn("Nothing was downloaded", self.stderr)

    def test_transfer_shorter_than_announced_fails_and_is_discarded(self):
        response = FakeResponse(self.BODY, {"Content-Length": str(len(self.BODY) + 1)})
        code, files = self.fetch({}, {("GET", self.URL): response})
        self.assertEqual((code, files), (1, []))
        self.assertIn("cut short", self.stderr)

    def test_table_of_names_is_fetched_with_a_resource_that_must_be_renamed(self):
        report = "https://example.org/genomes/toy_assembly_report.txt"
        code, files = self.fetch({"rename": {"report": report}}, {
            ("GET", self.URL): FakeResponse(self.BODY),
            ("GET", report): FakeResponse(b"# Sequence-Name\n1\n")})
        self.assertEqual((code, files), (0, ["thing.fetch.json", "thing.names.txt", "thing.vcf.gz"]),
                         self.stderr)
        names = self.record()["names_report"]
        self.assertEqual((names["url"], names["file"]), (report, "thing.names.txt"))
        self.assertEqual(names["file_sha256"], hashlib.sha256(b"# Sequence-Name\n1\n").hexdigest())


REPORT_HEAD = ("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
               "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\t"
               "UCSC-style-name\n")


def report(*rows):
    """An NCBI assembly report: rows of (name, GenBank, RefSeq, length, UCSC name)."""
    return "# Assembly name:  toy\n" + REPORT_HEAD + "".join(
        f"{name}\tassembled-molecule\t{name}\tChromosome\t{genbank}\t=\t{refseq}\t"
        f"Primary Assembly\t{length}\t{ucsc}\n"
        for name, genbank, refseq, length, ucsc in rows)


# The toy genome of SIZES, as NCBI would report it, plus a patch it does not have.
REPORT = report(("1", "CM1.1", "NC_1.1", 1000, "chr1"), ("2", "CM2.1", "NC_2.1", 800, "chr2"),
                ("MT", "J01415.2", "NC_012920.1", 100, "chrM"),
                ("HG1_PATCH", "KQ1.1", "NW_9.1", 50, "na"))


class PrepareCase(TempDirTestCase):
    """main(), given a registry entry and whatever the fetch step left behind."""

    def run_prepare(self, entry, files, record=None, sizes=SIZES):
        with open(self.path("chrom.sizes"), "w") as f:
            f.write("".join(f"{name}\t{length}\n" for name, length in sizes.items()))
        entry = dict({"name": "thing", "index": None, "min_matched": 0.95}, **entry)
        self.write_json("thing.json", entry)
        os.makedirs(self.path("fetched"), exist_ok=True)
        for filename, content in files.items():
            with open(self.path("fetched", filename), "wb") as f:
                f.write(content)
        self.write_json(os.path.join("fetched", "thing.fetch.json"), dict(
            {"resource": "thing", "variant": "toy", "version": "v1", "file": next(iter(files)),
             "file_sha256": "as-published"}, **(record or {})))
        argv = ["prepare_resource.py", "--registry", self.path("thing.json"), "--variant", "toy",
                "--genome", "toy_genome", "--fetched", self.path("fetched"),
                "--chrom-sizes", self.path("chrom.sizes"), "--outdir", self.path("prepared")]
        out, err, code = io.StringIO(), io.StringIO(), 0
        with mock.patch.object(sys, "argv", argv), redirect_stdout(out), redirect_stderr(err):
            try:
                prepare_resource.main()
            except SystemExit as e:
                code = e.code or 0
        self.stdout, self.stderr = out.getvalue(), err.getvalue()
        return code

    def provenance(self):
        with open(self.path("prepared", "provenance.json")) as f:
            return json.load(f)


class Renaming(PrepareCase):
    """dbSNP's case: sequence names rewritten to the genome's own before the gate."""

    RENAMED = {"kind": "vcf", "index": "tabix", "variants": {"toy": {
        "url": "https://example.org/x", "version": "v1",
        "rename": {"report": "https://example.org/report.txt"}, "naming": "RefSeq, renamed"}}}
    NAMES = {"names_report": {"url": "https://example.org/report.txt", "file": "thing.names.txt"}}

    def mapping(self, sizes):
        with open(self.path("report.txt"), "w") as f:
            f.write(REPORT)
        return prepare_resource.name_mapping(
            prepare_resource.read_names_report(self.path("report.txt")), sizes)

    def test_every_other_name_of_a_sequence_maps_to_the_one_the_genome_uses(self):
        mapping, summary = self.mapping(SIZES)
        self.assertEqual(mapping, {
            "CM1.1": "1", "NC_1.1": "1", "chr1": "1", "CM2.1": "2", "NC_2.1": "2", "chr2": "2",
            "J01415.2": "MT", "NC_012920.1": "MT", "chrM": "MT"})
        self.assertEqual((summary["sequences_in_report"], summary["sequences_this_genome_has"]),
                         (4, 3))
        # The same report serves a genome named the UCSC way.
        mapping, _ = self.mapping({"chr1": 1000, "chr2": 800, "chrM": 100})
        self.assertEqual((mapping["NC_1.1"], mapping["1"], mapping["NC_012920.1"]),
                         ("chr1", "chr1", "chrM"))

    def test_sequence_of_another_length_is_not_the_same_sequence(self):
        mapping, summary = self.mapping(dict(SIZES, **{"1": 999}))
        self.assertNotIn("NC_1.1", mapping)
        self.assertEqual(mapping["NC_2.1"], "2")
        self.assertEqual(summary["same_name_different_length"],
                         {"1": {"genome": 999, "report": 1000}})

    def published(self):
        text = vcf(("NC_1.1", 10), ("NC_1.1", 500), ("NC_2.1", 7), ("NC_012920.1", 5),
                   contig_lines="##contig=<ID=NC_1.1,length=1000>\n##contig=<ID=NC_2.1>\n")
        return text.replace("\t.\tA\tC", "\trs1;NC_1.1\tA\tC")

    def test_renamed_file_passes_the_gate_and_is_indexed_under_the_genomes_names(self):
        import pysam
        code = self.run_prepare(self.RENAMED, {"published.vcf.gz": gzip.compress(self.published().encode()),
                                               "thing.names.txt": REPORT.encode()}, self.NAMES)
        self.assertEqual(code, 0, self.stderr)
        data = self.path("prepared", "thing.data.gz")
        with gzip.open(data, "rt") as f:
            stored = f.read()
        # Names changed in the header and in the first column, and nowhere else:
        # an ID that happens to contain an accession is left alone.
        self.assertEqual(stored, self.published().replace("ID=NC_1.1", "ID=1").replace(
            "ID=NC_2.1", "ID=2").replace("\nNC_1.1\t", "\n1\t").replace(
            "\nNC_2.1\t", "\n2\t").replace("\nNC_012920.1\t", "\nMT\t"))
        self.assertIn("rs1;NC_1.1", stored)
        with pysam.TabixFile(data) as index:
            self.assertEqual(sorted(index.contigs), ["1", "2", "MT"])
            self.assertEqual(len(list(index.fetch("1"))), 2)
        record = self.provenance()
        self.assertEqual(record["compatibility"]["share_matched"], 1.0)
        self.assertFalse(record["stored_as_published"])
        self.assertEqual(record["file_sha256"], "as-published")
        self.assertEqual((record["renamed"]["records_renamed"], record["renamed"]["sequences_renamed"]),
                         (4, 3))
        self.assertEqual(record["renamed"]["examples"]["NC_1.1"], "1")
        self.assertEqual(record["renamed"]["report"], "https://example.org/report.txt")

    def test_record_on_a_sequence_the_genome_lacks_keeps_its_published_name(self):
        records = [("NC_1.1", n) for n in range(1, 40)] + [("NW_9.1", 5)]
        code = self.run_prepare(self.RENAMED, {"published.vcf.gz": gzip.compress(vcf(*records).encode()),
                                               "thing.names.txt": REPORT.encode()}, self.NAMES)
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual(self.provenance()["compatibility"]["contigs_absent_from_genome"],
                         {"NW_9.1": 1})
        with gzip.open(self.path("prepared", "thing.data.gz"), "rt") as f:
            self.assertIn("\nNW_9.1\t5\t", f.read())

    def test_renaming_does_not_rescue_a_file_for_another_build(self):
        # Right accession, but a position past the end of this genome's sequence.
        code = self.run_prepare(self.RENAMED, {"published.vcf.gz": gzip.compress(
            vcf(("NC_1.1", 10), ("NC_1.1", 1001)).encode()), "thing.names.txt": REPORT.encode()},
            self.NAMES)
        self.assertEqual(code, 1)
        self.assertIn("beyond the end", self.stderr)
        self.assertFalse(os.path.exists(self.path("prepared")))

    def test_file_stored_as_fetched_is_linked_and_one_to_be_recompressed_is_copied(self):
        plain = {"kind": "bed", "variants": {"toy": {"url": "https://example.org/x", "version": "v1"}}}
        bed = gzip.compress(b"1\t0\t500\tx\n")
        self.assertEqual(self.run_prepare(plain, {"published.bed.gz": bed}), 0, self.stderr)
        self.assertTrue(os.path.samefile(self.path("fetched", "published.bed.gz"),
                                         self.path("prepared", "thing.data.gz")))
        shutil.rmtree(self.path("prepared"))
        self.assertEqual(self.run_prepare(dict(plain, index="tabix"), {"published.bed.gz": bed}), 0,
                         self.stderr)
        self.assertFalse(os.path.samefile(self.path("fetched", "published.bed.gz"),
                                          self.path("prepared", "thing.data.gz")))
        with open(self.path("fetched", "published.bed.gz"), "rb") as f:
            self.assertEqual(f.read(), bed)        # the fetch step's file is untouched


def cache_tarball(members):
    """A .tar.gz holding {member name: bytes}, in the order given."""
    import tarfile
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name, content in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    return buffer.getvalue()


def cache(directories, root="toy_species/116_toy"):
    """A VEP cache with the given {sequence: [region file names]} and an info.txt."""
    members = {f"{root}/info.txt": b"species\ttoy_species\nassembly\ttoy\n"}
    for name, regions in directories.items():
        for region in regions:
            members[f"{root}/{name}/{region}"] = b"x"
        members[f"{root}/{name}/all_vars.gz"] = b"x"
    return cache_tarball(members)


class VepCache(PrepareCase):
    CACHE = {"kind": "vep_cache", "variants": {"toy": {
        "url": "https://example.org/release-{release}/x.tar.gz", "version": "annotation_release",
        "naming": "Ensembl-style"}}}
    FITS = {"1": ["1-1000000.gz"], "2": ["1-1000000.gz"], "MT": ["1-1000000.gz"]}

    def check(self, tarball, sizes=SIZES, min_matched=0.95):
        with open(self.path("cache.tar.gz"), "wb") as f:
            f.write(tarball)
        found = prepare_resource.inspect_vep_cache(self.path("cache.tar.gz"))
        return prepare_resource.gate_vep_cache(found, sizes, min_matched)

    def test_cache_with_a_directory_for_each_of_the_genomes_sequences_passes(self):
        reasons, summary = self.check(cache(dict(self.FITS, CHR_HG1_PATCH=["1-1000000.gz"])))
        self.assertEqual(reasons, [])
        self.assertEqual((summary["share_matched"], summary["sequence_directories"]), (1.0, 4))
        # A patch the genome leaves out is expected in a cache, and is listed.
        self.assertEqual(summary["directories_absent_from_genome"], {"CHR_HG1_PATCH": 2})

    def test_share_is_of_the_genomes_bases_not_of_its_sequence_count(self):
        # No directory for MT: one sequence of three, but 100 bases of 1,900.
        reasons, summary = self.check(cache({"1": ["1-1000000.gz"], "2": ["1-1000000.gz"]}),
                                      min_matched=0.94)
        self.assertEqual(reasons, [])
        self.assertEqual(summary["genome_bases_on_sequences_in_the_cache"], 1800)
        self.assertEqual(len(self.check(cache({"1": ["1-1000000.gz"]}), min_matched=0.94)[0]), 1)

    def test_differently_named_cache_is_refused(self):
        reasons, summary = self.check(cache({"chr1": ["1-1000000.gz"], "chr2": ["1-1000000.gz"]}))
        self.assertEqual(summary["share_matched"], 0.0)
        self.assertEqual(len(reasons), 1)
        self.assertIn("only 0.0% of this genome's bases", reasons[0])
        self.assertIn("chr1", reasons[0])

    def test_region_starting_beyond_the_end_of_its_sequence_means_another_build(self):
        reasons, _ = self.check(cache(dict(self.FITS, **{"2": ["1-1000000.gz", "1000001-2000000.gz"]})))
        self.assertEqual(len(reasons), 1)
        self.assertIn("starts at 1,000,001", reasons[0])
        self.assertIn("2 is 800 long", reasons[0])

    def test_member_that_would_land_outside_the_asset_is_refused(self):
        reasons, _ = self.check(cache_tarball({"toy_species/116_toy/1/1-1000000.gz": b"x",
                                               "../outside.txt": b"x"}))
        self.assertTrue(any("outside the asset" in reason for reason in reasons), reasons)

    def test_archive_that_is_not_a_cache_is_refused(self):
        reasons, _ = self.check(cache_tarball({"README.txt": b"x"}))
        self.assertIn("no sequence directories", reasons[0])

    def test_cache_that_fits_is_handed_on_packed_with_what_it_says_of_itself(self):
        tarball = cache(self.FITS)
        code = self.run_prepare(self.CACHE, {"x_116_toy.tar.gz": tarball}, {"version": "116"})
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual(sorted(os.listdir(self.path("prepared"))),
                         ["TAG", "provenance.json", "thing.data.tar.gz"])
        with open(self.path("prepared", "thing.data.tar.gz"), "rb") as f:
            self.assertEqual(f.read(), tarball)
        with open(self.path("prepared", "TAG")) as f:
            self.assertEqual(f.read().strip(), "116")
        record = self.provenance()
        self.assertEqual(record["cache_root"], "toy_species/116_toy")
        self.assertIn("assembly\ttoy", record["file_header"])
        self.assertTrue(record["stored_as_published"])
        self.assertEqual(record["compatibility"]["files"], 7)

    def test_cache_cut_short_in_transfer_is_refused(self):
        tarball = cache({name: [f"part{n}.gz" for n in range(60)] for name in ("1", "2", "MT")})
        code = self.run_prepare(self.CACHE, {"x.tar.gz": tarball[:len(tarball) // 2]})
        self.assertEqual(code, 1)
        self.assertIn("could not be read to its end", self.stderr)
        self.assertFalse(os.path.exists(self.path("prepared")))


if __name__ == "__main__":
    unittest.main()
