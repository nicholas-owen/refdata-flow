# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of the databases tree: data that belongs to no genome, kept beside the vault.

Nothing is fetched. The fetch tests replace fetch_resource.open_url, the one
function that touches the network, with canned responses. The filing tests are
given the directory a fetch would have left and a throwaway output directory.
tests/test_nextflow_stub.py runs the two steps together inside the workflow.
"""
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.insert(0, os.path.join(REPO, "bin"))
import fetch_database  # noqa: E402
import fetch_resource  # noqa: E402
import file_database  # noqa: E402

FASTA = gzip.compress(b">P1 toy protein\nMKV\n")


class FakeResponse:
    def __init__(self, body=b"", headers=None):
        self._body, self.headers = io.BytesIO(body), headers or {}

    def read(self, size=-1):
        return self._body.read(size)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def tarball(members):
    """A .tar.gz holding {member name: bytes}."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name, content in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    return buffer.getvalue()


class DatabaseTestCase(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="refdata-flow-test-")
        self.addCleanup(directory.cleanup)
        self.directory = directory.name
        self.outdir = self.path("out")

    def path(self, *parts):
        return os.path.join(self.directory, *parts)

    def run_main(self, module, argv):
        out, err, code = io.StringIO(), io.StringIO(), 0
        with mock.patch.object(sys, "argv", argv), redirect_stdout(out), redirect_stderr(err):
            try:
                module.main()
            except SystemExit as e:
                code = e.code or 0
        self.stdout, self.stderr = out.getvalue(), err.getvalue()
        return code

    def fetched(self, name, filename, content, version="1.0", unpack=None, **more):
        """The directory fetch_database.py leaves, written by hand."""
        directory = self.path(f"fetched-{name}-{version}")
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, filename), "wb") as f:
            f.write(content)
        record = dict({"database": name, "description": f"The {name} database.",
                       "terms": "https://example.org/terms", "unpack": unpack,
                       "url": f"https://example.org/{filename}", "version": version,
                       "retrieved_utc": "2026-11-15T12:00:00+00:00", "file": filename,
                       "file_bytes": len(content),
                       "file_sha256": hashlib.sha256(content).hexdigest(),
                       "checksum_verified": "md5 matches https://example.org/md5"}, **more)
        with open(os.path.join(directory, f"{name}.fetch.json"), "w") as f:
            json.dump(record, f)
        return directory

    def file(self, name, directory):
        return self.run_main(file_database, ["file_database.py", self.outdir, "--name", name,
                                             "--fetched", directory])

    def catalogue(self):
        with open(os.path.join(self.outdir, "databases", "databases.tsv"), newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))


class Fetch(DatabaseTestCase):
    URL = "https://example.org/pub/toy.fasta.gz"
    RELEASE = "https://example.org/pub/reldate.txt"

    def fetch(self, responses, held="", **definition):
        path = self.path("toy.json")
        with open(path, "w") as f:
            json.dump(dict({"name": "toy", "description": "A toy database.",
                            "terms": "https://example.org/terms", "url": self.URL,
                            "version": {"url": self.RELEASE, "pattern": r"Release (\S+)"},
                            "unpack": None}, **definition), f)
        self.calls = []

        def fake_open(url, method="GET"):
            self.calls.append((method, url))
            if (method, url) not in responses:
                raise OSError(f"no route to {url}")
            return responses[(method, url)]

        with mock.patch.object(fetch_resource, "open_url", fake_open):
            code = self.run_main(fetch_database, [
                "fetch_database.py", "--definition", path, "--held", held,
                "--outdir", self.path("fetched")])
        return code, sorted(os.listdir(self.path("fetched")))

    def test_database_is_fetched_at_the_release_its_publisher_states(self):
        code, files = self.fetch({("GET", self.RELEASE): FakeResponse(b"Toy Release 2026_04\n"),
                                  ("GET", self.URL): FakeResponse(FASTA)})
        self.assertEqual((code, files), (0, ["toy.fasta.gz", "toy.fetch.json"]), self.stderr)
        with open(self.path("fetched", "toy.fetch.json")) as f:
            record = json.load(f)
        self.assertEqual((record["database"], record["version"], record["unpack"]),
                         ("toy", "2026_04", None))
        self.assertEqual(record["file_sha256"], hashlib.sha256(FASTA).hexdigest())
        self.assertEqual(record["terms"], "https://example.org/terms")
        self.assertNotIn("resource", record)

    def test_release_already_in_the_tree_is_not_fetched_again(self):
        code, files = self.fetch({("GET", self.RELEASE): FakeResponse(b"Toy Release 2026_04\n")},
                                 held="2026_02,2026_04")
        self.assertEqual((code, files), (0, ["SKIP"]))
        self.assertEqual(self.calls, [("GET", self.RELEASE)])
        self.assertIn("  SKIP   toy: version 2026_04", self.stdout)

    def test_definition_that_is_not_a_databases_is_refused(self):
        code = self.run_main(fetch_database, ["fetch_database.py", "--definition",
                                              os.path.join(REPO, "resources", "clinvar.json"),
                                              "--outdir", self.path("fetched")])
        self.assertEqual(code, 1)
        self.assertIn("database definition", self.stderr)


class Filing(DatabaseTestCase):
    def test_database_is_filed_by_name_and_version_with_its_record_and_a_catalogue_line(self):
        code = self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA, version="2026_04"))
        self.assertEqual(code, 0, self.stderr)
        home = os.path.join(self.outdir, "databases", "toy", "2026_04")
        self.assertEqual(sorted(os.listdir(home)), ["provenance.json", "toy.fasta.gz"])
        with open(os.path.join(home, "toy.fasta.gz"), "rb") as f:
            self.assertEqual(f.read(), FASTA)
        with open(os.path.join(home, "provenance.json")) as f:
            self.assertFalse(json.load(f)["unpacked"])
        # Nothing half-made is left beside it.
        self.assertEqual(os.listdir(os.path.join(self.outdir, "databases", "toy")), ["2026_04"])
        self.assertEqual(self.catalogue(), [{
            "database": "toy", "version": "2026_04", "path": "databases/toy/2026_04/",
            "file": "toy.fasta.gz", "file_bytes": str(len(FASTA)),
            "file_sha256": hashlib.sha256(FASTA).hexdigest(), "unpacked": "no",
            "url": "https://example.org/toy.fasta.gz", "terms": "https://example.org/terms",
            "retrieved_utc": "2026-11-15T12:00:00+00:00",
            "checksum_verified": "md5 matches https://example.org/md5"}])

    def test_two_versions_sit_side_by_side_and_one_already_there_is_left_alone(self):
        self.assertEqual(self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA, "1.0")), 0)
        self.assertEqual(self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA + b"x", "2.0")), 0)
        self.assertEqual([row["version"] for row in self.catalogue()], ["1.0", "2.0"])
        # The same release offered again, with different bytes: the tree keeps its own.
        again = self.fetched("toy", "toy.fasta.gz", b"different", "1.0")
        self.assertEqual(self.file("toy", again), 0, self.stderr)
        self.assertIn("already in the tree", self.stdout)
        with open(os.path.join(self.outdir, "databases", "toy", "1.0", "toy.fasta.gz"), "rb") as f:
            self.assertEqual(f.read(), FASTA)

    def test_archive_is_unpacked_and_recorded_with_the_checksum_of_the_archive(self):
        packed = tarball({"hash.k2d": b"index", "taxo.k2d": b"taxonomy", "sub/names.txt": b"n"})
        code = self.file("idx", self.fetched("idx", "idx.tar.gz", packed, "20241228", unpack="tar"))
        self.assertEqual(code, 0, self.stderr)
        home = os.path.join(self.outdir, "databases", "idx", "20241228")
        self.assertEqual(sorted(os.listdir(home)), ["hash.k2d", "provenance.json", "sub", "taxo.k2d"])
        with open(os.path.join(home, "sub", "names.txt"), "rb") as f:
            self.assertEqual(f.read(), b"n")
        row = self.catalogue()[0]
        self.assertEqual((row["unpacked"], row["file"]), ("yes", "idx.tar.gz"))
        self.assertEqual(row["file_sha256"], hashlib.sha256(packed).hexdigest())
        with open(os.path.join(home, "provenance.json")) as f:
            self.assertEqual(json.load(f)["files_unpacked"], 3)

    def test_archive_that_reaches_outside_its_directory_is_refused_whole(self):
        packed = tarball({"hash.k2d": b"index", "../escaped.txt": b"x"})
        code = self.file("idx", self.fetched("idx", "idx.tar.gz", packed, "1", unpack="tar"))
        self.assertEqual(code, 1)
        self.assertIn("was NOT filed", self.stderr)
        self.assertIn("outside the database's directory", self.stderr)
        self.assertEqual(os.listdir(os.path.join(self.outdir, "databases", "idx")), [])
        self.assertFalse(os.path.exists(os.path.join(self.outdir, "databases", "escaped.txt")))

    def test_directory_an_interrupted_run_left_without_a_record_is_replaced(self):
        home = os.path.join(self.outdir, "databases", "toy", "1.0")
        os.makedirs(home)
        with open(os.path.join(home, "toy.fasta.gz"), "wb") as f:
            f.write(b"half a fi")
        self.assertEqual(self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA, "1.0")), 0)
        with open(os.path.join(home, "toy.fasta.gz"), "rb") as f:
            self.assertEqual(f.read(), FASTA)

    def test_catalogue_is_rebuilt_from_the_tree_not_appended_to(self):
        import shutil
        self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA, "1.0"))
        self.file("other", self.fetched("other", "other.fasta.gz", FASTA, "7"))
        self.assertEqual([row["database"] for row in self.catalogue()], ["other", "toy"])
        shutil.rmtree(os.path.join(self.outdir, "databases", "other"))
        code = self.run_main(file_database, ["file_database.py", self.outdir, "--catalogue-only"])
        self.assertEqual(code, 0, self.stderr)
        self.assertEqual([row["database"] for row in self.catalogue()], ["toy"])

    def test_crate_describes_each_database_though_the_vault_holds_no_genome(self):
        packed = tarball({"hash.k2d": b"index"})
        self.file("toy", self.fetched("toy", "toy.fasta.gz", FASTA, "2026_04"))
        self.file("idx", self.fetched("idx", "idx.tar.gz", packed, "20241228", unpack="tar"))
        with open(os.path.join(self.outdir, "ro-crate-metadata.json")) as f:
            graph = {entity["@id"]: entity for entity in json.load(f)["@graph"]}
        root = graph["./"]
        self.assertEqual(sorted(part["@id"] for part in root["hasPart"]),
                         ["databases/idx/20241228/", "databases/toy/2026_04/"])
        self.assertEqual(root["datePublished"], "2026-11-15T12:00:00+00:00")

        # Kept as one file: listed, with its size and the checksum it was fetched at.
        toy = graph["databases/toy/2026_04/"]
        self.assertEqual((toy["name"], toy["version"]), ("toy:2026_04", "2026_04"))
        self.assertEqual(toy["license"], {"@id": "https://example.org/terms"})
        data = graph["databases/toy/2026_04/toy.fasta.gz"]
        self.assertEqual(data["sha256"], hashlib.sha256(FASTA).hexdigest())
        self.assertEqual(data["contentSize"], str(len(FASTA)))
        action = graph["#fetch-database-toy-2026_04"]
        self.assertEqual(action["@type"], "CreateAction")
        self.assertEqual(action["object"], {"@id": "https://example.org/toy.fasta.gz"})
        self.assertEqual(action["result"], [{"@id": "databases/toy/2026_04/toy.fasta.gz"},
                                            {"@id": "databases/toy/2026_04/provenance.json"}])
        self.assertEqual(action["endTime"], "2026-11-15T12:00:00+00:00")

        # Unpacked: described as its directory, not file by file.
        idx = graph["databases/idx/20241228/"]
        self.assertEqual(idx["hasPart"], {"@id": "databases/idx/20241228/provenance.json"})
        self.assertNotIn("databases/idx/20241228/hash.k2d", graph)
        self.assertEqual(graph["#fetch-database-idx-20241228"]["result"][0],
                         {"@id": "databases/idx/20241228/"})
        for identifier, entity in graph.items():
            if entity["@type"] == "File":
                self.assertTrue(os.path.isfile(os.path.join(self.outdir, identifier)), identifier)


if __name__ == "__main__":
    unittest.main()
