# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests of bin/vault_crate.py: the RO-Crate that describes the vault.

The vault here is written by hand - a config and a few small files in the layout
refgenie keeps - so these need only PyYAML. tests/test_ingest_e2e.py and
tests/test_nextflow_stub.py check that real runs leave a crate behind and keep it
current; whether the crate satisfies the RO-Crate specification is a question for
the reference validator, which is not a dependency of this pipeline (see the
README for the command).
"""
import hashlib
import json
import logging
import os
import sys
import tempfile
import unittest

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import vault_crate  # noqa: E402

LOGGER = logging.getLogger("test_vault_crate")
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False

MOUSE, YEAST = "a" * 48, "b" * 48
FASTA = b">1\nACGTACGTAC\n"
GTF = b'1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n'
SHA_FASTA = hashlib.sha256(FASTA).hexdigest()


class HandBuiltVault(unittest.TestCase):
    """A vault on disk, built one asset at a time, and the crate made from it."""

    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="refdata-flow-test-")
        self.addCleanup(directory.cleanup)
        self.outdir = directory.name
        self.config_file = os.path.join(self.outdir, "refgenie_export.yaml")
        self.genomes = {}
        self.save()

    def save(self):
        with open(self.config_file, "w") as f:
            yaml.safe_dump({"config_version": 0.4, "genome_folder": self.outdir,
                            "genomes": self.genomes}, f)

    def add(self, digest, asset, tag, files, parents=(), aliases=None, on_disk=True):
        """Add <digest>/<asset>:<tag> holding `files`: {seek key: (filename, content)}."""
        genome = self.genomes.setdefault(digest, {"aliases": [], "assets": {}})
        if aliases:
            genome["aliases"] = list(aliases)
        entry = genome["assets"].setdefault(
            asset, {"asset_description": f"Toy {asset}", "default_tag": tag, "tags": {}})
        entry["tags"][tag] = {
            "seek_keys": dict({key: name for key, (name, _) in files.items()}, dir="."),
            "asset_parents": list(parents),
        }
        if on_disk:
            directory = os.path.join(self.outdir, "data", digest, asset, str(tag))
            os.makedirs(directory, exist_ok=True)
            for name, content in files.values():
                with open(os.path.join(directory, name), "wb") as f:
                    f.write(content if isinstance(content, bytes)
                            else json.dumps(content).encode())
        self.save()

    def genome_with_download(self, digest=MOUSE, release="116", record=None, **extra):
        """A genome, its annotation at one release, and the record of the download."""
        self.add(digest, "fasta", "default",
                 {"fasta": (f"{digest}.fa", FASTA),
                  "fai": (f"{digest}.fa.fai", b"1\t10\t3\t10\t11\n"),
                  "chrom_sizes": (f"{digest}.chrom.sizes", b"1\t10\n"),
                  "dict": (f"{digest}.dict", b"@HD\tVN:1.0\n")},
                 aliases=["toyAsm_ensembl", "mouse"])
        self.add(digest, "ensembl_gtf", release,
                 {"ensembl_gtf": (f"{digest}.gtf.gz", GTF)},
                 parents=[f"{digest}/fasta:default"])
        if record is None:
            record = {
                "schema": 2, "provider": "Ensembl", "species": "Mus_musculus",
                "assembly": "toyAsm", "annotation_release": release,
                "retrieved_utc": "2026-11-09T10:00:00+00:00",
                "genome_url": "https://ftp.example.org/toyAsm.dna.fa.gz",
                "annotation_url": f"https://ftp.example.org/toyAsm.{release}.gtf.gz",
                "tool_versions": {"genomepy": "0.16.4", "python": "3.14.0"},
                "genome": {"file_sha256": "f" * 64, "content_sha256": SHA_FASTA},
                "annotation": {"file_sha256": "9" * 64, "content_sha256": "8" * 64},
            }
        record.update(extra)
        self.add(digest, "provenance", release,
                 {"provenance": (f"{digest}.provenance.json", record)})

    def transcriptome(self, digest=MOUSE, release="116", tool_version="0.12.7"):
        self.add(digest, "fasta_txome", release,
                 {"fasta_txome": (f"{digest}.fa", b">t1\nACGT\n"),
                  "provenance": (f"{digest}.provenance.json",
                                 {"asset": "fasta_txome", "genome": "toyAsm_ensembl",
                                  "annotation": f"ensembl_gtf:{release}", "tool": "gffread",
                                  "tool_version": tool_version})},
                 parents=[f"{digest}/fasta:default", f"{digest}/ensembl_gtf:{release}"])

    def resource(self, digest=MOUSE, stored_as_published=True):
        self.add(digest, "clinvar", "20261109",
                 {"clinvar": (f"{digest}.vcf.gz", b"vcf"),
                  "tbi": (f"{digest}.vcf.gz.tbi", b"tbi"),
                  "provenance": (f"{digest}.provenance.json", {
                      "resource": "clinvar", "variant": "GRCh38", "version": "20261109",
                      "description": "Variants and their clinical significance.",
                      "url": "https://ftp.example.org/clinvar.vcf.gz",
                      "terms": "https://www.example.org/terms",
                      "retrieved_utc": "2026-11-10T12:00:00+00:00",
                      "file_sha256": "c" * 64,
                      "stored_as_published": stored_as_published})},
                 parents=[f"{digest}/fasta:default"])

    # -- reading the crate -------------------------------------------------------

    def crate(self):
        return vault_crate.build_crate(self.outdir, self.config_file)

    def entities(self):
        return {entity["@id"]: entity for entity in self.crate()["@graph"]}

    def actions(self, prefix):
        return [e for e in self.crate()["@graph"]
                if e["@type"] == "CreateAction" and e["@id"].startswith(prefix)]

    @staticmethod
    def ids(value):
        return [item["@id"] for item in (value if isinstance(value, list) else [value])]


class Structure(HandBuiltVault):
    def test_descriptor_and_root_carry_what_the_specification_requires(self):
        self.genome_with_download()
        crate = self.crate()
        # The second context is the one that defines `sha256`.
        self.assertEqual(crate["@context"], ["https://w3id.org/ro/crate/1.1/context",
                                             "https://w3id.org/ro/terms/workflow-run/context"])
        descriptor, root = crate["@graph"][0], crate["@graph"][1]
        self.assertEqual(descriptor["@id"], "ro-crate-metadata.json")
        self.assertEqual(descriptor["about"], {"@id": "./"})
        self.assertEqual(descriptor["conformsTo"], {"@id": "https://w3id.org/ro/crate/1.1"})
        self.assertEqual((root["@id"], root["@type"]), ("./", "Dataset"))
        for required in ("name", "description", "datePublished", "license"):
            self.assertTrue(root[required], required)
        self.assertEqual(root["conformsTo"], {"@id": "https://w3id.org/ro/wfrun/process/0.5"})

    def test_every_reference_leads_to_an_entity_in_the_crate(self):
        self.genome_with_download()
        self.transcriptome()
        self.resource()
        entities = self.entities()

        def references(value):
            if isinstance(value, dict):
                if set(value) == {"@id"}:
                    yield value["@id"]
                else:
                    for item in value.values():
                        yield from references(item)
            elif isinstance(value, list):
                for item in value:
                    yield from references(item)

        described_elsewhere = {"https://w3id.org/ro/crate/1.1",
                               "http://schema.org/CompletedActionStatus"}
        for entity in entities.values():
            for target in references({k: v for k, v in entity.items() if k != "@id"}):
                if target not in described_elsewhere:
                    self.assertIn(target, entities, f"{entity['@id']} -> {target}")

    def test_vault_genome_asset_and_file_nest_by_haspart(self):
        self.genome_with_download()
        entities = self.entities()
        self.assertEqual(self.ids(entities["./"]["hasPart"]), [f"data/{MOUSE}/"])
        genome = entities[f"data/{MOUSE}/"]
        self.assertEqual((genome["identifier"], genome["alternateName"]),
                         (MOUSE, ["mouse", "toyAsm_ensembl"]))
        self.assertEqual(self.ids(genome["hasPart"]),
                         [f"data/{MOUSE}/ensembl_gtf/116/", f"data/{MOUSE}/fasta/default/",
                          f"data/{MOUSE}/provenance/116/"])
        annotation = entities[f"data/{MOUSE}/ensembl_gtf/116/"]
        self.assertEqual((annotation["name"], annotation["version"]), ("ensembl_gtf:116", "116"))
        self.assertEqual(self.ids(annotation["hasPart"]),
                         [f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"])
        self.assertEqual(self.ids(annotation["isBasedOn"]), [f"data/{MOUSE}/fasta/default/"])

    def test_lone_value_is_written_as_the_value_not_a_list_of_one(self):
        self.genome_with_download()
        entities = self.entities()
        self.assertEqual(entities["./"]["hasPart"], {"@id": f"data/{MOUSE}/"})
        self.assertEqual(entities[f"data/{MOUSE}/ensembl_gtf/116/"]["isBasedOn"],
                         {"@id": f"data/{MOUSE}/fasta/default/"})
        self.assertIsInstance(entities[f"data/{MOUSE}/fasta/default/"]["hasPart"], list)

    def test_file_has_its_size_and_format(self):
        self.genome_with_download()
        entities = self.entities()
        fasta = entities[f"data/{MOUSE}/fasta/default/{MOUSE}.fa"]
        self.assertEqual(fasta["@type"], "File")
        self.assertEqual(fasta["contentSize"], str(len(FASTA)))
        self.assertEqual(fasta["encodingFormat"],
                         ["text/plain", {"@id": "http://edamontology.org/format_1929"}])
        gtf = entities[f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"]
        self.assertEqual(gtf["encodingFormat"],
                         ["application/gzip", {"@id": "http://edamontology.org/format_2306"}])
        self.assertEqual(entities["http://edamontology.org/format_2306"]["name"], "GTF")
        # refgenie's `dir` seek key is the asset's directory, not a file in it.
        self.assertFalse([i for i in entities if i.endswith("/.")])

    def test_only_what_is_on_disk_is_listed(self):
        self.genome_with_download()
        # In the config, as refgenie leaves a skipped build, but never written.
        self.add(MOUSE, "fasta_txome", "116", {"fasta_txome": (f"{MOUSE}.fa", b"")},
                 on_disk=False)
        os.remove(os.path.join(self.outdir, "data", MOUSE, "fasta", "default", f"{MOUSE}.dict"))
        entities = self.entities()
        self.assertNotIn(f"data/{MOUSE}/fasta_txome/116/", entities)
        self.assertNotIn(f"data/{MOUSE}/fasta/default/{MOUSE}.dict", entities)
        self.assertIn(f"data/{MOUSE}/fasta/default/{MOUSE}.fa", entities)

    def test_empty_vault_is_still_a_crate(self):
        crate = self.crate()
        root = crate["@graph"][1]
        self.assertEqual(root["hasPart"], [])
        self.assertNotIn("mentions", root)
        self.assertRegex(root["datePublished"], r"^\d{4}-\d{2}-\d{2}$")

    def test_unchanged_vault_gives_the_same_crate_byte_for_byte(self):
        self.genome_with_download()
        self.genome_with_download(YEAST, "110")
        self.transcriptome()
        self.resource()
        vault_crate.write_vault_crate(self.outdir, self.config_file, LOGGER)
        path = os.path.join(self.outdir, "ro-crate-metadata.json")
        with open(path, "rb") as f:
            first = f.read()
        os.utime(self.config_file, (1, 1))
        vault_crate.write_vault_crate(self.outdir, self.config_file, LOGGER)
        with open(path, "rb") as f:
            self.assertEqual(f.read(), first)
        self.assertEqual(json.loads(first), self.crate())
        self.assertEqual([n for n in os.listdir(self.outdir) if n.startswith(".ro-crate-")], [])

    def test_date_published_is_the_latest_retrieval(self):
        self.genome_with_download()
        self.resource()
        self.assertEqual(self.crate()["@graph"][1]["datePublished"],
                         "2026-11-10T12:00:00+00:00")


class Checksums(HandBuiltVault):
    """A checksum is given only where a recorded one applies to the stored bytes."""

    def test_genome_fasta_carries_the_hash_of_the_downloads_content(self):
        self.genome_with_download()
        fasta = self.entities()[f"data/{MOUSE}/fasta/default/{MOUSE}.fa"]
        self.assertEqual(fasta["sha256"], SHA_FASTA)
        with open(os.path.join(self.outdir, fasta["@id"]), "rb") as f:
            self.assertEqual(hashlib.sha256(f.read()).hexdigest(), fasta["sha256"])

    def test_annotation_carries_the_hash_of_the_delivered_file(self):
        self.genome_with_download()
        gtf = self.entities()[f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"]
        self.assertEqual(gtf["sha256"], "9" * 64)

    def test_files_nothing_recorded_a_hash_for_have_none(self):
        self.genome_with_download()
        self.transcriptome()
        entities = self.entities()
        for name in (f"data/{MOUSE}/fasta/default/{MOUSE}.fa.fai",
                     f"data/{MOUSE}/fasta_txome/116/{MOUSE}.fa",
                     f"data/{MOUSE}/provenance/116/{MOUSE}.provenance.json"):
            self.assertNotIn("sha256", entities[name], name)

    def test_resource_stored_as_published_carries_the_publishers_hash(self):
        self.genome_with_download()
        self.resource(stored_as_published=True)
        self.assertEqual(
            self.entities()[f"data/{MOUSE}/clinvar/20261109/{MOUSE}.vcf.gz"]["sha256"], "c" * 64)

    def test_resource_recompressed_for_indexing_does_not(self):
        # The recorded hash is of the publisher's file; the vault's is a different
        # compression of the same content, so the hash would not match it.
        self.genome_with_download()
        self.resource(stored_as_published=False)
        self.assertNotIn(
            "sha256", self.entities()[f"data/{MOUSE}/clinvar/20261109/{MOUSE}.vcf.gz"])

    def test_record_without_hashes_adds_none(self):
        self.genome_with_download(record={"assembly": "toyAsm", "annotation_release": "116"})
        entities = self.entities()
        self.assertNotIn("sha256", entities[f"data/{MOUSE}/fasta/default/{MOUSE}.fa"])
        self.assertNotIn("sha256", entities[f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"])


class Actions(HandBuiltVault):
    def test_download_names_the_tool_its_sources_and_what_it_produced(self):
        self.genome_with_download()
        entities = self.entities()
        action, = self.actions("#download-")
        self.assertEqual(action["endTime"], "2026-11-09T10:00:00+00:00")
        self.assertEqual(action["actionStatus"], {"@id": "http://schema.org/CompletedActionStatus"})
        tool = entities[action["instrument"]["@id"]]
        self.assertEqual((tool["@type"], tool["name"], tool["version"]),
                         ("SoftwareApplication", "genomepy", "0.16.4"))
        self.assertTrue(tool["url"].startswith("https://"))
        self.assertEqual(self.ids(action["object"]),
                         ["https://ftp.example.org/toyAsm.dna.fa.gz",
                          "https://ftp.example.org/toyAsm.116.gtf.gz"])
        self.assertEqual(self.ids(action["result"]),
                         [f"data/{MOUSE}/fasta/default/{MOUSE}.fa",
                          f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz",
                          f"data/{MOUSE}/provenance/116/{MOUSE}.provenance.json"])
        self.assertIn(action["@id"], self.ids(entities["./"]["mentions"]))

    def test_each_release_of_a_genome_is_its_own_download(self):
        self.genome_with_download(release="116")
        self.genome_with_download(release="110")
        actions = self.actions("#download-")
        self.assertEqual(len(actions), 2)
        self.assertEqual(sorted(self.ids(a["result"])[1].split("/")[3] for a in actions),
                         ["110", "116"])

    def test_record_from_before_providers_were_named_still_finds_its_annotation(self):
        self.genome_with_download(record={"assembly": "toyAsm", "annotation_release": "116"})
        action, = self.actions("#download-")
        self.assertIn(f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz", self.ids(action["result"]))
        self.assertNotIn("endTime", action)
        self.assertNotIn("object", action)
        self.assertNotIn("version", self.entities()[action["instrument"]["@id"]])

    def test_source_that_is_not_an_address_is_described_not_linked(self):
        # genomepy builds UCSC annotations from a database table, and says so in
        # place of an address.
        self.genome_with_download(annotation_url="ucsc table: ncbiRefSeq")
        action, = self.actions("#download-")
        source = self.entities()[self.ids(action["object"])[1]]
        self.assertTrue(source["@id"].startswith("#source-"))
        self.assertEqual(source["name"], "ucsc table: ncbiRefSeq")
        self.assertNotIn("url", source)

    def test_sequence_only_download_has_no_annotation_in_it(self):
        self.genome_with_download(annotation_release="")
        action, = self.actions("#download-")
        self.assertEqual(self.ids(action["object"]), ["https://ftp.example.org/toyAsm.dna.fa.gz"])
        self.assertNotIn(f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz", self.ids(action["result"]))

    def test_transcriptome_is_derived_from_the_genome_and_its_annotation_release(self):
        self.genome_with_download()
        self.transcriptome()
        entities = self.entities()
        action, = self.actions("#derive-")
        tool = entities[action["instrument"]["@id"]]
        self.assertEqual((tool["name"], tool["version"]), ("gffread", "0.12.7"))
        self.assertEqual(self.ids(action["object"]),
                         [f"data/{MOUSE}/fasta/default/{MOUSE}.fa",
                          f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"])
        self.assertEqual(self.ids(action["result"]),
                         [f"data/{MOUSE}/fasta_txome/116/{MOUSE}.fa",
                          f"data/{MOUSE}/fasta_txome/116/{MOUSE}.provenance.json"])
        self.assertEqual(self.ids(entities[f"data/{MOUSE}/fasta_txome/116/"]["isBasedOn"]),
                         [f"data/{MOUSE}/fasta/default/", f"data/{MOUSE}/ensembl_gtf/116/"])

    def derived(self, name, tag, files, parents, **fields):
        """An asset derived from others, with the record a derive script leaves."""
        record = {"asset": name, "genome": "toyAsm_ensembl", "tool": "sometool",
                  "tool_version": "1.0", "derived_utc": "2026-11-14T09:00:00+00:00"}
        record.update(fields)
        files = dict(files, provenance=(f"{MOUSE}.provenance.json", record))
        self.add(MOUSE, name, tag, files, parents=[f"{MOUSE}/{p}" for p in parents])

    def test_any_asset_whose_record_names_itself_and_a_tool_is_a_derivation(self):
        # Nothing in the builder knows the name "tx2gene".
        self.genome_with_download()
        self.derived("tx2gene", "116", {"tx2gene": (f"{MOUSE}.tx2gene.tsv", b"t1\tg1\tG1\n")},
                     ["ensembl_gtf:116"], tool="refdata-flow", tool_version="1.2.3")
        entities = self.entities()
        action, = self.actions("#derive-tx2gene-")
        self.assertEqual(action["name"], "Derive tx2gene from ensembl_gtf:116")
        self.assertEqual(action["endTime"], "2026-11-14T09:00:00+00:00")
        tool = entities[action["instrument"]["@id"]]
        self.assertEqual((tool["name"], tool["version"]), ("refdata-flow", "1.2.3"))
        # Read from its one parent, the annotation, and not from the genome.
        self.assertEqual(self.ids(action["object"]),
                         [f"data/{MOUSE}/ensembl_gtf/116/{MOUSE}.gtf.gz"])
        self.assertEqual(self.ids(action["result"]),
                         [f"data/{MOUSE}/tx2gene/116/{MOUSE}.tx2gene.tsv",
                          f"data/{MOUSE}/tx2gene/116/{MOUSE}.provenance.json"])
        table = entities[f"data/{MOUSE}/tx2gene/116/{MOUSE}.tx2gene.tsv"]
        self.assertEqual(table["encodingFormat"], ["text/tab-separated-values",
                                                   {"@id": "http://edamontology.org/format_3475"}])

    def test_derivation_from_the_genome_alone_reads_the_fasta(self):
        self.genome_with_download()
        self.derived("twobit", "default", {"twobit": (f"{MOUSE}.2bit", b"\x43\x27\x41\x1a")},
                     ["fasta:default"], tool="faToTwoBit", tool_version="502")
        entities = self.entities()
        action, = self.actions("#derive-twobit-")
        self.assertEqual(self.ids(action["object"]), [f"data/{MOUSE}/fasta/default/{MOUSE}.fa"])
        self.assertEqual(entities[f"data/{MOUSE}/twobit/default/{MOUSE}.2bit"]["encodingFormat"],
                         {"@id": "http://edamontology.org/format_3009"})
        self.assertTrue(entities[action["instrument"]["@id"]]["url"].startswith("https://"))

    def test_an_index_made_with_the_file_is_a_result_too(self):
        self.genome_with_download()
        self.derived("gtf_indexed", "116",
                     {"gtf_indexed": (f"{MOUSE}.sorted.gtf.gz", GTF),
                      "tbi": (f"{MOUSE}.sorted.gtf.gz.tbi", b"tbi")},
                     ["ensembl_gtf:116"], tool="htslib", tool_version="1.25")
        action, = self.actions("#derive-gtf_indexed-")
        self.assertEqual(self.ids(action["result"]),
                         [f"data/{MOUSE}/gtf_indexed/116/{MOUSE}.sorted.gtf.gz",
                          f"data/{MOUSE}/gtf_indexed/116/{MOUSE}.sorted.gtf.gz.tbi",
                          f"data/{MOUSE}/gtf_indexed/116/{MOUSE}.provenance.json"])
        index = self.entities()[f"data/{MOUSE}/gtf_indexed/116/{MOUSE}.sorted.gtf.gz.tbi"]
        self.assertEqual(index["encodingFormat"], {"@id": "http://edamontology.org/format_3616"})

    def test_latest_derivation_counts_towards_date_published(self):
        self.genome_with_download()                       # retrieved 2026-11-09
        self.derived("tx2gene", "116", {"tx2gene": (f"{MOUSE}.tx2gene.tsv", b"")},
                     ["ensembl_gtf:116"])                 # derived 2026-11-14
        self.assertEqual(self.crate()["@graph"][1]["datePublished"],
                         "2026-11-14T09:00:00+00:00")

    def test_record_that_names_no_tool_or_another_asset_is_not_a_derivation(self):
        self.genome_with_download()
        self.derived("tx2gene", "116", {"tx2gene": (f"{MOUSE}.tx2gene.tsv", b"x")},
                     ["ensembl_gtf:116"], tool="")
        self.derived("bed12", "116", {"bed12": (f"{MOUSE}.bed", b"x")},
                     ["ensembl_gtf:116"], asset="something_else")
        self.assertEqual(self.actions("#derive-"), [])
        self.assertIn(f"data/{MOUSE}/bed12/116/{MOUSE}.bed", self.entities())

    def test_resource_fetch_links_its_source_and_the_publishers_terms(self):
        self.genome_with_download()
        self.resource()
        entities = self.entities()
        action, = self.actions("#fetch-")
        self.assertEqual(action["endTime"], "2026-11-10T12:00:00+00:00")
        self.assertEqual(self.ids(action["object"]),
                         ["https://ftp.example.org/clinvar.vcf.gz",
                          f"data/{MOUSE}/fasta/default/{MOUSE}.chrom.sizes"])
        self.assertEqual(self.ids(action["result"]),
                         [f"data/{MOUSE}/clinvar/20261109/{MOUSE}.vcf.gz",
                          f"data/{MOUSE}/clinvar/20261109/{MOUSE}.vcf.gz.tbi",
                          f"data/{MOUSE}/clinvar/20261109/{MOUSE}.provenance.json"])
        dataset = entities[f"data/{MOUSE}/clinvar/20261109/"]
        self.assertEqual(dataset["license"], {"@id": "https://www.example.org/terms"})
        self.assertEqual(entities[action["instrument"]["@id"]]["name"], "refdata-flow")

    def test_provenance_files_are_listed_as_describing_their_asset(self):
        self.genome_with_download()
        self.transcriptome()
        entities = self.entities()
        self.assertEqual(
            entities[f"data/{MOUSE}/fasta_txome/116/{MOUSE}.provenance.json"]["about"],
            {"@id": f"data/{MOUSE}/fasta_txome/116/"})
        self.assertEqual(entities[f"data/{MOUSE}/provenance/116/"]["about"],
                         {"@id": f"data/{MOUSE}/"})

    def test_unreadable_record_costs_the_action_not_the_crate(self):
        self.genome_with_download()
        with open(os.path.join(self.outdir, "data", MOUSE, "provenance", "116",
                               f"{MOUSE}.provenance.json"), "w") as f:
            f.write("{ not json")
        self.assertEqual(self.actions("#"), [])
        self.assertIn(f"data/{MOUSE}/fasta/default/{MOUSE}.fa", self.entities())


if __name__ == "__main__":
    unittest.main()
