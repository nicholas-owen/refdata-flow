<p align="center">
  <img src="refdata-flow_logo.svg" alt="refdata-flow logo" width="600">
</p>


This project provides an automated workflow for downloading bioinformatics reference datasets - genomic DNA sequence (FASTA) and gene annotation (GTF) - for multiple species, cataloguing them in a [refgenie](http://refgenie.databio.org/) vault, and deriving from the two the files analyses ask for beside a genome. It also fetches the third-party resources you name, such as ClinVar, the GATK known-sites files, dbSNP or a VEP cache, and registers each only after checking it against its genome; and it files databases that belong to no genome, such as Swiss-Prot or a Kraken 2 index, beside the vault.

It runs in **four stages**, wrapped by a single script:
1. **Resolve** - an interactive Python step (`resolve.py`) turns generic queries (like `grch38` or `mouse`) into exact assembly requests.
2. **Download** - a **Nextflow** pipeline asks the vault what it already holds, then fetches the rest in parallel via the **genomepy** engine.
3. **Catalogue** - once every download has finished, the same Nextflow run registers the results in a refgenie vault with human-readable aliases (the `REGISTER` process, which runs `update_refgenie.py`).
4. **Derive and check** - for each genome with an annotation, a transcriptome is cut from the vault's own FASTA and GTF with **gffread**, and the small derived files (`tx2gene`, an indexed GTF, BED12, refFlat, `.2bit`) are made beside it (see [Transcriptome](#transcriptome)). Requested resources are fetched, checked against their genome and registered (see [Third-party resources](#third-party-resources-clinvar-the-encode-blacklist)).

![The pipeline as a metro map](refdata-flow_map.svg)

## Prerequisites

- **Python >= 3.10** with `pip` and the `venv` module. The wrapper self-installs the exact, hash-locked dependencies from `requirements.txt` into a hidden virtual environment; `refgenie` 0.13.0 requires Python 3.10+.
  On Debian/Ubuntu the `venv` module's bootstrap (`ensurepip`) ships in a **separate package**: without it `python3 -m venv` creates a directory with an interpreter but no `pip`, and the wrapper fails with a confusing `.venv/bin/pip: No such file or directory`:
  ```bash
  sudo apt install python3-venv          # or python3.12-venv, matching `python3 --version`
  ```
- **Nextflow >= 25.04**, enforced by `manifest.nextflowVersion` in `nextflow.config` so an incompatible version fails immediately with a clear message rather than a parser error. Tested against 26.04.x; CI lints with 26.04.4 under the strict syntax parser. Nextflow 25.04 and later require **Java 17+**; earlier releases run on Java 11 but are not supported here.
- **Conda (required).** `run_pipeline.sh` runs the download stage under `-profile conda`, so `conda` must be on the `PATH` of the shell that launches the pipeline. See below for which distribution to use.

### Installing Conda

**Miniforge is the recommended distribution**, for three reasons specific to this pipeline:

- It defaults to **conda-forge**, which is what Nextflow asks for: it invokes `conda create ... -c conda-forge -c bioconda`. Miniconda defaults to Anaconda's `defaults` channel instead, and mixing that with conda-forge is a common source of solver failures.
- **Licensing.** Anaconda's terms restrict commercial and larger-institution use of the `defaults` channel. Miniforge never uses it.
- **It ships `mamba`.** The pipeline's conda environment is an eight-package bioconda spec, and solving that on the classic solver can take many minutes.

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
$HOME/miniforge3/bin/conda init bash
exec bash -l
command -v conda && conda --version    # must print a path and a version
```

That last check matters: Nextflow resolves `conda` from the launching shell's `PATH`, so if `conda init` hasn't taken effect the pipeline fails on its first task.

With Miniforge installed, `mamba` is available and you can set `conda.useMamba = true` in the `conda` profile to speed up environment creation considerably. It is **off by default** because it fails outright when mamba is absent.

> **First run builds the environment.** Nextflow caches it under `work/conda/`, so later runs are fast, but `run_pipeline.sh --clean` deletes `work/` and therefore destroys that cache along with the resume cache. `conda.createTimeout` is set to 60 minutes in `nextflow.config` because the default of 20 can expire mid-solve.

**Docker / Singularity:** a `docker` profile exists but its image tag is **unverified and untested**: confirm the tag against [quay.io](https://quay.io/repository/biocontainers/genomepy?tab=tags) before relying on it. There is currently **no** Singularity/Apptainer profile. Conda is the only supported and tested execution path.

> **Annotation requires the UCSC tools.** genomepy shells out to `genePredToGtf`, `gtfToGenePred`, `genePredToBed`, `bedToGenePred` and `gff3ToGenePred` to convert gene annotations, and to `bgzip`/`samtools` for compression and indexing. **The `conda` profile declares all of these explicitly**, so on the supported path you need do nothing.
>
> They are *not* pip-installable, so they are absent from `requirements.txt` and from the wrapper's virtual environment. If you bypass the conda profile - for example running `nextflow run main.nf -profile standard` directly - sequence downloads will succeed and every row with `annotation=true` will fail. `download.py` checks **before** downloading anything, rather than failing after a multi-gigabyte FASTA has already been fetched. If `genePredToGtf` is missing it exits with code 3 naming the fix; if only the other converters are absent it warns and proceeds. That gate matches genomepy's own `check_ucsc_tools()`, which tests the same single tool, so the preflight can never reject a run genomepy would have accepted. To install them by hand:
> ```bash
> conda install -c bioconda ucsc-genepredtogtf ucsc-gtftogenepred ucsc-genepredtobed ucsc-bedtogenepred ucsc-gff3togenepred
> ```
> Requests with `annotation=false` need none of this.

Dependency versions are pinned and hash-locked in `requirements.txt` for reproducibility. To upgrade, edit `requirements.in` and regenerate with `pip-compile --generate-hashes --output-file=requirements.txt requirements.in`.

### Running the tests

The tests need no network and no extra packages. Run them with the pipeline's own interpreter, after `run_pipeline.sh` has been run once so that `.venv` exists:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

- `tests/test_resolve.py` drives the resolver with the provider search faked.
- `tests/test_update_refgenie.py` covers the vault naming, alias and default-tag rules.
- `tests/test_preflight.py` covers which requests the vault already satisfies, with the provider's "current release" supplied by the test.
- `tests/test_download.py` covers the download record, with the download itself replaced.
- `tests/test_resources.py` covers third-party resources: the registry, the fetch step with the network replaced, the check of a resource against a genome, the renaming of one that names its sequences differently, and the check of a VEP cache.
- `tests/test_databases.py` covers the databases that belong to no genome: the fetch step with the network replaced, and the filing of a database beside the vault.
- `tests/test_vault_crate.py` covers the vault's RO-Crate, built from a hand-written vault.
- `tests/test_derive.py` runs the two derive scripts with the real tools on a hand-written annotation and genome. Skipped where the tools are missing; see below.
- `tests/test_ingest_e2e.py` builds a real refgenie vault from 40-base toy genomes: one assembly from two providers, plus an older pinned release. It is skipped, with the reason printed, where refgenie or `samtools` is missing or the filesystem has no symlinks.
- `tests/test_nextflow_stub.py` executes `main.nf`. Most of it uses `nextflow run -stub`, which replaces the download with a toy assembly directory and leaves everything around it real: request validation, the duplicate guard, the publish layout. The `Registration` tests run the whole thing, download and catalogue, in one Nextflow run, with and without `--cleanup`. The retry tests run the real script with a `python3` that only fails. It is skipped where `nextflow` is not on `PATH`, and takes about four minutes.

Some tests run the real derivation tools and are skipped where those are not on `PATH`; the pipeline's Python environment does not include them. `RealDerivation` (two tests) needs `gffread`. `tests/test_derive.py` needs `gtfToGenePred`, `genePredToBed`, `bgzip` and `tabix` for the annotation assets, and `faToTwoBit` and `twoBitInfo` for the `.2bit`. To run them all, put the tools on `PATH` first, for example from a conda environment:

```bash
conda create -n refdata-tools -c bioconda gffread=0.12.7 htslib ucsc-gtftogenepred ucsc-genepredtobed ucsc-fatotwobit ucsc-twobitinfo
```

You can run the stub by hand too. Nothing is downloaded, and the `provenance.json` it writes is marked `"stub": true`:

```bash
nextflow run main.nf -stub -profile standard --input requests/requests_resolved.csv --outdir /tmp/stub-out
```

## Directory Structure

```
├── bin/
│   ├── derive_annotation.sh # Makes tx2gene, the indexed GTF, BED12 and refFlat (DERIVE_ANNOTATION)
│   ├── derive_twobit.sh     # Writes a genome as .2bit and checks it (DERIVE_TWOBIT)
│   ├── download.py          # Fetches one genome and its annotation (FETCH_GENOME)
│   ├── fetch_database.py    # Fetches one database that belongs to no genome (FETCH_DATABASE)
│   ├── fetch_resource.py    # Fetches one third-party resource (FETCH_RESOURCE)
│   ├── file_database.py     # Files a database beside the vault and rewrites its catalogue (FILE_DATABASE)
│   ├── preflight.py         # Lists the requests the vault does not already hold (PREFLIGHT)
│   ├── prepare_resource.py  # Renames a resource if it must be, checks it against its genome and indexes it (PREPARE_RESOURCE)
│   ├── register_asset.py    # Registers a transcriptome or a resource (REGISTER_ASSET)
│   ├── resource_registry.py # Reads the descriptions in resources/ and databases/
│   ├── resolve.py           # Interactive resolver: generic queries to exact assemblies
│   ├── update_refgenie.py   # Registers downloads in the refgenie vault (REGISTER)
│   └── vault_crate.py       # Describes the vault as an RO-Crate (ro-crate-metadata.json)
├── conf/
│   └── sge.config           # Scheduler profile for the UCL clusters
├── data/
│   └── references/          # Default output directory for genomes
├── databases/               # Databases that belong to no genome and can be requested, one file each
├── main.nf                  # The Nextflow pipeline
├── nextflow.config          # Nextflow configuration profiles
├── README.md
├── recipes/                 # refgenie recipes this pipeline builds assets with
├── resources/               # Third-party resources that can be requested, one file each
├── tests/                   # Offline test suite
├── requests/
│   ├── requests_raw.csv.example   # Tracked template, copied on first run
│   ├── requests_raw.csv           # The file you edit (not tracked)
│   ├── databases.csv.example      # Tracked template for the list of databases
│   ├── databases.csv              # Your list of databases, if you want any (not tracked)
│   └── requests_resolved.csv      # Resolver output and cache (not tracked)
└── run_pipeline.sh          # Unified wrapper script to run everything
```

> **Neither request file is tracked by git.** `requests_raw.csv` is your genome list, so tracking it would make every `git pull` a merge conflict against whatever the repository happens to ship. `requests_resolved.csv` is the resolver's cache. The first run copies `requests_raw.csv.example` into place and stops so you can edit it.

## How to use

### Step 1: Write Generic Requests
Edit the `requests/requests_raw.csv` file to list the generic species you want. 
Example `requests/requests_raw.csv`:
```csv
species,query,provider,annotation,aliases,release
Homo_sapiens,grch38,ensembl,true,human,116
Mus_musculus,mouse,,true,mm39;mm,
```
*(Note: If you leave `provider` empty like for Mouse, it will search all providers. You can also supply multiple custom aliases separated by semicolons in the `aliases` column.)*

> **Name the provider when you know it.** A row with a provider searches only that provider, which takes seconds. A row without one searches all four, and NCBI's assembly summaries alone take several minutes to download. It also keeps the row independent of the others' availability: if the provider you named is unreachable the row fails saying so, and an outage anywhere else is not your problem. A provider name genomepy does not know (`ensemble`) is refused before anything is searched.

> **`release` is optional, and only Ensembl can honour it.** Leave it empty to take whatever the provider currently publishes. Supplying it pins the download to that Ensembl release, and the request is checked before anything is fetched: the release must exist *and* must contain the assembly you asked for, so a typo or a release predating the assembly fails immediately with the valid alternatives listed.
>
> **It is rejected, not ignored, for other providers.** genomepy accepts a release for GENCODE, UCSC and NCBI and then silently downloads the current files anyway, which would record a version that was never used. A `release` on any provider but Ensembl therefore fails the row rather than running unpinned. UCSC cannot be pinned even in principle: it patches assemblies in place and builds annotations live from its MySQL server.
>
> The annotation version is **recorded either way**, pinned or not - see *Annotation release tagging* below.

#### A genome that no provider lists: request it by address

Some reference genomes cannot be found by searching, because they are not in any provider's catalogue: NCBI's GRCh38 "analysis sets", the Broad's `Homo_sapiens_assembly38.fasta` that the GATK resource bundle is built against, a consortium's own release. They are published as a file at an address. Give the row provider `url` and put the address in a `url` column:

```csv
species,query,provider,annotation,aliases,release,url,annotation_url
Homo_sapiens,GRCh38_no_alt,url,false,hs38_noalt,,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz,
```

Such a row is taken as written. Nothing is searched and nothing is chosen for you:

- **`query` is the genome's name.** With no provider to supply one, the name is yours. It is registered as `<name>_url`, so the example is reachable as `GRCh38_no_alt_url`, `GRCh38_no_alt` and `hs38_noalt`.
- **An annotation must be given by address too**, in `annotation_url`, as a GTF, GFF or BED file. It is registered as `url_gtf`. Leave `annotation` false to take sequence alone.
- **A release cannot be pinned.** An address is one file; a newer one is a different address.
- **The address is restricted** to `http`, `https` or `ftp` and to letters, digits and `. _ ~ / % + = : @ , -`. It is passed to a shell command, so query strings and anything a shell would act on are refused.
- **A name means one address.** If `GRCh38_no_alt` is already in the vault from one address and a later request gives another, the run stops before downloading and tells you to use a different name. Otherwise the vault would hand back the first file under the second request's name.

This matters most for variant calling. Your Ensembl human genome names its contigs `1`, `2`, `X`. The GATK known-sites files, the ENCODE blacklist and gnomAD all use `chr1`, `chr2`, `chrX`. They match the analysis-set and Broad genomes as they are, and do not match an Ensembl one.

The columns are optional. A request file without them works as before.

#### Third-party resources: ClinVar, the ENCODE blacklist

A genome and its annotation come from a provider. Some things you need alongside them come from somebody else entirely: a variant database, a list of regions to exclude. Ask for them in a `resources` column, as `name:build` pairs separated by semicolons:

```csv
species,query,provider,annotation,aliases,release,url,annotation_url,resources
Homo_sapiens,GRCh38,ensembl,TRUE,human,116,,,clinvar:GRCh38
Homo_sapiens,GRCh38_no_alt,url,false,hs38_noalt,,https://.../GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz,,blacklist:hg38
```

```bash
refgenie seek human/clinvar              # the VCF, with its .tbi beside it
refgenie seek hs38_noalt/blacklist       # the BED file
refgenie seek human/clinvar.provenance   # where it came from, and what the check found
```

| Resource | Builds | What you get |
|---|---|---|
| `clinvar` | `GRCh38`, `GRCh37` | NCBI's VCF, verified against its published md5, with a tabix index. Named `1`, `2`, `X`. About 200 MB, replaced weekly. |
| `blacklist` | `hg38`, `hg19`, `mm10`, `dm6`, `dm3`, `ce11`, `ce10` | The ENCODE blacklist v2, stored as published. Named `chr1`, `chr2`, `chrX`. A few kB. |

The GATK known-sites files, dbSNP and the VEP cache are requested the same way and have a section of their own, [below](#variant-calling-the-gatk-known-sites-files-dbsnp-a-vep-cache).

**You name the build because the pipeline cannot work it out.** Two builds of one species can share every contig name, and only you know which you hold.

**Nothing is registered until it has been checked against the genome.** A resource is a list of positions on named sequences, and whether those mean anything depends on the genome beside them. Tools do not complain when they do not match: a known-sites file that says `chr1` beside a genome that says `1` simply matches nothing, and the analysis runs to the end. So each file is compared with the genome's own sequence names and lengths first, and refused if:

- too few of its records are on sequences the genome has (different naming);
- any record lies beyond the end of its sequence (right names, wrong build);
- it declares sequence lengths that differ from the genome's.

A refusal says which, and the run exits non-zero with the vault untouched:

```text
ERROR: blacklist:hg38 was NOT registered for GRCh38_p14_ensembl:
  - only 0.0% of its 636 records are on sequences this genome has (at least 100% is required).
    The genome names its sequences 1, 10, 11, 12, ...; the resource names them chrX, chr3, chr5, chr2, ....
  UCSC: chr1, chr2, chrX. Fits a UCSC-named GRCh38 genome, such as an analysis set; not an Ensembl one.
```

That example is the real result of asking for the blacklist on an Ensembl human genome. It belongs on a UCSC-named one, which is what requesting a genome by address is for.

The check cannot prove a match. Two builds can agree on every name and length that a small file happens to touch. It catches what can be caught, and each thing it catches is certain.

Other things to know:

- **The version is the tag.** ClinVar is tagged with the date NCBI last changed the file, so weekly releases sit side by side and the newest is the default. The blacklist is `v2`.
- **It is only fetched when it has changed.** Each run asks the publisher for the current version first and skips the download if the vault has it. With `--check_releases false`, whatever is held is kept.
- **The publisher is asked on every run, including a resumed one.** `run_pipeline.sh` always resumes, and Nextflow replays a step whose inputs have not changed without running it. A download is replayed that way; the answer "still current" never is.
- **A few records may sit on sequences your genome lacks.** ClinVar has about fifty on unplaced scaffolds that a primary assembly leaves out. That passes, and they are listed in the asset's provenance.
- **Resources follow the genome, not the row.** Two rows for one genome that both ask for ClinVar fetch it once.

**Adding a resource** is two small files, with no code: a description in `resources/<name>.json` (its address for each build, how it is versioned, whether to index it) and a recipe in `recipes/<name>.json`. `resources/clinvar.json` and `resources/blacklist.json` are the patterns to copy, and `bin/resource_registry.py` lists every field. A file that names its sequences differently from any genome can be given a `rename` entry (see dbSNP below). Kinds other than VCF, BED and a VEP cache are not handled.

#### Variant calling: the GATK known-sites files, dbSNP, a VEP cache

These go in the same `resources` column. What differs is which genome each one fits.

| Resource | Builds | What you get |
|---|---|---|
| `gatk_dbsnp138`, `gatk_known_indels`, `gatk_mills` | `hg38` | The three known-sites files GATK's base recalibration is documented with, from the Broad's resource bundle. 1.6 GB, 60 MB and 20 MB. |
| `gatk_1000g_snps`, `gatk_hapmap`, `gatk_omni` | `hg38` | The training and truth sets for variant quality score recalibration, from the same bundle. 1.9 GB, 60 MB and 50 MB. |
| `dbsnp` | `GRCh38`, `GRCh37` | dbSNP build 157 from NCBI, with its sequence names rewritten to your genome's. About 30 GB. |
| `vep_cache` | `GRCh38`, `GRCh37`, `GRCm39`, `R64-1-1` | The Ensembl VEP cache for the release your annotation is at, unpacked. Human is about 28 GB packed. |

**The GATK files fit the genome they were made for.** They name sequences `chr1`, `chr2`, `chrX`, and are refused on an Ensembl genome. Request the Broad's genome by address and ask for them on that row:

```csv
species,query,provider,annotation,aliases,release,url,annotation_url,resources
Homo_sapiens,hg38_broad,url,false,hs38,,https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta,,gatk_dbsnp138:hg38;gatk_known_indels:hg38;gatk_mills:hg38
```

```bash
refgenie seek hs38/fasta        # the reference, with its .fai and .dict beside it
refgenie seek hs38/gatk_mills   # a known-sites VCF, with its .tbi beside it
```

Each is tagged `v0`, the Broad's own name for this edition of the bundle, and verified against the md5 that Google Cloud Storage sends with the download.

**dbSNP fits any GRCh38 genome, because it is renamed on arrival.** NCBI publishes it with RefSeq accessions (`NC_000001.11`), which no genome in a vault uses. So the pipeline also fetches NCBI's assembly report for the build, which lists every sequence under each of its names with its length, and rewrites each name to the one your genome uses: `1` for an Ensembl genome, `chr1` for the Broad's.

- **A sequence is renamed only where your genome's sequence of that name has the length the report gives.** The right name on another build is not a match.
- **Nothing else in a record is changed**, and the renamed file then has to pass the same check as any other resource.
- **About 6% of its records keep their RefSeq names.** They are on patches and alternate loci that a primary assembly does not have. They are counted, by sequence, in the asset's provenance file.
- **It is pinned to build 157** in dbSNP's archive and tagged `b157`. A newer build is a change to `resources/dbsnp.json`.
- **It needs room and time.** The download is about 30 GB; the renamed copy is as much again in Nextflow's `work/`, and the vault's copy a third. Renaming and indexing run on the machine you launched from and should be expected to take hours; this has not been timed on the full file.

**A VEP cache follows your annotation.** Ask for `vep_cache:GRCh38` on a row with an Ensembl annotation, and the cache fetched is the one for the release that annotation is at. It is tagged with that release, like the transcriptome.

```bash
vep --offline --cache --dir_cache "$(refgenie seek human/vep_cache)" --cache_version 116 \
    --species homo_sapiens --assembly GRCh38 ...
```

- **It is unpacked in the vault**, so the path `refgenie seek` prints is the directory VEP takes as `--dir_cache`.
- **A genome held at two releases gets a cache for each.** A row pinned to 110 and a row at 116 are two caches, tagged `110` and `116`.
- **It is only for a row with an Ensembl annotation.** On any other row it is refused before anything is fetched, saying so: there is no release to follow.
- **It is checked against the genome in the form a cache allows.** A cache is a directory per sequence, holding a file per million bases. At least 95% of your genome's bases must be on sequences the cache has a directory for, and no file may start beyond the end of its sequence.
- **Ensembl posts no md5 for it.** The archive is compared with the size Ensembl announces and read to its end before it is registered, which catches a transfer that was cut short.
- **Use it with the VEP of the same number.** A cache is tied to its release.
- **Yeast is numbered differently.** Ensembl publishes yeast through Ensembl Genomes, whose releases have their own numbers: 63 alongside Ensembl 116. An Ensembl yeast annotation is therefore tagged `63`, and so is its cache, which is fetched from Ensembl Genomes and used with `--cache_version 63`.
- **A release with no cache is not retried.** If the publisher has no file for the release your annotation is at, the fetch stops once and says so.

**The size is checked before a large download starts.** Each fetch compares the size the publisher announces with the space free in Nextflow's `work/` and stops, without retrying, if the file cannot fit. That is a floor, not a budget: a resource is held once as fetched and once in the vault, and a renamed one a third time in between, until `work/` is cleaned.

#### Databases that belong to no genome

Swiss-Prot, Pfam, SILVA and a Kraken 2 index are reference data too, but they are not written in any assembly's coordinates. So they have no row in the request file, which is a list of genomes, and no place in refgenie, where every asset hangs off a genome. List them in a file of their own:

```bash
cp requests/databases.csv.example requests/databases.csv
```

```csv
database
swissprot
pfam
```

`run_pipeline.sh` reads `requests/databases.csv` on every run if it exists. `--databases <file>` names another list.

| Database | What you get |
|---|---|
| `swissprot` | UniProtKB/Swiss-Prot as one FASTA file, tagged with the UniProt release (`2026_04`). About 90 MB. |
| `pfam` | `Pfam-A.hmm.gz`, tagged with the Pfam release (`39.0`) and verified against EBI's md5. About 660 MB. |
| `silva` | SILVA 138.2 SSU Ref NR99 with taxonomy, verified against its md5. About 200 MB. |
| `kraken2_viral`, `kraken2_standard` | A prebuilt Kraken 2 index of 2024-12-28, unpacked. 570 MB and 69 GB packed. |

They are filed beside the vault, by name and version:

```text
<outdir>/databases/
├── databases.tsv                  <-- one line per database and version
├── swissprot/2026_04/
│   ├── uniprot_sprot.fasta.gz
│   └── provenance.json
└── kraken2_viral/20241228/        <-- the directory to give kraken2 as --db
    ├── hash.k2d, opts.k2d, taxo.k2d, ...
    └── provenance.json
```

- **The version is worked out before anything is fetched**, from the publisher's own release file where there is one, and a version already in the tree is not fetched again.
- **A new release is a new directory beside the old one.** Nothing is replaced or removed, so an analysis can keep naming the version it used. There is no "current" link; `databases.tsv` lists what is held.
- **There is no genome to check it against.** It is verified against the publisher's checksum where one is posted (Pfam, SILVA). UniProt and the Kraken 2 indexes post none the pipeline can use, so those are compared with the announced size only, and the record says so.
- **A version appears in one step.** It is assembled out of sight and renamed into place, so a reader never finds half a database.
- **`databases.tsv` is rewritten from the tree on every run**, like the vault's index, and each version is listed in the vault's RO-Crate.
- **It needs the registration steps**, so it cannot be combined with `--register false`.

**Adding a database** is one file, `databases/<name>.json`: its address, how it is versioned, where its checksum is, and `"unpack": "tar"` if it is published as an archive. The five that ship are the patterns to copy.

> **Three things are checked before anything is downloaded**, and each fails the row with the line number and the reason:
>
> - **`annotation=true` on an assembly with no annotation.** The provider search already says whether one is published, so the request stops here instead of after the genome has been fetched.
> - **Two rows that resolve to the same download.** `GRCh38` and `human` both resolve to `GRCh38.p14`; that is one genome, so put it on one row with both sets of aliases. The same assembly at two *releases*, or from two *providers*, is two different downloads and is allowed.
> - **One alias on two different genomes.** An alias names exactly one genome. Sharing an alias between two releases of the same assembly is fine, because that is the same genome.

> **Aliases are labels, not conversions.** An alias is a second name for whatever you downloaded; it does not change the data. Aliasing an Ensembl assembly as `hg38` is legal and tempting, but `hg38` is UCSC's name and implies UCSC conventions: `chr1`, `chr2`, `chrM`. An Ensembl genome has contigs named `1`, `2`, `MT`, so `refgenie seek hg38/fasta` would hand a downstream tool exactly the naming it did not expect. Keep provider-specific nicknames for genomes from that provider.

### Step 2: Run the Pipeline Wrapper
Instead of memorizing the Python and Nextflow commands, simply run the provided shell script:

```bash
bash run_pipeline.sh requests/requests_raw.csv --outdir /path/to/your/custom/storage
```
*(If you omit `--outdir`, it defaults to `data/references`)*

This script will automatically:
1. **Resolve**: Parse your `requests/requests_raw.csv`. It connects to the databases and interactively prompts you if there are multiple matches for your generic queries.
2. **Download**: Automatically launch the Nextflow pipeline on the resolved targets using Conda. It first asks the vault which requests it already holds, and downloads only the rest (see [Genomes already in the vault](#-genomes-already-in-the-vault-are-not-downloaded-again)).
3. **refgenie catalogue**: Finally, in the same Nextflow run, it automatically initializes a `refgenie_export.yaml` for you (if one doesn't exist yet), computes the sequence digests, and links your aliases. The downloaded source directories are **kept** unless you pass `--cleanup` (see below).
4. **Transcriptome**: Derive a transcriptome for every genome that has an annotation and does not have one yet (see [Transcriptome](#transcriptome); `--no-txome` skips this).
5. **Small derived assets**: Derive a transcript-to-gene table, an indexed GTF, BED12 and refFlat from each annotation, and a `.2bit` from each genome, wherever the vault lacks them (see [Small derived assets](#small-derived-assets); `--no-derived` skips this).
6. **Resources and databases**: Fetch the third-party resources your rows ask for and register each against its genome once it has been checked, and file any databases listed in `requests/databases.csv` beside the vault (see [Databases that belong to no genome](#databases-that-belong-to-no-genome)).

> **Cataloguing is part of the Nextflow run.** It used to be a separate step the wrapper ran after Nextflow returned. It is now the `REGISTER` process, with three properties worth knowing:
>
> - **One task, in request order.** Genomes are registered in the order they appear in your request file, not the order their downloads finish. That order decides who gets a bare assembly name when two providers share it, so it must not vary from run to run.
> - **It always runs on the machine that launched Nextflow**, whatever executor a profile selects, using the Python environment `run_pipeline.sh` builds. It is never cached, so every run re-checks the vault; for a genome already there this takes seconds.
> - **A download failure stops before it.** If any download fails after its retries, nothing is registered in that run. Re-running picks up the completed downloads from the resume cache.
>
> To download without cataloguing, run Nextflow directly with `--register false`.

### 💡 Pro-Tip: Skipping the Wizard & Caching
The interactive wizard remembers your choices. If you run the pipeline again to add a new species, it will look at `requests/requests_resolved.csv`. If an older species has already been resolved, **it will automatically skip the slow network search and interactive prompt** for that species to save time. 

If you ever *want* to change the assembly for a species, simply delete that specific row from `requests/requests_resolved.csv` (or delete the file entirely) and the wizard will prompt you again.

### 🗄️ Reclaiming space: `--cleanup` (opt-in)

Once a genome has been ingested into the refgenie vault, its original download directory is redundant. With `--cleanup` the sources are not kept under `--outdir`, and copies that an earlier run left there are removed:
```bash
bash run_pipeline.sh requests/requests_raw.csv --cleanup
```
Two safeguards apply, because deleting downloaded reference data is not always reversible: providers archive and reorganise old releases, so re-fetching the identical files months later is not guaranteed:

- **Nothing is removed unless the ingest is verified.** The FASTA must have been registered, and if you requested annotation, the GTF must have been registered too. A partial ingest leaves the source directory in place and the run exits non-zero.
- **Provenance survives cleanup.** Before removal, `provenance.json`, genomepy's `README.txt` (source URLs) and, where the provider supplies an assembly accession, `assembly_report.txt` (sequence-name mappings between UCSC / Ensembl / GenBank accessions) are copied to `<outdir>/provenance/<provider>/<species>/<assembly>/`, which is never deleted. The remaining files are the sequence itself and indexes that refgenie regenerates in the vault. `provenance.json` and `assembly_report.txt` are also in the vault itself, as the `provenance` and `assembly_report` assets, whether or not you clean up.

Nextflow's own copy of each download stays in `work/` either way, which is what lets a later run resume. It goes when you run `bash run_pipeline.sh --clean`.

Prior to 0.9.1 this cleanup ran unconditionally, including on partially-failed ingests, and took the provenance files with it.

### 🧹 Pro-Tip: Storage Optimization
Nextflow can generate a lot of background cache in the `.nextflow/` and `work/` directories over time. To instantly purge all logs and pipeline cache to recover your disk space, simply run:
```bash
bash run_pipeline.sh --clean
```

### 🔍 Pro-Tip: Dry-Run Mode
If you want to parse your `requests/requests_raw.csv` and resolve exactly which datasets and providers will be used without actually downloading gigabytes of data, run the pipeline in dry-run mode:
```bash
bash run_pipeline.sh --dry-run
```

### 📈 Pro-Tip: Debug / Progress Mode
By default Nextflow shows only a single status line per genome while it downloads. To stream live per-genome progress messages (`[1/3] Starting`, `[2/3] Downloading…`, `[3/3] Writing provenance…`) to the console, add `--debug`:
```bash
bash run_pipeline.sh --debug
```
Output from genomes running in parallel is interleaved, with each line tagged by its assembly. For the detailed byte-level download bar, tail an individual task's log instead: the hash in the status line (e.g. `de/f44c46`) is its work directory:
```bash
tail -f work/de/f44c46*/.command.log
```

### 🧭 Genomes already in the vault are not downloaded again

Every run starts with a preflight: it reads the vault and lists each request as `SKIP` or `FETCH`, with the reason.

```text
Preflight: checking 3 request(s) against the vault at /shared/reference_datasets
  SKIP   GRCh38.p14 from Ensembl, release 116: genome and ensembl_gtf:116 are in the vault
  SKIP   sacCer3 from UCSC: genome and ucsc_gtf:20210210 are in the vault, and 20210210 is current
  FETCH  GRCm39 from Ensembl: Ensembl now publishes ensembl_gtf release 117; the vault has 116
Preflight: 1 to download, 2 already in the vault.
```

This is what makes re-running cheap, and it does not depend on Nextflow's `work/` directory. You can run `--clean`, lose `work/`, or run on an HPC that never keeps it, and a finished genome is still not fetched twice: the vault itself is the record.

A request is skipped when the vault really holds what it asks for:

- **The genome**, with its sequence dictionary on disk.
- **The annotation, at the right release**, if you asked for one. For a pinned request that is the pinned release. For an unpinned request it is the release the provider publishes today.

So an unpinned row follows the provider. When Ensembl moves from 116 to 117, the next run fetches that row again and the vault gains a `117` tag beside `116`. Pin the release if you do not want that.

Three things to know:

- **Unreachable providers never cause a refetch.** If the current release cannot be checked and the vault holds some release of the annotation, the row is skipped and the report says so.
- **To check nothing online**, run Nextflow with `--check_releases false`. Any release already in the vault then satisfies an unpinned row.
- **To download everything again regardless**, pass `--refetch`:
  ```bash
  bash run_pipeline.sh requests/requests_raw.csv --refetch
  ```

Rows that are skipped still have their aliases applied, so adding an alias to a genome you already hold takes seconds.

### 🔁 Pro-Tip: Resuming an Interrupted Run
The preflight above covers genomes that are *in the vault*. Nextflow's resume cache covers the gap before that: a download that completed but had not been registered when the run stopped. `run_pipeline.sh` always launches Nextflow with `-resume`, so simply re-running continues where the previous run left off:
```bash
bash run_pipeline.sh requests/requests_raw.csv        # automatically resumes the previous session
```
A couple of things to be aware of:

- Resume works at the **per-genome** level: it skips genomes whose download *completed successfully* in a prior run. A genome that was still downloading when you interrupted it did not complete, so it restarts from the beginning (there is no partial/byte-level resume). Re-running is safe: finished genomes are skipped, only unfinished ones re-run.
- Resume relies on the `work/` and `.nextflow/` directories. Running `bash run_pipeline.sh --clean` (or deleting those directories) **erases the resume cache**. Genomes already in the vault are unaffected, because the preflight skips them; what is lost is any download that had finished but was not yet registered.
- On a resumed run, skipped genomes are shown as `cached` in the Nextflow output, e.g. `… FETCH_GENOME (…) [100%] 1 of 1, cached: 1 ✔`.

### 🔂 Failed downloads are retried

A download that fails is tried again, twice, waiting 30 seconds and then 60. Most failures are a provider that is briefly unreachable or a dropped connection, and the request itself is fine. Two cases are handled differently:

- **A missing tool is never retried.** `download.py` exits with code 3 when the UCSC annotation tools are absent, and trying again cannot fix that.
- **When the attempts run out, other downloads finish first.** The run still exits non-zero and nothing is ingested, but genomes that were mid-download complete and are kept for the next `-resume`.

At most four downloads run at once, as a courtesy to the providers. The three settings are `max_retries`, `retry_wait` (seconds) and `max_downloads` in `nextflow.config`. `run_pipeline.sh` does not forward them, so change them there, or pass them to a bare `nextflow run` as `--max_retries 0`.

## Output & Provenance

For each requested assembly, Nextflow will generate an output directory structured by provider, species, and assembly.

Example:
`/path/to/your/custom/storage/ensembl/Homo_sapiens/GRCh38.p14/`

A request pinned to a release gets one extra level, so two releases of the same assembly can be downloaded side by side without one overwriting the other:
`/path/to/your/custom/storage/Ensembl/Homo_sapiens/release-116/GRCh38.p14/`

The `provider` and `species` directory names are used **exactly as they appear in your resolved CSV** - so a `provider` of `ensembl` produces `ensembl/`, not `Ensembl/`. Keep the casing consistent between runs, or you will end up with parallel directory trees for the same provider.

Inside each assembly directory, you will find:
- The downloaded compressed `.fa.gz` (and `.gtf.gz` if requested).
- A `README.txt` generated natively by genomepy detailing the source URLs.
- A `provenance.json` describing the download: see *The download record* below.

### The download record (`provenance.json`)

Each download writes a record of what was asked for, where it came from, and what arrived. It also goes into the vault, so it outlives the download and can be found by someone who did not build the vault:

```bash
refgenie seek human/provenance          # the record for the newest annotation release held
refgenie seek human/provenance:115      # for a particular release
```

| It records | From |
|---|---|
| Species, assembly, provider, the release you pinned (if any) | your request |
| Source URLs for the genome and the annotation, assembly accession, taxon ID | genomepy's `README.txt`, parsed into fields |
| What genomepy changed: contigs it removed, how it adjusted the annotation | genomepy's `README.txt` |
| The masking and alt-contig settings the genome was fetched with | the arguments this pipeline passed |
| SHA-256 of each delivered file, and of the text inside it | the files |
| Sequence count and total length; for the annotation, counts by feature type, gene and transcript counts, the contigs and biotypes present | the files |
| Versions of genomepy, Python, `bgzip` and `samtools`; the pipeline version | the task's environment |
| When it was retrieved, in UTC | the clock |

Two checksums are kept for each file because they answer different questions. `file_sha256` is the compressed file exactly as stored, and changes if it is recompressed. `content_sha256` is the text inside, and is the one you can check: refgenie stores the genome decompressed, so this must match it.

```bash
sha256sum "$(refgenie seek human/fasta)"     # equals genome.content_sha256 in the record
```

Two limits to read the annotation's checksum with:

- **It is the GTF as delivered, not as published.** genomepy removes the header lines and matches contig names to the genome, so the checksum cannot be compared with Ensembl's `CHECKSUMS` or GENCODE's `MD5SUMS`. It identifies what is in your vault and lets two vaults be compared.
- **A UCSC annotation has no published file at all.** It is generated at download time from UCSC's database tables. The record's retrieval time and table date are the only version it has.

The counts are there because a checksum that differs tells you nothing else. Two annotations of the same genome with different gene and transcript counts tell you how they differ.

If part of the record cannot be collected, the download still succeeds and the missing part is named under `provenance_errors`.

`<outdir>/provenance/index.tsv` is a plain-text table of the whole vault, one line per record: each genome's digest and aliases, its provider, assembly, annotation release and accession, when it was retrieved, the two content checksums, and every asset and tag it holds. It is rewritten from the vault after every registration, so it answers "what is in here?" without refgenie:

```bash
column -t -s $'\t' /path/to/outdir/provenance/index.tsv | less -S
```

### The vault as an RO-Crate (`ro-crate-metadata.json`)

The records above are written for this pipeline. `<outdir>/ro-crate-metadata.json` says the standard part of them in a vocabulary other tools already read: [RO-Crate](https://www.researchobject.org/ro-crate/), which is schema.org as JSON-LD. A repository, a FAIR assessment or someone else's script can find out what the vault holds and what made it without knowing anything about refdata-flow.

| In the crate | Is |
|---|---|
| The vault | the root dataset |
| Each genome | a `Dataset`, `data/<digest>/`, with its digest as identifier and its aliases as alternate names |
| Each asset at each tag | a `Dataset`, `data/<digest>/<asset>/<tag>/`, linked to the assets it was built from |
| Each file | a `File`, with its size, its format, and a checksum where one was recorded |
| Each database at each version | a `Dataset`, `databases/<name>/<version>/`, with its file if it is kept as one |
| Each download, each transcriptome derivation, each resource fetch | a `CreateAction`: the tool and its version, what it read, what it produced, when it finished |

An asset that is a whole tree, such as an unpacked VEP cache or Kraken 2 index, is described as its directory and not file by file.

It does not replace the `provenance.json` files. They stay as they are, each is listed in the crate as a file describing its asset, and they remain the only place for what has no standard term: the summary of an annotation, what the check of a resource against its genome found.

Three things to know when reading it:

- **It is rebuilt from the vault after every registration and never edited**, like the index. Only files that are on disk are listed, and an unchanged vault gives the same file byte for byte.
- **A missing `sha256` means "not recorded", not "not checked".** A checksum is given only where the pipeline recorded one that applies to the stored bytes: the genome FASTA, the annotation as delivered, and a resource stored exactly as published. Nothing is hashed to fill the field in. A genome registered before the download record held checksums has none, nor source addresses or a retrieval time, and `--refetch` does not add them: the vault keeps a genome's files and its record from the download that first supplied them. Only a vault built afresh gets the full record for such a genome.
- **It declares RO-Crate 1.1 with the Process Run Crate profile 0.5.** Both have newer editions. These are the ones the reference validator can check, and the crate passes every required check of both:

```bash
pip install roc-validator          # not a dependency of the pipeline; any environment will do
rocrate-validator validate -p process-run-crate-0.5 -l required /path/to/outdir
```

At `-l recommended` the validator reports more, and none of it is fixed by guessing. It asks for things the pipeline does not know: who ran each step (`agent`), who authored and publishes the vault, an end time for steps that recorded none, and a web address that identifies each tool at its exact version. It also reports each step's `actionStatus`, which the crate writes as a link, the way the RO-Crate specification's own examples do, and which the validator compares with a piece of text.

> **Which variant of the assembly you get.** The pipeline takes genomepy's defaults and does not currently expose a switch for either: sequence is **soft-masked** (repeats in lowercase) where the provider offers a masked file, and for providers that distinguish them it is the **primary assembly** rather than toplevel, so alt haplotypes and patch scaffolds are excluded. For Ensembl human that means `GRCh38.p14` names the assembly correctly while the file itself holds 194 sequences: 25 chromosomes plus unplaced scaffolds, no alts. This matters if you intended alt-aware alignment, or if a tool treats lowercase as masked-out. Both settings are recorded in `provenance.json`, under `install_arguments`, with the number of sequences that resulted.

## refgenie vault layout

refgenie organises the genomes under your `--outdir` into a vault keyed by a **sequence-collection digest**, computed over each record's name, length and sequence in file order. Two requests that resolve to the same sequences share a single physical copy.

> **Scope of the digest.** It is computed over the genome FASTA only. The annotation has no digest of that kind; it is identified by its release tag and by the SHA-256 in the download record (see *The download record*).
>
> The digest is *related to* the GA4GH standards without conforming to them, in two independent ways.
>
> **Encoding.** refgenie hashes each sequence with the truncated-SHA-512 construction [refget](https://ga4gh.github.io/refget/) specifies, then hex-encodes it in the legacy `TRUNC512` form. refget's canonical identifier is `SQ.` followed by a base64url `sha512t24u`, and the spec says TRUNC512 "usage SHOULD be discouraged". So the identifier differs in form regardless of the sequence.
>
> **Normalisation.** refget requires sequences be uppercased before hashing; refgenie does not. Where a genome is soft-masked, the digest therefore also differs in value. Masking follows the provider: genomepy requests soft-masked by default, so Ensembl and GENCODE genomes generally are, but not every provider ships a masked file - UCSC `sacCer3`, for instance, arrives entirely uppercase.
>
> The collection digest additionally uses the legacy Henge serialisation (`name>length>sequence_digest`, comma-joined) rather than the form specified by Refget Sequence Collections v1.0.0. Note that it covers each record's **name and length as well as its sequence, in file order** - so the same assembly from two providers, or at two masking levels, is two vault entries. See `refgenconf/seqcol.py` as shipped with refgenie 0.13.0.
>
> Treat the digest as a reliable identifier for *the sequences in this vault*, not as a portable, standards-based fingerprint.

### Annotation release tagging

The digest above covers the genome FASTA and nothing else. Two vaults with the same digest could therefore hold annotations from different releases, and `refgenie seek human/ensembl_gtf` would return a confident path either way. The annotation is tagged to close that gap.

Each provider is tagged with the most specific stable identifier it actually publishes, so the vocabulary differs by design:

| Provider | Tag | Source |
|---|---|---|
| Ensembl | release number, e.g. `116` | the release in the download URL |
| GENCODE | release, e.g. `44` | the release in the download path |
| NCBI | `YYYYMMDD` | the file's `Last-Modified` |
| UCSC | `YYYYMMDD` | when UCSC last changed the annotation table |

A date here is **upstream's** modification date, not the date you ran the pipeline, so re-running next year produces the same tag. The UCSC date comes from the MySQL table genomepy actually reads rather than from the published file dump, because the two disagree.

Ingesting a newer release adds a tag rather than replacing one. The new release becomes what `refgenie seek <genome>/<asset>` returns, and the previous one stays reachable by naming it:

```bash
refgenie seek human/ensembl_gtf:115
```

"Newer" means the higher release, not the one ingested most recently. Pinning release 110 for one analysis when 116 is already in the vault adds `:110` as a tag and leaves the default at 116.

If the release cannot be identified - an unreachable provider, say - the asset is tagged `default` and a warning is printed. Nothing fails, but two releases could later collide under that tag.

### Example Directory Structure

```text
<outdir>/
├── refgenie_export.yaml            <-- The master database index
├── ro-crate-metadata.json          <-- The vault described as an RO-Crate
├── databases/                      <-- Databases that belong to no genome, if you asked for any
│
├── data/                           <-- The vault. Files are named by digest, not by assembly.
│   └── 89c2a7e1.../                <-- sequence-collection digest (truncated here for width)
│       ├── fasta/default/
│       │   ├── 89c2a7e1....fa          <-- note: decompressed by the fasta recipe
│       │   ├── 89c2a7e1....fa.fai
│       │   ├── 89c2a7e1....dict        <-- sequence dictionary, for GATK and Picard
│       │   └── 89c2a7e1....chrom.sizes
│       └── ensembl_gtf/116/        <-- named for the provider recipe; tagged with the release
│           ├── 89c2a7e1....gtf.gz
│           ├── 89c2a7e1..._ensembl_TSS.bed
│           └── 89c2a7e1..._ensembl_gene_body.bed
│
└── alias/                          <-- Human-readable shortcuts. Symlinks are per FILE.
    ├── GRCh38.p14/                 <-- provider's assembly name
    │   ├── fasta/default/
    │   │   ├── GRCh38.p14.fa     -> ../../../../data/89c2a7e1.../fasta/default/89c2a7e1....fa
    │   │   ├── GRCh38.p14.fa.fai -> ...
    │   │   ├── GRCh38.p14.dict   -> ...
    │   │   └── GRCh38.p14.chrom.sizes -> ...
    │   └── ensembl_gtf/116/
    │       ├── GRCh38.p14.gtf.gz -> ../../../../data/89c2a7e1.../ensembl_gtf/116/...
    │       └── GRCh38.p14_ensembl_TSS.bed -> ...
    │
    ├── GRCh38/                     <-- assembly stem
    ├── GRCh38_p14/                 <-- the same name with "." as "_"
    ├── GRCh38_p14_ensembl/         <-- the vault name: <assembly>_<provider>, always present
    └── human/                      <-- your custom alias from the CSV
                                        (each with the same fasta/ and ensembl_gtf/ trees)
```

> **The same assembly from two providers.** Ensembl and GENCODE both publish `GRCh38.p14` and `GRCm39`, with different contig names. Each genome is built under `<assembly>_<provider>` (`GRCm39_ensembl`, `GRCm39_gencode`), which is always unambiguous. The bare names (`GRCm39`) go to whichever was ingested first and stay there; the second is logged as reachable under its vault name. Use the vault name, or your own alias, whenever you hold both.

Three things this layout is easy to get wrong:

- **Every alias points straight at the digest directory.** Aliases are not chained through one another, so removing one does not affect the rest.
- **The symlinks are per file, not per directory**, and each is renamed to the alias. So `alias/human/fasta/default/` contains `human.fa`, and a tool given that path sees a sensible filename rather than a 48-character digest.
- **The FASTA is stored uncompressed.** refgenie's `fasta` recipe runs `gzip -df` before indexing, so budget for the decompressed size in the vault on top of the retained `.fa.gz` download.

### Sequence dictionary (GATK, Picard)

Every genome is built with a sequence dictionary, made by `samtools dict`. It is part of the `fasta` asset, so it sits beside the FASTA and its index under every alias, with the same base name:

```text
alias/human/fasta/default/human.fa
alias/human/fasta/default/human.fa.fai
alias/human/fasta/default/human.dict
```

That placement is the point. GATK and Picard do not take the dictionary's path as an argument: given `-R human.fa` they look for `human.dict` in the same directory. So the path from `refgenie seek` works as a reference directly:

```bash
gatk HaplotypeCaller -R "$(refgenie seek human/fasta)" ...
refgenie seek human/fasta.dict        # the dictionary itself, if you need its path
```

The dictionary describes the sequence and nothing else: contig names, lengths and MD5s. Its `UR` field names the FASTA file beside it, not a path on this machine, so the same sequence gives the same dictionary wherever the vault lives.

To build it, the genome is registered with `recipes/fasta.json`, which is refgenie's stock `fasta` recipe plus that one command. The FASTA, its index and the chromosome sizes file are produced exactly as before.

> **Genomes already in your vault get one on their next ingest.** A genome registered before this existed has no dictionary. The next time it is ingested from a download, its `fasta` asset is rebuilt once to add it; for a human genome that is a few minutes. The sequence digest, aliases and annotations are unchanged. A genome whose download is no longer available is left as it is, with a warning.

### Transcriptome

Every run also gives each requested genome that has an annotation a transcriptome, as the asset `fasta_txome`. There is nothing to ask for; it is part of a run. To leave it out, pass `--no-txome`:

```bash
bash run_pipeline.sh requests/requests_raw.csv --no-txome
```

```bash
refgenie seek human/fasta_txome          # transcripts for the newest annotation release held
refgenie seek human/fasta_txome:115      # for a particular release
```

It is **derived, not downloaded**: `gffread -w` cuts the transcripts out of the vault's own genome using the vault's own GTF. That is a deliberate choice over fetching a provider's cDNA file:

- **Every transcript in the annotation is there.** Ensembl's cDNA file leaves out non-coding genes, so a transcriptome built from it has no lncRNAs, miRNAs or snoRNAs.
- **Every transcript comes from a contig the genome has.** Provider cDNA is cut from an assembly with more contigs than the primary assembly this pipeline stores. That matters for a decoy-aware salmon index, where the transcriptome and the genome must be the same set of sequences.
- **It works for every provider**, including UCSC, which publishes no transcriptome at all.

What you get, and what to know about it:

- **The tag is the annotation release.** `fasta_txome:116` was derived from `ensembl_gtf:116`, and the vault records both the genome and that GTF as the asset's parents. Hold two annotation releases and you get two transcriptomes.
- **Sequence is upper case.** The genomes here are soft-masked, so the transcripts would otherwise carry lower-case repeats that no provider's cDNA has.
- **Nothing is filtered.** Restricting to protein-coding, or dropping retained introns, is an analysis decision and is left to the analysis.
- **Headers are transcript IDs**, as in the GTF, followed by `CDS=start-end` for coding transcripts.
- **`<alias>.provenance.json`** in the asset records the gffread version, the command, the transcript count and when it was derived. A transcriptome derived before the time was recorded has none, and keeps its record; `refgenie remove <alias>/fasta_txome:<release> -f` followed by an ordinary run derives it again with one.
- **The file is called `<alias>.fa`**, the same name as the genome's FASTA one directory over (`alias/human/fasta_txome/116/human.fa` and `alias/human/fasta/default/human.fa`). Tell them apart by the directory.

It runs after registration, so it also covers genomes already in the vault: a vault built before transcriptomes existed, or with `--no-txome`, gains them on its next ordinary run. Nothing is downloaded for that, and only missing transcriptomes are derived. A genome requested with `annotation=false` gets none. If an annotation's contig names do not match its genome, the derivation fails saying so, and nothing is registered.

No aligner or quantifier indexes are built. An index depends on a tool version and on choices such as k-mer size, which belong to an analysis, not to a reference.

### Small derived assets

Every run also makes the small files that pipelines ask for alongside a genome and its GTF, and that otherwise get regenerated, slightly differently, in each project. Like the transcriptome they are derived from what the vault already holds, need no download, and are part of a run unless you pass `--no-derived`.

| Asset | File | What it is | Read by |
|---|---|---|---|
| `tx2gene` | `<alias>.tx2gene.tsv` | transcript_id, gene_id, gene_name; tab-separated, no header | tximport, Salmon's gene-level summary |
| `gtf_indexed` | `<alias>.sorted.gtf.gz` + `.tbi` | the same GTF, sorted by position, bgzipped, tabix-indexed | VEP `--gtf`, IGV, pysam |
| `bed12` | `<alias>.bed` | one BED12 line per transcript | RSeQC, UCSC tools |
| `refflat` | `<alias>.refFlat.txt` | the same gene models in refFlat format | Picard `CollectRnaSeqMetrics` |
| `twobit` | `<alias>.2bit` | the genome in UCSC's `.2bit` format | deepTools, BLAT, rtracklayer |

```bash
refgenie seek human/tx2gene              # for the newest annotation release held
refgenie seek human/tx2gene:115          # for a particular release
refgenie seek human/gtf_indexed.tbi      # the index beside the sorted GTF
refgenie seek human/twobit
```

The four made from an annotation are tagged with its release and record it as their parent, exactly as the transcriptome does. `twobit` is made once per genome, including one requested with `annotation=false`, and records the `fasta` asset as its parent. Each has an `<alias>.provenance.json` with the tool, its version, the command, and when it was made.

What to know before using them:

- **The transcript names are the transcriptome's.** `tx2gene` takes each `transcript_id` as the GTF writes it, which is also what `fasta_txome` names its records by. On Ensembl human 116, all 646,577 names match.
- **`tx2gene` never has an empty gene name.** Where a gene has none, the third column repeats the gene_id. That is 35,483 of the 78,941 genes in Ensembl human 116.
- **`refflat` leads with the gene_id, not the gene name.** Picard groups transcripts into genes by the first column, and gene names are not unique: the same release has 758 genes called `Y_RNA`. Use `tx2gene` to turn an id into a name.
- **In `refflat`, a gene that is in more than one place is named with its place.** Picard skips, whole and without a warning, any gene whose transcripts are not all on one sequence and one strand. Ensembl and GENCODE gene_ids never are in two places. RefSeq's gene_id is the gene symbol, and in NCBI's GRCh38.p14 annotation 1,899 of 43,897 genes are: the same gene on X and on Y, or on a chromosome and on a patch scaffold. Those are written `ABO@9` and `ABO@HG2030_PATCH`, with the strand added where one sequence has the gene on both (`GENE@7(+)`). For those genes the first column is therefore not the gene_id that `tx2gene` has; strip from the `@` to get it. The asset's `provenance.json` counts them.
- **`bed12` follows the UCSC convention for non-coding transcripts**: thickStart equals thickEnd. It is made with UCSC's `gtfToGenePred` and `genePredToBed`, not `gffread --bed`, which marks a whole non-coding transcript as coding.
- **In a UCSC or NCBI annotation, a RefSeq accession placed more than once is renamed from the second place on.** RefSeq puts some transcripts at several places: a pseudo-autosomal gene is on X and on Y, and many are repeated on alternate sequences. genomepy builds these GTFs with UCSC's `genePredToGtf`, which keeps the accession for the first place and adds a suffix after that: `NM_000451.4` on chrX, `NM_000451.4_2` on chrY. In UCSC's hg38 `ncbiRefSeq` that is 11,456 of 197,654 transcripts. They are separate transcripts in every asset here, the transcriptome included, so the names agree with each other; they are not all plain accessions.
- **An annotation that does repeat one transcript ID at two places**, which one given by address may, gets a `bed12` and a `refflat` line for each place under that ID, and one `tx2gene` line. The asset's `provenance.json` counts them. Left alone, `gtfToGenePred` stops at the first such ID.
- **`gtf_indexed` has the provider's records, every one and unchanged; only their order differs.** The sort is stable, so a gene stays ahead of its transcripts and a transcript ahead of its exons. Where order does not matter, the provider's own asset (`ensembl_gtf` and its like) is still the annotation as delivered.
- **`twobit` is checked, not assumed.** After writing it, the sequence names and lengths are read back out of the `.2bit` and must equal the genome's own, or nothing is registered. Soft-masking is kept.

On the human genome the four annotation assets take about three and a half minutes and 2.4 GB of memory in one task. Only what the vault lacks is registered, so removing one asset and running again rebuilds that one.

### Running on an SGE cluster

`--profile sge` adds a scheduler profile for the UCL Research Computing clusters (Myriad, Kathleen, Young) on top of conda:

```bash
bash run_pipeline.sh requests/requests_raw.csv --profile sge
```

Only the derivations are submitted as jobs. Everything else stays on the machine you launched from, for a reason each:

| Step | Runs | Why |
|---|---|---|
| Downloads (`FETCH_GENOME`, `FETCH_RESOURCE`, `FETCH_DATABASE`) | where you launched | need outbound network, which compute nodes often lack |
| Preflight, registration, the check of a resource, the filing of a database | where you launched | they read and edit the vault's single config, or write beside it |
| Transcriptome (`DERIVE_TXOME`) | compute node | needs neither; 1 core, 4 GB, 2 h |
| Small derived assets (`DERIVE_ANNOTATION`, `DERIVE_TWOBIT`) | compute node | need neither; 1 core, 6 GB and 8 GB, 2 h |

So the downloading and compressing, which is the slow part, still happens on the launching machine. So does the renaming of dbSNP and the unpacking of a VEP cache or a Kraken 2 index, each of which reads tens of gigabytes. If your site does not allow that on a login node, launch the pipeline from an interactive session or a transfer node that has network access.

The scheduler request format (memory per core, node-local scratch) is the same as rnaseq-flow's `conf/sge.config`. For another site, copy `conf/sge.config` and adjust it. The profile has been checked to load and resolve; it has **not** yet submitted a job to a real scheduler.

### Benefits for HPC Research
1. **Zero Duplication:** If two researchers run an analysis (one explicitly asking for `GRCh38.p14`, the other asking for `grch38`), they are both funnelled through symlinks to the *exact same physical file* in the `data/` vault, rather than each maintaining a private copy. The alias layer is what makes that work: one download stays reachable under every name people naturally reach for, so nobody re-downloads a genome simply because their preferred name is not in the vault.
2. **Sequence-level Provenance:** Researchers can query refgenie to get the sequence digest (e.g., `2230c53b...`), which identifies the exact set of DNA sequences used, so a later run can be checked against it. The digest covers the FASTA, not the annotation, and does not by itself say which upstream release the data came from; the annotation's release tag and the download record beside it do.

### Using the Vault (querying with refgenie)

The pipeline writes its master config to `refgenie_export.yaml` inside your `--outdir` (default `data/references/refgenie_export.yaml`). You must point `refgenie` at this file, either per-command or via an environment variable.

If you ran the pipeline without installing `refgenie` system-wide, the copy it created lives in the project's virtual environment at `.venv/bin/refgenie` (or run `source .venv/bin/activate` to put `refgenie` on your `PATH`).

**Option A - point at the config per command with `-c`:**
```bash
.venv/bin/refgenie list -c data/references/refgenie_export.yaml
.venv/bin/refgenie seek -c data/references/refgenie_export.yaml GRCh38.p14/fasta
```

**Option B - set the `REFGENIE` environment variable once (recommended):**
```bash
export REFGENIE=/absolute/path/to/data/references/refgenie_export.yaml
refgenie list
refgenie seek GRCh38.p14/fasta
```
Add the `export` line to your `~/.bashrc` to make it permanent (see `walkthrough.md` §6, "HPC Configuration"). Any tool built on `refgenconf` picks up the same variable. nf-core pipelines do **not** - they take reference paths as parameters, or use iGenomes - so pass them explicit paths from `refgenie seek`.

**Common queries:**
- `refgenie list` - list all genomes and their aliases (`human`, `mm39`, `yeast`, …).
- `refgenie seek GRCh38.p14/fasta` (or `human/fasta`) - print the absolute path to the FASTA in the vault.
- `refgenie seek GRCm39/ensembl_gtf` - print the path to an annotation. **The asset name reflects the provider the annotation actually came from**, so use the one matching your request: `ensembl_gtf`, `gencode_gtf`, `ucsc_gtf` or `ncbi_gtf`, or `url_gtf` for an annotation you gave by address.

> **Annotation asset naming.** Earlier versions filed every non-Ensembl annotation under `gencode_gtf`, which asserted a GENCODE origin the data did not have - GENCODE publishes human and mouse only, so a UCSC yeast annotation was recorded as GENCODE. Annotations are now registered under a recipe named for their real source. `gencode_gtf` is refgenie's built-in, used unchanged. `ucsc_gtf`, `ncbi_gtf` and `ensembl_gtf` are shipped in `recipes/` and supplied via `refgenie build --recipe`.

> `ucsc_gtf` and `ncbi_gtf` deliberately only store the GTF, deriving nothing: the TSS and gene-body transformations in the Ensembl recipe assume Ensembl attribute ordering and would corrupt a UCSC or NCBI file.

> **`ensembl_gtf` diverges from refgenie's built-in, on purpose.** The stock recipe pipes its derived BEDs through `sed 's/^/chr/'`, because refgenie's Ensembl recipes are written for pipelines that pair an Ensembl GTF with a UCSC- or GENCODE-named genome. refdata-flow takes the FASTA and the GTF from the same provider, so an Ensembl genome has contigs named `1`/`2`/`X`/`MT` and that rewrite produced `ensembl_tss` and `ensembl_gene_body` BEDs referencing `chr1`/`chr2` - matching nothing in the vault, and silently returning no results from tools like `bedtools getfasta`. The recipe in `recipes/ensembl_gtf.json` is refgenie's, with those two `sed` calls removed and nothing else changed; the GTF asset itself is identical, and the BEDs differ from upstream only in the contig-name column.
>
> A vault built before this change may hold a `gencode_gtf` asset that is actually UCSC- or NCBI-sourced. New runs create the correctly named asset and leave the old one in place; remove it manually if you want the mislabelled entry gone.

### Securing the Vault on HPC

To stop users from editing or deleting the master config and the genome vault, the protection has to be enforced by **ownership and permissions**, not by the config alone. Two points to understand first:

- On Unix you cannot protect a file from the account that *owns* it - an owner can always change the permissions back. So the config/vault must be owned by a **separate curator/service account (or admin group)**, not by ordinary users.
- *Deletion* is governed by the **parent directory's** write bit, not the file's. A user who can write to the containing directory can `rm` a read-only file they don't own - so you must lock down the directory, not just `refgenie_export.yaml`.

**Baseline - owned by a curator account, read-only to everyone else:**
```bash
# run as the curator/service account (or via your HPC admins)
chown -R refdata-svc:hpcusers /shared/reference_datasets
find /shared/reference_datasets -type d -exec chmod 0755 {} \;   # traverse+read, no write for group/other
chmod -R a-w /shared/reference_datasets                          # remove write everywhere
chmod 0444  /shared/reference_datasets/refgenie_export.yaml
```
Users can `refgenie list`/`seek` but cannot edit, rename, or delete anything.

**Finer control with POSIX ACLs** (if the filesystem supports them):
```bash
setfacl -R -m g:hpcusers:rX /shared/reference_datasets
```

**Strongest file-level lock - the immutable flag (with a caveat):**
```bash
sudo chattr +i /shared/reference_datasets/refgenie_export.yaml   # cannot be changed/deleted until -i
```
This is reliable on local **ext4/xfs**, but **most parallel HPC filesystems (Lustre, GPFS, BeeGFS) do not support `chattr +i`** - confirm with your admins before relying on it.

**Best practice:** host the vault in your centre's curated, admin-owned, read-only reference-data area (often exported read-only and surfaced via `module load`), so permissions are enforced centrally and the `REFGENIE` variable can be set in a shared module file.

Two operational notes:

- **Updates must run as the owner.** Once the vault is read-only, re-running the pipeline (Step 3 edits the YAML and removes temp dirs) will fail for anyone else. The intended pattern: the curator account runs `run_pipeline.sh` to add genomes, then re-applies the lock.
- **Protection is not backup.** Permissions prevent edits, not corruption from a bad run. The YAML is tiny - keep it under version control or filesystem snapshots, and optionally record a `sha256sum refgenie_export.yaml` alongside it so tampering is detectable.

## Version History

- **Unreleased**: Correctness fixes from the 2026-11-10 review. Behaviour changes - read before upgrading.
  - **Variant-calling resources.** The six GATK known-sites files for hg38, dbSNP build 157 and the Ensembl VEP cache can be requested in the `resources` column. dbSNP is renamed from RefSeq accessions to the genome's own sequence names through NCBI's assembly report before it is checked; a VEP cache follows the release of the genome's annotation and is unpacked in the vault. See *Variant calling*.
  - **Databases that belong to no genome** (Swiss-Prot, Pfam, SILVA, Kraken 2 indexes) are listed in `requests/databases.csv` and filed under `<outdir>/databases/<name>/<version>/`, with a `databases.tsv` catalogue and an entry each in the RO-Crate. See *Databases that belong to no genome*.
  - **A resumed run asks each publisher again.** A fetch step that had answered "the vault's copy is current" was replayed from Nextflow's cache on every later run for as long as `work/` survived, so a new ClinVar release was not seen until the cache was cleared. Once a genome holds some version of a resource, its fetch step now runs on every run.
  - **A fetch checks the space free before it downloads**, compares what arrived with the size announced, and can verify against a checksum sent with the download as well as one posted in a file. A resource stored exactly as fetched is linked into place in `work/` instead of copied.
  - **Vault genomes are named `<assembly>_<provider>`.** Previously the assembly name alone, which Ensembl and GENCODE share for `GRCh38.p14` and `GRCm39`: the second provider's annotation was filed under the first provider's genome. Bare assembly names remain as aliases, first come first served. A vault built earlier picks up the new name when its download directory is re-ingested; if that directory was already removed with `--cleanup`, the row is reported with "no genome named ...".
  - **Two releases of one assembly can be requested together.** A pinned request downloads into `<provider>/<species>/release-<n>/<assembly>/`, and the resolver cache is keyed on the release as well.
  - **The default annotation tag only moves forward.** Ingesting an older release no longer takes the default from a newer one.
  - **Three checks moved to resolve time**: `annotation=true` on an unannotated assembly, two rows resolving to the same download, and one alias on two different genomes. All previously surfaced after downloading, or not at all.
  - **`requests_resolved.csv` gains an `aliases` column**, and `bin/update_refgenie.py` no longer takes the raw request CSV as an argument.
  - **Failed downloads are retried**, twice by default, with at most four downloads running at once. Previously the first failure ended the run and killed every other download in progress. See *Failed downloads are retried*.
  - **Cataloguing moved into Nextflow.** The wrapper's third step is now the `REGISTER` process, which starts when every download has finished and reads each one as a task input instead of from the copy published under `--outdir`. `--cleanup` therefore no longer publishes sources and then deletes them; it does not publish them. `bin/update_refgenie.py` can still be run by hand against an existing download tree, and gains `--source` and `--manifest`.
  - **Every genome has a sequence dictionary**, beside its FASTA under every alias, so `refgenie seek <genome>/fasta` is usable as a GATK or Picard reference directly. It is a fourth file in the `fasta` asset (`recipes/fasta.json`), and genomes already in a vault gain one on their next ingest. See *Sequence dictionary*.
  - **Third-party resources can be requested** in a `resources` column, starting with ClinVar and the ENCODE blacklist. Each is fetched only when its published version changes, verified against the publisher's checksum where there is one, and registered only after its sequence names and coordinates have been checked against the genome's own. See *Third-party resources*.
  - **A genome can be requested by address** (provider `url`), for the ones no provider's catalogue lists: analysis-set genomes, the Broad's assembly38. The row names the genome and gives the file's address; an annotation is given by address too. See *A genome that no provider lists*.
  - **The download record says what was downloaded.** `provenance.json` gains source URLs, the assembly accession and taxon, the masking and alt-contig settings, SHA-256 checksums of each file and of its contents, a summary of the annotation, tool versions and a UTC timestamp. It used to hold the four request fields and a local-time stamp. See *The download record*.
  - **The record lives in the vault**, as the `provenance` asset, tagged by annotation release and reachable by every alias. The NCBI assembly report, which maps each sequence's name between naming schemes, is kept beside it as `assembly_report`. `<outdir>/provenance/index.tsv` lists the whole vault as a plain table.
  - **The vault describes itself as an RO-Crate.** `<outdir>/ro-crate-metadata.json` lists every genome, asset and file, and each download, derivation and fetch as a `CreateAction` (RO-Crate 1.1, Process Run Crate 0.5). The `provenance.json` files are unchanged and are listed in it. The crate and `index.tsv` are now rewritten after every registration; the index used to be written once, before transcriptomes and resources were added, so its `assets` column missed whatever the same run registered afterwards. See *The vault as an RO-Crate*.
  - **Every run derives the small files pipelines ask for** (`--no-derived` to skip): `tx2gene`, `gtf_indexed`, `bed12` and `refflat` from each annotation, tagged with its release, and `twobit` from each genome. Derived from what the vault holds, each with its own provenance record and its parents recorded. See *Small derived assets*.
  - **Every run derives a transcriptome** for each annotated genome (`--no-txome` to skip). Derived from the vault's own genome and GTF with gffread, registered as `fasta_txome` and tagged with the annotation release, with both inputs recorded as parents. Never a provider's cDNA download. See *Transcriptome*.
  - **An `sge` profile (`--profile sge`)** for the UCL clusters, ported from rnaseq-flow. Only the transcriptome derivation is submitted; downloads and registration stay on the launching machine. See *Running on an SGE cluster*.
  - **A preflight asks the vault before downloading.** Requests the vault already satisfies are skipped, so re-runs no longer depend on Nextflow's `work/` cache surviving. An unpinned annotation is fetched again when the provider publishes a newer release, and gains a tag. `--refetch` downloads everything regardless. See *Genomes already in the vault are not downloaded again*.
  - **An alias, once granted, stays.** Re-ingesting a genome used to leave it with only the aliases of the request that ingested it last: `refgenie build` resets a genome's alias list on every call. An alias given by an earlier request file vanished from the config, with its links still on disk. The ingest now records the list before each build and restores it.
  - **Download tasks use the conda environment's own tools.** The wrapper used to put its virtual environment at the front of `PATH` before launching Nextflow, which left it ahead of the conda environment inside every task: downloads ran under the venv's Python and compressed each FASTA with the pysam stand-in for `bgzip`. The wrapper no longer exports it; only the `REGISTER` step uses the venv.
  - **Boolean and numeric parameters given on the command line are read correctly.** Nextflow passes them through as text, and Groovy treats the text `false` as true. This only became reachable with the new `--cleanup` and `--register` parameters, where `--cleanup false` would have switched cleanup on.
  - **`provenance.json` records `release_requested`**, so a download directory can be matched to the request it answers.
  - **Documents.** The user guide, the talk (`refdata-flow deck.html`) and the workflow map are regenerated for the pipeline as it now is. The talk of 2026-05-27 claimed that refgenie's digest is a GA4GH refget hash, that nf-core pipelines read `REFGENIE`, and that the vault keeps genomes compressed; none of that was true, and it is withdrawn.
  - **Tests.** `tests/` holds an offline suite covering the resolver, the ingest against real refgenie, and `main.nf` itself through a new `stub:` block. CI runs all of it and fails if any test is skipped.
- **0.9.1**: Harm-reduction release. Behaviour changes - read before upgrading.
  - **Source data is no longer deleted by default.** Cleanup of the downloaded assembly directories is now opt-in via `run_pipeline.sh --cleanup`, and even then runs only after the FASTA (and the GTF, if annotation was requested) is verified into the vault. Previously it ran unconditionally, including after a partially-failed ingest.
  - **Provenance survives cleanup.** `provenance.json`, genomepy's `README.txt` and `assembly_report.txt` are copied to `<outdir>/provenance/<provider>/<species>/<assembly>/` before any removal. Previously all three were destroyed by the cleanup step that immediately followed their creation.
  - **The resolver cache no longer returns the wrong genome.** It is keyed on `(species, query)` instead of `species` alone, and re-resolves when a request's provider filter conflicts with the cached result. Previously, editing `query` for an already-resolved species (e.g. `GRCh38` → `GRCh37`) silently returned the *old* assembly. The resolved CSV gains a `query` column to support this; rows written by an earlier version have no query and are simply re-resolved once.
  - **Failures are reported and are non-zero.** `resolve.py` exits 2 and names every request it could not resolve, and the wrapper stops rather than downloading a partial set. `update_refgenie.py` isolates each assembly (one bad row no longer abandons the rest), treats a missing FASTA or a requested-but-absent GTF as an error rather than a warning, and exits non-zero with a summary.
  - **Conda profile completed.** The `conda` profile previously declared only genomepy, leaving the UCSC annotation tools, `samtools` and `htslib` as implicit transitive dependencies - or, outside conda, simply absent. All are now named explicitly, and `conda.createTimeout` is raised to 60 minutes because the 20-minute default can expire mid-solve on a spec this size. Prerequisites now state that conda is required and which distribution to use.
  - **Annotation preflight.** `download.py` now verifies the UCSC annotation tools are on `PATH` *before* downloading, and exits 3 with the exact install command if they are missing. Previously genomepy only discovered this at the annotation step - after the genome FASTA had been fetched and filtered, which for a human assembly is several GB and many minutes wasted. See the note under Prerequisites.
  - **Licensing:** the project is now released under the **GNU General Public License v3.0 or later** - anyone who distributes refdata-flow or a modified version must make their source available under the same terms. `CITATION.cff` is populated so the tool can be cited.
  - **Docs:** corrected the "transcriptomes" claim (the pipeline fetches FASTA + GTF), scoped the sequence-digest and reproducibility claims to what is actually computed, fixed `genome_config.yaml` → `refgenie_export.yaml`, and documented that provider/species directory names use the CSV casing verbatim.
- **0.9beta20260625**: Bug fix, reproducibility, and CI release.
  - **Bug fixes (`main.nf`):** two issues seen running under recent Nextflow are resolved. (1) The input-validation helper `safe()` was promoted from a local closure to a script-level function, so it resolves under both the legacy and strict syntax parsers (fixes `safe is not defined`). (2) The `publishDir` path is now a closure rather than a plain string, so the `provider`/`species` input variables resolve in the task scope (fixes the runtime `No such variable: provider` error; see nextflow-io/nextflow#2811).
  - **Reproducibility:** the Python toolchain is now pinned and hash-locked in `requirements.txt` (genomepy 0.16.4, refgenie 0.13.0, pysam 0.24.0, pyyaml 6.0.2, plus all transitive dependencies); `requirements.in` is the source for regeneration. `run_pipeline.sh` installs with `pip --require-hashes`, adds a Python >= 3.10 preflight check and a virtualenv success marker, and provides `samtools`/`bgzip` (required by `refgenie build` and genomepy) via pysam-backed shims placed on the venv `PATH`. The `nextflow.config` conda and docker profiles are pinned to genomepy 0.16.4 to match.
  - **CI:** added a GitHub Actions workflow (`.github/workflows/smoke-test.yml`) that runs on every push/PR - a `--dry-run` smoke test (with static shell/Python checks and a pre-seeded resolver cache so it needs no network) plus a `nextflow lint` job (Nextflow 26.04.4, strict parser) that statically validates `main.nf` and `nextflow.config`. Sample fixtures live in `examples/`.
  - **Usability:** new `run_pipeline.sh --debug` flag streams live per-genome progress messages to the console (via the process `debug` directive and phase logging in `download.py`); off by default.
  - **Repo tidying:** `.gitignore` now excludes `.venv/`, `work/`, `.nextflow/`, `.nextflow.log*`, and `data/`; added placeholder `LICENSE` and `CITATION.cff` files (to be completed after Zenodo deposit).
- **0.9beta20260604**: Security hardening - shell injection guard via Groovy input validation in `main.nf`; path traversal protection and safe directory removal with full error logging in `update_refgenie.py`; GTF build failures now surface as errors; structured append-mode logging to `update_refgenie.log`.
- **0.9beta20260527**: Initial complete beta release featuring full refgenie 0.13 integration (config schema 0.4), automated configuration mapping, sequence-collection digests, and direct YAML-based alias management.
