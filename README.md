<p align="center">
  <img src="refdata-flow_logo.svg" alt="refdata-flow logo" width="600">
</p>


This project provides an automated workflow for downloading bioinformatics reference datasets - genomic DNA sequence (FASTA) and gene annotation (GTF) - for multiple species, and cataloguing them in a [refgenie](http://refgenie.databio.org/) vault.

It runs in **three stages**, wrapped by a single script:
1. **Resolve** - an interactive Python step (`resolve.py`) turns generic queries (like `grch38` or `mouse`) into exact assembly requests.
2. **Download** - a **Nextflow** pipeline fetches the files in parallel via the **genomepy** engine, suitable for HPC compute nodes.
3. **Catalogue** - `update_refgenie.py` registers the results in a refgenie vault with human-readable aliases.

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

## Directory Structure

```
├── bin/
│   ├── download.py          # Core Python script used by Nextflow
│   ├── resolve.py           # Interactive Python pre-flight wizard
│   └── update_refgenie.py   # Post-processing script to automate the refgenie catalogue
├── data/
│   └── references/          # Default output directory for genomes
├── main.nf                  # The Nextflow pipeline
├── nextflow.config          # Nextflow configuration profiles
├── README.md
├── requests_raw.csv         # The configuration file you edit
├── requests_resolved.csv    # Generated by resolve.py, used by Nextflow
└── run_pipeline.sh          # Unified wrapper script to run everything
```

## How to use

### Step 1: Write Generic Requests
Edit the `requests_raw.csv` file to list the generic species you want. 
Example `requests_raw.csv`:
```csv
species,query,provider,annotation,aliases
Homo_sapiens,grch38,ensembl,true,human
Mus_musculus,mouse,,true,mm39;mm
```
*(Note: If you leave `provider` empty like for Mouse, it will search all providers. You can also supply multiple custom aliases separated by semicolons in the `aliases` column.)*

> **Aliases are labels, not conversions.** An alias is a second name for whatever you downloaded; it does not change the data. Aliasing an Ensembl assembly as `hg38` is legal and tempting, but `hg38` is UCSC's name and implies UCSC conventions: `chr1`, `chr2`, `chrM`. An Ensembl genome has contigs named `1`, `2`, `MT`, so `refgenie seek hg38/fasta` would hand a downstream tool exactly the naming it did not expect. Keep provider-specific nicknames for genomes from that provider.

### Step 2: Run the Pipeline Wrapper
Instead of memorizing the Python and Nextflow commands, simply run the provided shell script:

```bash
bash run_pipeline.sh requests_raw.csv --outdir /path/to/your/custom/storage
```
*(If you omit `--outdir`, it defaults to `data/references`)*

This script will automatically:
1. **Resolve**: Parse your `requests_raw.csv`. It connects to the databases and interactively prompts you if there are multiple matches for your generic queries.
2. **Download**: Automatically launch the Nextflow pipeline on the resolved targets using Conda. It always passes `-resume`, so genomes that completed in a previous run are skipped (see [Resuming an Interrupted Run](#-pro-tip-resuming-an-interrupted-run) for the limits).
3. **refgenie catalogue**: Finally, it automatically initializes a `refgenie_export.yaml` for you (if one doesn't exist yet), computes the sequence digests, and links your aliases. The downloaded source directories are **kept** unless you pass `--cleanup` (see below).

### 💡 Pro-Tip: Skipping the Wizard & Caching
The interactive wizard remembers your choices. If you run the pipeline again to add a new species, it will look at `requests_resolved.csv`. If an older species has already been resolved, **it will automatically skip the slow network search and interactive prompt** for that species to save time. 

If you ever *want* to change the assembly for a species, simply delete that specific row from `requests_resolved.csv` (or delete the file entirely) and the wizard will prompt you again.

### 🗄️ Reclaiming space: `--cleanup` (opt-in)

Once a genome has been ingested into the refgenie vault, its original download directory is redundant. Passing `--cleanup` removes it:
```bash
bash run_pipeline.sh requests_raw.csv --cleanup
```
Two safeguards apply, because deleting downloaded reference data is not always reversible: providers archive and reorganise old releases, so re-fetching the identical files months later is not guaranteed:

- **Nothing is removed unless the ingest is verified.** The FASTA must have been registered, and if you requested annotation, the GTF must have been registered too. A partial ingest leaves the source directory in place and the run exits non-zero.
- **Provenance survives cleanup.** Before removal, `provenance.json`, genomepy's `README.txt` (source URLs) and, where the provider supplies an assembly accession, `assembly_report.txt` (sequence-name mappings between UCSC / Ensembl / GenBank accessions) are copied to `<outdir>/provenance/<provider>/<species>/<assembly>/`, which is never deleted. The remaining files are the sequence itself and indexes that refgenie regenerates in the vault.

Prior to 0.9.1 this cleanup ran unconditionally, including on partially-failed ingests, and took the provenance files with it.

### 🧹 Pro-Tip: Storage Optimization
Nextflow can generate a lot of background cache in the `.nextflow/` and `work/` directories over time. To instantly purge all logs and pipeline cache to recover your disk space, simply run:
```bash
bash run_pipeline.sh --clean
```

### 🔍 Pro-Tip: Dry-Run Mode
If you want to parse your `requests_raw.csv` and resolve exactly which datasets and providers will be used without actually downloading gigabytes of data, run the pipeline in dry-run mode:
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

### 🔁 Pro-Tip: Resuming an Interrupted Run
You do **not** need to pass a resume flag - `run_pipeline.sh` always launches Nextflow with `-resume`, so simply re-running the pipeline continues where the previous run left off:
```bash
bash run_pipeline.sh requests_raw.csv        # automatically resumes the previous session
```
A couple of things to be aware of:

- Resume works at the **per-genome** level: it skips genomes whose download *completed successfully* in a prior run. A genome that was still downloading when you interrupted it did not complete, so it restarts from the beginning (there is no partial/byte-level resume). Re-running is safe: finished genomes are skipped, only unfinished ones re-run.
- Resume relies on the `work/` and `.nextflow/` directories. Running `bash run_pipeline.sh --clean` (or deleting those directories) **erases the resume cache**, so the next run re-downloads everything. Only `--clean` once you're finished.
- On a resumed run, skipped genomes are shown as `cached` in the Nextflow output, e.g. `… FETCH_GENOME (…) [100%] 1 of 1, cached: 1 ✔`.

## Output & Provenance

For each requested assembly, Nextflow will generate an output directory structured by provider, species, and assembly.

Example:
`/path/to/your/custom/storage/ensembl/Homo_sapiens/GRCh38.p14/`

The `provider` and `species` directory names are used **exactly as they appear in your resolved CSV** - so a `provider` of `ensembl` produces `ensembl/`, not `Ensembl/`. Keep the casing consistent between runs, or you will end up with parallel directory trees for the same provider.

Inside each assembly directory, you will find:
- The downloaded compressed `.fa.gz` (and `.gtf.gz` if requested).
- A `README.txt` generated natively by genomepy detailing the source URLs.
- A custom `provenance.json` recording exactly when the dataset was downloaded and what arguments were used, providing an audit trail.

> **Which variant of the assembly you get.** The pipeline takes genomepy's defaults and does not currently expose a switch for either: sequence is **soft-masked** (repeats in lowercase) where the provider offers a masked file, and for providers that distinguish them it is the **primary assembly** rather than toplevel, so alt haplotypes and patch scaffolds are excluded. For Ensembl human that means `GRCh38.p14` names the assembly correctly while the file itself holds 194 sequences: 25 chromosomes plus unplaced scaffolds, no alts. This matters if you intended alt-aware alignment, or if a tool treats lowercase as masked-out. Neither setting is recorded in `provenance.json`.

## refgenie vault layout

refgenie organises the genomes under your `--outdir` into a vault keyed by a **sequence-collection digest**, computed over each record's name, length and sequence in file order. Two requests that resolve to the same sequences share a single physical copy.

> **Scope of the digest.** It is computed over the genome FASTA only - the annotation (GTF) is not currently digested, and neither artifact is checksummed at download time.
>
> The digest is *related to* the GA4GH standards without conforming to them, in two independent ways.
>
> **Encoding.** refgenie hashes each sequence with the truncated-SHA-512 construction [refget](https://ga4gh.github.io/refget/) specifies, then hex-encodes it in the legacy `TRUNC512` form. refget's canonical identifier is `SQ.` followed by a base64url `sha512t24u`, and the spec says TRUNC512 "usage SHOULD be discouraged". So the identifier differs in form regardless of the sequence.
>
> **Normalisation.** refget requires sequences be uppercased before hashing; refgenie does not. Where a genome is soft-masked, the digest therefore also differs in value. Masking follows the provider: genomepy requests soft-masked by default, so Ensembl and GENCODE genomes generally are, but not every provider ships a masked file - UCSC `sacCer3`, for instance, arrives entirely uppercase.
>
> The collection digest additionally uses the legacy Henge serialisation (`name>length>sequence_digest`, comma-joined) rather than the form specified by Refget Sequence Collections v1.0.0. Note that it covers each record's **name and length as well as its sequence, in file order** - so the same assembly from two providers, or at two masking levels, is two vault entries. See `refgenconf/seqcol.py` as shipped with refgenie 0.13.0.
>
> The pipeline also does not pin the upstream provider release, so the same request re-run months later may resolve to newer files. Treat the digest as a reliable identifier for *the sequences in this vault*, not as a portable, standards-based fingerprint.

### Example Directory Structure

```text
<outdir>/
├── refgenie_export.yaml            <-- The master database index
│
├── data/                           <-- The vault. Files are named by digest, not by assembly.
│   └── 89c2a7e1.../                <-- sequence-collection digest (truncated here for width)
│       ├── fasta/default/
│       │   ├── 89c2a7e1....fa          <-- note: decompressed by the fasta recipe
│       │   ├── 89c2a7e1....fa.fai
│       │   └── 89c2a7e1....chrom.sizes
│       └── ensembl_gtf/default/    <-- asset named for the provider recipe, not "annotation"
│           ├── 89c2a7e1....gtf.gz
│           ├── 89c2a7e1..._ensembl_TSS.bed
│           └── 89c2a7e1..._ensembl_gene_body.bed
│
└── alias/                          <-- Human-readable shortcuts. Symlinks are per FILE.
    ├── GRCh38.p14/                 <-- provider's assembly name
    │   ├── fasta/default/
    │   │   ├── GRCh38.p14.fa     -> ../../../../data/89c2a7e1.../fasta/default/89c2a7e1....fa
    │   │   ├── GRCh38.p14.fa.fai -> ...
    │   │   └── GRCh38.p14.chrom.sizes -> ...
    │   └── ensembl_gtf/default/
    │       ├── GRCh38.p14.gtf.gz -> ../../../../data/89c2a7e1.../ensembl_gtf/default/...
    │       └── GRCh38.p14_ensembl_TSS.bed -> ...
    │
    ├── GRCh38/                     <-- assembly stem
    ├── GRCh38_p14/                 <-- refgenie's internal name ("." becomes "_")
    └── human/                      <-- your custom alias from the CSV
                                        (each with the same fasta/ and ensembl_gtf/ trees)
```

Three things this layout is easy to get wrong:

- **Every alias points straight at the digest directory.** Aliases are not chained through one another, so removing one does not affect the rest.
- **The symlinks are per file, not per directory**, and each is renamed to the alias. So `alias/human/fasta/default/` contains `human.fa`, and a tool given that path sees a sensible filename rather than a 48-character digest.
- **The FASTA is stored uncompressed.** refgenie's `fasta` recipe runs `gzip -df` before indexing, so budget for the decompressed size in the vault on top of the retained `.fa.gz` download.

### Benefits for HPC Research
1. **Zero Duplication:** If two researchers run an analysis (one explicitly asking for `GRCh38.p14`, the other asking for `grch38`), they are both funnelled through symlinks to the *exact same physical file* in the `data/` vault, rather than each maintaining a private copy. The alias layer is what makes that work: one download stays reachable under every name people naturally reach for, so nobody re-downloads a genome simply because their preferred name is not in the vault.
2. **Sequence-level Provenance:** Researchers can query refgenie to get the sequence digest (e.g., `2230c53b...`), which identifies the exact set of DNA sequences used, so a later run can be checked against it. Note the limits described above - the digest covers the FASTA, not the annotation, and does not by itself record which upstream release the data came from.

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
- `refgenie seek GRCm39/ensembl_gtf` - print the path to an annotation. **The asset name reflects the provider the annotation actually came from**, so use the one matching your request: `ensembl_gtf`, `gencode_gtf`, `ucsc_gtf` or `ncbi_gtf`.

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
