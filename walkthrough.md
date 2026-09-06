# refdata-flow: Project Walkthrough

The pipeline is now completely operational and ready for HPC deployment! Here is a breakdown of the architecture we built today.

## 1. The Three-Step Architecture

### The 3-Step Execution
1. **Resolution (`resolve.py`)**: Interactively searches Ensembl/UCSC/NCBI and binds `requests/requests_raw.csv` to exact assemblies into a fast, cached `requests/requests_resolved.csv` ledger.
2. **Download (`Nextflow`)**: Mass-parallelizes Conda-based `genomepy` jobs to heavily compress and download the `fasta` and `gtf` files into temporary, provider-specific staging directories.
3. **Refgenie Ingestion (`update_refgenie.py`)**:
    - Sweeps over the Nextflow output directories.
    - Triggers `refgenie build` to parse, index, and compute immutable cryptographic hashes of the sequence data, storing them in the permanent `data/references/data/` vault.
    - **Optional cleanup (`--cleanup`)**: Removes the Nextflow staging downloads once their contents are verified into the vault, to save disk space. Opt-in since 0.9.1: upstream providers archive old releases, so a deleted download may not be re-fetchable, and `provenance.json`, `README.txt` and `assembly_report.txt` are preserved to `<outdir>/provenance/` before anything is removed.
    - **Robust Aliasing**: Reads the custom `aliases` column from your CSV and uses Python's native `yaml` library to instantly and reliably inject all human-readable names (e.g. `grch38`, `human`) into the master config, completely bypassing Windows OS execution quirks!

## 2. Strict Refgenie Hashes (`bin/update_refgenie.py`)

To ensure absolute reproducibility when publishing papers, the pipeline integrates natively with **refgenie 0.13**:
- **Automated Initialization**: It detects if the HPC has a config file yet, and automatically runs `refgenie init` if it's the very first run.
- **Cryptographic Hashing**: It executes `refgenie build` on the downloaded FASTA and GTF files. This forces Refgenie to read every DNA letter, compute the globally-agreed **Refget hash** (e.g. `2230c53b...`), and ingest it into the secure `data/` vault.
- **Human Aliases**: It applies the user's custom CSV queries (like `mouse`) as symlinked aliases so researchers don't have to memorize the hashes.

## 3. Storage Optimization

Storage is expensive. To prevent the pipeline from duplicating the 3GB genomes:
1. `genomepy` downloads the files compressed natively as `.fa.gz` and `.gtf.gz` via the `bgzip` flag.
2. Nextflow uses `-resume` to instantly skip any genome it successfully downloaded in a previous run.
3. After `refgenie build` securely copies the sequences into its cryptographic `data/` vault, the Python script executes `shutil.rmtree` to completely **wipe the Nextflow temporary downloads**, achieving zero disk duplication!

## 4. Dry-Run Verification

Before launching the pipeline to download massive genomic datasets, users can verify exactly what will happen using the `--dry-run` flag. This securely parses the CSV, probes the databases to find exact assembly matches, caches the ledger, and aborts before Nextflow takes over:

```bash
$ bash run_pipeline.sh --dry-run

==========================================
  refdata-flow Pipeline
==========================================

[Step 1] Running interactive reference resolver...
--- Species 'Homo_sapiens' already resolved to GRCh38.p14. Skipping search. ---
--- Species 'Mus_musculus' already resolved to GRCm39. Skipping search. ---
--- Species 'Saccharomyces_cerevisiae' already resolved to sacCer3. Skipping search. ---
Successfully resolved 3 requests. Written to requests/requests_resolved.csv.

==========================================
  DRY RUN COMPLETE
  The following datasets would be downloaded:
==========================================
  - Homo_sapiens: GRCh38.p14 (Provider: Ensembl, Annotations: TRUE)
  - Mus_musculus: GRCm39 (Provider: GENCODE, Annotations: TRUE)
  - Saccharomyces_cerevisiae: sacCer3 (Provider: UCSC, Annotations: FALSE)
==========================================
  Exiting without downloading.
```

## 5. Storage Cleanup Mode

Nextflow pipelines inherently generate significant background metadata over time, including `.nextflow.log.*` files, `.nextflow/` local caches, and the massive `work/` directory where the background jobs process the genomes. 

To safely sweep away all of this pipeline footprint after Refgenie has successfully secured your data, users can simply run the `--clean` flag:

```bash
$ bash run_pipeline.sh --clean

==========================================
  refdata-flow Pipeline
==========================================

Cleaning up Nextflow cache, logs, and work directories...
Clean complete! Storage optimized.
```
This guarantees that your infrastructure stays incredibly lean, storing only the final cryptographic reference vault.

## 6. HPC Configuration (Environment Variables)

To let the `refgenie` CLI locate your master genome vault without an explicit `-c` flag on every command, configure the `REFGENIE` environment variable on your HPC. Any tool built on `refgenconf` reads it too. Note that nf-core pipelines do not: they take reference paths as parameters, or use iGenomes, so pass them explicit paths obtained from `refgenie seek`.

### Step-by-Step Instructions

1. **Locate your Master Config:**
   Identify the absolute path to the `refgenie_export.yaml` file generated by the pipeline. For example: `/mnt/hpc/shared_data/reference_datasets/data/references/refgenie_export.yaml`

2. **Update your Bash Profile:**
   Add the export command to your user or system-wide `~/.bashrc` (or `~/.bash_profile`) file so it loads automatically on every login:
   ```bash
   echo "export REFGENIE=/mnt/hpc/shared_data/reference_datasets/data/references/refgenie_export.yaml" >> ~/.bashrc
   ```

3. **Apply the Changes:**
   Refresh your current terminal session to apply the variable:
   ```bash
   source ~/.bashrc
   ```

4. **Verify the Configuration:**
   Test that Refgenie globally recognizes your vault from any directory:
   ```bash
   refgenie list
   ```
   *If configured correctly, this command will instantly list all of your hashed genomes and custom aliases (e.g., `hg38`, `GRCm39`) without needing to be inside the project folder!*

> [!TIP]
> **System-Wide Modules:** If your HPC uses an Environment Module system (e.g., Lmod), you should append `setenv REFGENIE /path/to/refgenie_export.yaml` to the shared `refgenie` module file so all researchers inherit the vault path instantly when they type `module load refgenie`.

## 7. Generating the User Guide

The pipeline's documentation can be exported as a beautifully formatted Microsoft Word document (`.docx`), generated programmatically via a Node.js script. This ensures the documentation is always perfectly styled without manual Word formatting.

### Build Instructions

You must have [Node.js](https://nodejs.org/) installed on your machine. To build the document, run the following commands from the root of the repository:

```bash
cd docs
npm install          # Only required the first time to install dependencies
npm run build        # Generates the .docx file
```

This will instantly create `refdata-flow_User_Guide.docx` in the root of the repository!
