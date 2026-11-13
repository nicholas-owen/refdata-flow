#!/usr/bin/env nextflow
// Copyright (C) 2026 Nicholas Owen
// SPDX-License-Identifier: GPL-3.0-or-later

nextflow.enable.dsl=2

// Validate a value against a strict allowlist before it is interpolated into a
// shell command. Only alphanumerics, dots, underscores and hyphens are allowed.
// Declared as a script-level function (not a local closure) so the call resolves
// under both the legacy and the strict Nextflow language parsers - the strict
// parser treats `name(args)` as a function call, so a local closure variable is
// reported as "not defined". error() aborts the run immediately, preventing
// shell injection.
def safe(String val, String field) {
    if( !(val ==~ /^[A-Za-z0-9._-]+$/) )
        error "Unsafe characters detected in ${field}: '${val}'"
    return val
}

// An address for a genome requested by URL, or '' when the row has none.
//
// The same rule as URL_RE in bin/resolve.py, which has already applied it; this
// is the guard for a resolved CSV written by hand. An address is interpolated into
// a shell command inside double quotes, so the characters a shell still acts on
// there - $, backtick, double quote, backslash - must not get through, and the
// allowlist is kept narrow enough to read at a glance.
def safeUrl(val, String field) {
    def text = (val ?: '').toString().trim()
    if( text && !(text ==~ /^(https?|ftp):\/\/[A-Za-z0-9._~\/%+=:@,-]+$/) )
        error "Unusable ${field}: '${text}'. It must be an http, https or ftp address using only letters, digits and . _ ~ / % + = : @ , -"
    return text
}

// Read a parameter as a boolean or a whole number, whatever type it arrived as.
//
// A value from nextflow.config is a Boolean or an Integer. The same parameter
// given on the command line is a String (verified against Nextflow 26.04:
// `--cleanup false` arrives as the text 'false'), and Groovy counts every
// non-empty string as true. Used as it comes, `--cleanup false` therefore switched
// cleanup ON, and `--max_retries 10` could not be compared with a number at all.
// Every boolean and numeric parameter in this file is read through these two -
// except in process directives, where a function defined here cannot be called
// (Nextflow reports "Unknown process directive: `isTrue`" at run time, though
// `nextflow lint` accepts it). Directives repeat the same conversion inline.
//
// A directive is also evaluated BEFORE the workflow block runs, so it must not be
// the thing that trips over a mistyped number: `--max_retries lots` would then die
// with a bare 'For input string: "lots"'. The numeric directives fall back to a
// harmless value instead, and the asInt() checks at the top of the workflow block
// stop the run with a proper message before any task is created.
def isTrue(value) {
    return value.toString().trim().toLowerCase() == 'true'
}

def asInt(value, String name) {
    def text = value.toString().trim()
    if( !(text ==~ /^[0-9]+$/) )
        error "--${name} must be a whole number, got '${text}'"
    return text.toInteger()
}

// Refuse a request list in which two rows name the same download.
//
// bin/resolve.py already rejects these, so this only fires for a resolved CSV that
// was written by hand and run with a bare `nextflow run`. It is here because the
// failure is otherwise silent: two FETCH_GENOME tasks with one publish directory
// simply overwrite each other, and when the rows disagree on `annotation` the
// provenance.json that survives can belong to the sequence-only task.
//
// Each row is the tuple built in the workflow below: species, assembly, provider,
// annotation, release, then the two addresses, which play no part here. Provider is compared case-insensitively because genomepy
// reports 'Ensembl' while a hand-written file may say 'ensembl'.
//
// Call this from a `map`, not a `flatMap`. Under Nextflow 26.04 an error() raised
// inside flatMap still stops the run, but reaches the user as "Unexpected error
// [InvocationTargetException]" with the message only in .nextflow.log; from map it
// is printed as written. Found by tests/test_nextflow_stub.py the first time this
// guard was actually executed.
def rejectDuplicates(List rows) {
    def names = rows.collect { row ->
        "${row[1]} from '${(row[2] ?: '').toString().toLowerCase()}' (release: ${row[4] ?: 'current'})".toString()
    }
    def repeated = names.findAll { name -> names.count(name) > 1 }.unique()
    if( repeated )
        error "Duplicate request for ${repeated.join('; ')}. Each provider, assembly and release may appear once; merge the rows."
    return rows
}

process FETCH_GENOME {
    tag "${species}_${assembly}"

    // Needs outbound network, which compute nodes often do not have. A scheduler
    // profile uses this label to keep the task on the machine Nextflow was
    // launched from (see conf/sge.config).
    label 'network'

    // When params.debug is true, stream the task's stdout to the Nextflow console
    // so download.py's phase messages ([1/3] Starting, [2/3] Downloading…,
    // [3/3] Writing provenance…) appear live instead of only the single status
    // line. Enabled via `run_pipeline.sh --debug`. Output from parallel tasks
    // interleaves (each line is tagged with its assembly); for the detailed
    // byte-level progress bar, tail the task's work/<hash>/.command.log.
    debug(params.debug.toString().toLowerCase() == 'true')

    // A download fails for reasons that have nothing to do with the request: a
    // provider that is briefly unreachable, a connection dropped part-way through a
    // multi-gigabyte file. Without a strategy the first such failure ended the run,
    // and Nextflow's default ('terminate') killed every other download still in
    // progress along with it.
    //
    //   * exit 3 is bin/download.py's "a required tool is missing". Trying again
    //     cannot fix that, so it is never retried.
    //   * anything else is retried, up to params.max_retries further attempts.
    //   * when the attempts run out, 'finish' rather than 'terminate': downloads
    //     that are already running complete (and stay cached for -resume) before
    //     the run exits non-zero.
    //
    // The wait before a retry lives here, not as a `sleep` in the script block, so
    // the script - which is part of what -resume matches a finished task on - is
    // the same on every attempt. The cost of sleeping here is that Nextflow's task
    // bookkeeping pauses for those seconds; the downloads themselves carry on.
    errorStrategy {
        if( task.exitStatus == 3 || task.attempt > asInt(params.max_retries, 'max_retries') )
            return 'finish'
        sleep(asInt(params.retry_wait, 'retry_wait') * 1000L * task.attempt)
        return 'retry'
    }
    maxRetries(params.max_retries.toString().trim().isInteger() ? params.max_retries.toString().trim().toInteger() : 0)

    // Every task is a download from one of four public providers. Unbounded, a long
    // request list opens one connection per row at once, which is unkind to them and
    // tends to get throttled.
    maxForks(params.max_downloads.toString().trim().isInteger() ? params.max_downloads.toString().trim().toInteger() : 1)

    // Closure form (not a plain string) so the input variables `provider` and
    // `species` are resolved in the task scope - a bare "${provider}" string
    // raises "No such variable: provider" (see nextflow-io/nextflow#2811).
    //
    // A pinned release gets its own directory level. Without it, two releases of one
    // assembly published to the same <provider>/<species>/<assembly>/ and the second
    // overwrote the first, so the ingest tagged both rows from whichever
    // provenance.json landed last. An unpinned request keeps the original layout.
    // bin/update_refgenie.py builds the same path; keep the two in step.
    //
    // `release` is safe to interpolate here: the workflow block has already
    // required it to be empty or all digits.
    //
    // This copy is the retained source, not the ingest's input: REGISTER reads the
    // task output directly (see below). So `--cleanup true` simply switches the
    // publishing off, where it used to publish and then delete.
    publishDir {
        release ? "${params.outdir}/${provider}/${species}/release-${release}" : "${params.outdir}/${provider}/${species}"
    }, mode: 'copy', enabled: params.cleanup.toString().toLowerCase() != 'true'

    input:
    tuple val(species), val(assembly), val(provider), val(annotation), val(release), val(url), val(annotation_url)

    // The three values travel with the directory so REGISTER can tell which request
    // it answers: a directory called GRCm39 does not say whose GRCm39 it is.
    output:
    tuple val(provider), val(assembly), val(release), path("${assembly}")

    script:
    // Validate inputs against the strict allowlist (see safe() above) before any
    // shell interpolation, to prevent shell injection from CSV-supplied values.
    def sp  = safe(species,  'species')
    def asm = safe(assembly, 'assembly')
    def prv = safe(provider, 'provider')
    // annotation_flag is safe: derived from a ternary with two hardcoded strings,
    // no user-supplied data is interpolated into it.
    def annotation_flag = annotation.toString().toLowerCase() == 'true' ? '--annotation' : ''

    // release is optional, so it cannot go through safe(), whose pattern requires at
    // least one character. Its allowlist is tighter than safe()'s in any case: only
    // Ensembl can be pinned and Ensembl releases are integers, so digits only.
    // bin/resolve.py has already rejected a release on an unpinnable provider; this
    // is the injection guard, not the semantic check.
    def rel = release?.toString()?.trim() ?: ''
    if( rel && !(rel ==~ /^[0-9]+$/) )
        error "Unsafe or non-numeric release: '${rel}'"
    def release_flag = rel ? "--release ${rel}" : ''

    // A genome requested by address (provider 'url'). Both values went through
    // safeUrl() in the workflow block, and are empty for every other provider.
    def url_flag = url ? "--url \"${url}\"" : ''
    def annotation_url_flag = annotation_url ? "--annotation-url \"${annotation_url}\"" : ''

    // download.py is called without an explicit path - Nextflow automatically adds
    // the project bin/ directory to PATH inside every process (fixes issue 3.1).
    """
    download.py \\
        --species "${sp}" \\
        --assembly "${asm}" \\
        --provider "${prv}" \\
        --pipeline-version "${workflow.manifest.version}" \\
        ${annotation_flag} \\
        ${release_flag} \\
        ${url_flag} \\
        ${annotation_url_flag}
    """

    // `nextflow run -stub` runs this instead of the script above: no network, no
    // genomepy, no conda. It writes the same directory download.py leaves behind,
    // at toy size, so the parts of this workflow that only ever ran during a real
    // download - the request checks, the duplicate guard, the release-level publish
    // path - can be executed in CI in seconds. tests/test_nextflow_stub.py drives it
    // and hands the result to bin/update_refgenie.py.
    //
    // Contigs are chr-prefixed for the providers that name them that way. That is
    // not decoration: one assembly requested from Ensembl and from GENCODE must
    // yield two different sequences, as it does in reality, or refgenie computes one
    // digest for both and the test of keeping them apart tests nothing. The third
    // contig carries the assembly name for the same reason, one level up: two
    // different assemblies must not come out as the same sequence either.
    //
    // provenance.json is marked "stub": true so the output can never be mistaken
    // for a record of a real download.
    stub:
    def sp  = safe(species,  'species')
    def asm = safe(assembly, 'assembly')
    def prv = safe(provider, 'provider')
    def want_gtf = annotation.toString().toLowerCase() == 'true'
    def chr = ['ucsc', 'gencode'].contains(prv.toLowerCase()) ? 'chr' : ''
    // An unpinned request stands for "the current release". 999 plays that part:
    // a number, like every real release, and above any a test would pin, so the
    // newest-release-is-the-default rule behaves here as it does on real data.
    def rel_tag = release ?: '999'
    def attrs = 'gene_id "g1"; gene_version "1"; gene_name "G1"; gene_source "stub"; gene_biotype "protein_coding";'
    def exon_attrs = 'gene_id "g1"; gene_version "1"; transcript_id "t1"; transcript_version "1"; exon_number "1"; gene_name "G1"; gene_source "stub"; gene_biotype "protein_coding";'
    def release_json = want_gtf ? ", \"annotation_release\": \"${rel_tag}\", \"annotation_release_kind\": \"stub\"" : ''
    """
    mkdir -p "${asm}"
    printf '>${chr}1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n>${chr}2\\nTTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAA\\n>${chr}Un_${asm}\\nGATTACAGATTACAGATTACA\\n' \\
        | gzip > "${asm}/${asm}.fa.gz"
    if [ "${want_gtf}" = "true" ]; then
        printf '${chr}1\\tstub\\tgene\\t1\\t30\\t.\\t+\\t.\\t${attrs}\\n${chr}1\\tstub\\texon\\t1\\t30\\t.\\t+\\t.\\t${exon_attrs}\\n' \\
            | gzip > "${asm}/${asm}.annotation.gtf.gz"
    fi
    printf '%s\\n' '{"stub": true, "species": "${sp}", "assembly": "${asm}", "provider": "${prv}", "release_requested": "${release ?: ''}", "annotation_requested": ${want_gtf}${release_json}}' \\
        > "${asm}/provenance.json"
    """
}

// Ask the vault which requests it already satisfies, so only the rest are fetched.
//
// Until this existed the only thing that kept a finished genome from being
// downloaded again was the resume cache in work/ - which `--clean` deletes and an
// HPC does not keep. The vault's config is kept, and already records what has been
// built, so it is consulted first. bin/preflight.py has the rules; the short
// version is that a genome is skipped when it and the requested annotation release
// are really there, and fetched whenever that is in doubt.
//
// Emits the keys of the requests to download, one per line, in the same format
// REGISTER uses to name its sources.
//
// `ready` is not used by the script. It is the number of rows that passed the
// checks in the workflow block, and taking it as an input is what makes this task
// wait for those checks: without it the preflight could start reading a request
// file that is about to be rejected.
//
// Like REGISTER it runs on the launching machine with the wrapper's Python
// environment, is never cached (the vault can change between runs), and has no
// stub block, so `-stub` runs it for real. Its report is always shown: which
// genomes a run is about to skip is not debugging detail.
process PREFLIGHT {
    tag "${ready} request(s)"
    cache false
    debug true

    input:
    path resolved_csv
    val vault
    val ready

    output:
    path 'fetch_keys.txt', emit: keys

    script:
    def release_flag = isTrue(params.check_releases) ? '' : '--no-release-check'
    """
    export PATH="${params.venv_bin}:\$PATH"
    preflight.py "${resolved_csv}" "${vault}" --fetch-keys fetch_keys.txt ${release_flag}
    """
}

// Catalogue every download in the refgenie vault at params.outdir.
//
// ONE task for the whole request list, not one per genome, and that is deliberate:
//
//   * The vault has a single config file, so registrations cannot overlap.
//   * The order must be the request file's order. The first genome to claim a bare
//     assembly name keeps it (GRCm39 from Ensembl and from GENCODE cannot both be
//     'GRCm39'), so registering in whichever order the downloads happened to
//     finish would hand that name to a different genome from one run to the next.
//     bin/update_refgenie.py walks the resolved CSV top to bottom, as it always has.
//
// Each download arrives as a staged input, not by its published path under
// params.outdir. publishDir copies in the background and promises nothing about
// when it finishes, so reading the published copy from a downstream task is a race.
// `stageAs: 'sources/?/*'` gives each directory its own numbered parent, because
// two requests for one assembly produce two directories with the same name.
//
// Not cached. Registration is idempotent and cheap once a genome is in the vault,
// and a cached "already registered" would be wrong the moment the vault was edited
// or rebuilt outside this workflow.
//
// It runs on the machine that launched Nextflow, whatever executor the profile
// selects (see `withName: REGISTER` in nextflow.config), using the Python
// environment run_pipeline.sh builds: that is where refgenie is installed.
//
// There is no stub block, so `nextflow run -stub` runs this for real against the
// stubbed downloads - which is how tests/test_nextflow_stub.py exercises it.
process REGISTER {
    tag "${keys.size()} downloaded"
    maxForks 1
    cache false
    debug(params.debug.toString().toLowerCase() == 'true')

    input:
    path resolved_csv
    val vault
    tuple val(keys), path(sources, stageAs: 'sources/?/*')

    output:
    path 'registered.tsv', emit: manifest

    script:
    // keys[i] and sources[i] are the same download: both lists come from one list of
    // pairs, split in the workflow block. Each key is built there from values that
    // have been through safe(), and each source renders as its staged path.
    def dirs = sources instanceof List ? sources : [sources]
    def source_args = (0..<keys.size()).collect { i -> "--source '${keys[i]}=${dirs[i]}'" }.join(' \\\n        ')
    def cleanup_flag = isTrue(params.cleanup) ? '--cleanup' : ''
    """
    export PATH="${params.venv_bin}:\$PATH"
    update_refgenie.py "${resolved_csv}" "${vault}" \\
        --manifest registered.tsv ${cleanup_flag} \\
        ${source_args}
    """
}

// Derive a transcriptome from a genome and an annotation release already in the
// vault: `gffread -w` over the vault's own FASTA and GTF.
//
// This is the first process that does real computing rather than fetching or
// filing, so it is the first that goes wherever the profile sends it - a compute
// node under -profile sge. It needs no network and does not touch the vault's
// config: its inputs are staged from the vault's data/ tree and its outputs go to
// REGISTER_ASSET, which is the only thing that writes them back.
//
// Derived, never downloaded, and nothing is filtered. A provider's cDNA file
// leaves out non-coding genes and is cut from an assembly with more contigs than
// the vault's genome has; deriving gives the annotation's transcripts, all of
// them, from the sequence actually held. See recipes/fasta_txome.json.
//
// Three details that are easy to undo by accident:
//
//   * `fai` is an input although the script never names it. gffread looks for
//     <fasta>.fai beside the FASTA and, finding none, tries to write one - into
//     the vault. Staging both puts them side by side in the task directory.
//   * Sequence is upper-cased. The vault's genomes are soft-masked, so gffread
//     emits lower-case repeats; no provider's cDNA looks like that, and how each
//     quantifier treats them is not something a reference should leave open.
//   * An empty result is an error. gffread writes nothing, and exits 0, when the
//     GTF's contig names do not match the FASTA's - the one case where a GTF and
//     a genome that do not belong together would otherwise be filed as a pair.
//
// The provenance file has no timestamp on purpose: the same inputs and the same
// gffread give byte-identical outputs, and so the same asset digest.
process DERIVE_TXOME {
    tag "${genome}:${release}"
    label 'derive'
    cpus 1
    memory '4 GB'
    time '2h'

    input:
    tuple val(genome), val(annotation_asset), val(release), path(fasta), path(fai), path(gtf)

    output:
    tuple val(genome), val(annotation_asset), val(release), path("${genome}.transcripts.fa.gz"), path("${genome}.provenance.json")

    script:
    """
    if ! command -v gffread >/dev/null 2>&1; then
        echo "[${genome}] ERROR: gffread is not on PATH. The conda profile supplies it; outside conda, install it with: conda install -c bioconda gffread" >&2
        exit 3
    fi

    gzip -dcf "${gtf}" > annotation.gtf
    gffread -w transcripts.fa -g "${fasta}" annotation.gtf

    n=\$(grep -c '^>' transcripts.fa || true)
    if [ "\${n:-0}" -eq 0 ]; then
        echo "[${genome}] ERROR: gffread wrote no transcripts from ${annotation_asset}:${release}. The usual cause is contig names in the GTF that do not match the genome." >&2
        exit 1
    fi

    awk '/^>/ {print; next} {print toupper(\$0)}' transcripts.fa | gzip -nc > "${genome}.transcripts.fa.gz"

    version=\$(gffread --version 2>&1 | head -n 1)
    printf '{"asset": "fasta_txome", "genome": "%s", "annotation": "%s:%s", "tool": "gffread", "tool_version": "%s", "command": "gffread -w transcripts.fa -g GENOME.fa ANNOTATION.gtf", "sequence_case": "upper", "transcripts": %s}\\n' \\
        "${genome}" "${annotation_asset}" "${release}" "\$version" "\$n" > "${genome}.provenance.json"
    """

    // For `-stub`: one short record per transcript_id in the GTF, with no gffread.
    // Enough for the tests to follow a transcriptome from the manifest into the
    // vault; the real derivation is tested separately, where gffread is installed.
    stub:
    """
    gzip -dcf "${gtf}" | grep -o 'transcript_id "[^"]*"' | sort -u \\
        | sed 's/transcript_id "\\(.*\\)"/>\\1/' \\
        | while read -r header; do echo "\$header"; echo ACGTACGTACGTACGTACGT; done \\
        | gzip -nc > "${genome}.transcripts.fa.gz"
    n=\$(gzip -dc "${genome}.transcripts.fa.gz" | grep -c '^>' || true)
    printf '{"stub": true, "asset": "fasta_txome", "genome": "%s", "annotation": "%s:%s", "tool": "stub", "transcripts": %s}\\n' \\
        "${genome}" "${annotation_asset}" "${release}" "\$n" > "${genome}.provenance.json"
    """
}

// Fetch a third-party resource - ClinVar, the ENCODE blacklist - for a genome
// already in the vault. What can be fetched, and from where, is described by the
// files in resources/; see bin/resource_registry.py.
//
// bin/fetch_resource.py works out the published version before downloading and
// fetches only when the vault does not already hold it, so this task is cheap on
// every run but the first. It then verifies the file against the publisher's
// checksum where one is posted. `fetched/` always comes out: holding either the
// file and its fetch record, or a file called SKIP, which the workflow filters on.
//
// Like FETCH_GENOME it needs the network, carries the 'network' label, and is
// retried the same way. Its report is always shown.
//
// The stub writes a two-line file on the genome's first sequence, so that the
// real compatibility check which follows has something true to pass; it reports
// SKIP when the genome already holds the resource, as the real script would.
process FETCH_RESOURCE {
    tag "${name}:${variant} for ${genome}"
    label 'network'
    debug true

    errorStrategy {
        if( task.attempt > asInt(params.max_retries, 'max_retries') )
            return 'finish'
        sleep(asInt(params.retry_wait, 'retry_wait') * 1000L * task.attempt)
        return 'retry'
    }
    maxRetries(params.max_retries.toString().trim().isInteger() ? params.max_retries.toString().trim().toInteger() : 0)

    input:
    tuple val(genome), val(name), val(variant), val(held), path(chrom_sizes), path(registry)

    output:
    tuple val(genome), val(name), val(variant), path(chrom_sizes), path(registry), path('fetched')

    script:
    def release_flag = isTrue(params.check_releases) ? '' : '--no-release-check'
    """
    fetch_resource.py --registry "${registry}" --variant "${variant}" --genome "${genome}" \\
        --held "${held}" --outdir fetched ${release_flag}
    """

    stub:
    """
    mkdir -p fetched
    if [ -n "${held}" ]; then
        echo "stub: ${name}:${variant} is already held (${held})" > fetched/SKIP
        echo "  SKIP   ${name}:${variant} for ${genome}: held (${held})"
        exit 0
    fi
    contig=\$(head -n 1 "${chrom_sizes}" | cut -f 1)
    if grep -q '"kind": *"vcf"' "${registry}"; then
        file=${name}.vcf.gz
        printf '##fileformat=VCFv4.1\\n##source=stub\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n%s\\t5\\t.\\tA\\tC\\t.\\t.\\t.\\n%s\\t9\\t.\\tG\\tT\\t.\\t.\\t.\\n' "\$contig" "\$contig" \\
            | gzip -nc > "fetched/\$file"
    else
        file=${name}.bed.gz
        printf '%s\\t1\\t5\\tstub\\n' "\$contig" | gzip -nc > "fetched/\$file"
    fi
    printf '{"stub": true, "resource": "%s", "variant": "%s", "url": "stub", "version": "stub1", "file": "%s"}\\n' \\
        "${name}" "${variant}" "\$file" > "fetched/${name}.fetch.json"
    echo "  FETCH  ${name}:${variant} for ${genome} (stub)"
    """
}

// Decide whether a fetched resource belongs to the genome it was requested for,
// and index it if its registry entry says to. bin/prepare_resource.py explains
// the three checks; in short, the resource's sequence names and positions are
// compared with the genome's own chrom.sizes, and a resource that does not fit is
// refused here, loudly, and never reaches the vault.
//
// A refusal fails the task. 'finish' rather than the default, so the other
// resources and transcriptomes in the run are still completed and registered
// before the run exits non-zero.
//
// It runs on the launching machine with the wrapper's Python environment, whose
// pysam builds the tabix index. No stub block: under -stub this is the real check.
process PREPARE_RESOURCE {
    tag "${name}:${variant} for ${genome}"
    errorStrategy 'finish'
    debug true

    input:
    tuple val(genome), val(name), val(variant), path(chrom_sizes), path(registry), path(fetched)

    output:
    tuple val(genome), val(name), path('prepared')

    script:
    """
    export PATH="${params.venv_bin}:\$PATH"
    prepare_resource.py --registry "${registry}" --variant "${variant}" --genome "${genome}" \\
        --fetched "${fetched}" --chrom-sizes "${chrom_sizes}" --outdir prepared
    """
}

// File something in the vault against a genome already there: a transcriptome
// derived from it, or a resource fetched for it.
//
// The counterpart of REGISTER for everything that is not the download itself. One
// task per asset, one at a time: like REGISTER it edits the vault's single config,
// runs on the launching machine with the wrapper's Python environment, and is
// never cached. Order does not matter here - no two of these compete for a name -
// so there is no need to batch them into one task as REGISTER does.
//
// What to register arrives as data: the asset's name, which is also the name of
// its recipe in recipes/; the tag; the files, each with the key its recipe knows
// it by; and any parent assets to record. file_keys[i] names files[i].
process REGISTER_ASSET {
    tag "${genome}/${asset}:${asset_tag}"
    maxForks 1
    cache false
    debug(params.debug.toString().toLowerCase() == 'true')

    input:
    val vault
    tuple val(genome), val(asset), val(asset_tag), val(file_keys), path(files), val(parents)

    script:
    def staged = files instanceof List ? files : [files]
    def file_args = (0..<file_keys.size()).collect { i -> "--file '${file_keys[i]}=${staged[i]}'" }.join(' ')
    def parent_args = parents.collect { parent -> "--parent '${parent}'" }.join(' ')
    """
    export PATH="${params.venv_bin}:\$PATH"
    register_asset.py "${vault}" --genome "${genome}" --asset "${asset}" --tag "${asset_tag}" \\
        ${file_args} ${parent_args}
    """
}

workflow {
    // Where the vault goes, as an absolute path: REGISTER runs in its own work
    // directory, so a relative --outdir would otherwise point inside it.
    def vault = params.outdir.toString().startsWith('/') ? params.outdir.toString() : "${launchDir}/${params.outdir}".toString()

    // Fail now, not after hours of downloading, if there is nothing to register
    // with. run_pipeline.sh always builds this environment before it gets here.
    def do_register = isTrue(params.register)

    // Checked once here so a mistyped number stops the run with its own message,
    // before any task exists, and not from inside a directive later on.
    asInt(params.max_retries, 'max_retries')
    asInt(params.retry_wait, 'retry_wait')
    if( asInt(params.max_downloads, 'max_downloads') < 1 )
        error "--max_downloads must be at least 1"

    if( do_register && !file("${params.venv_bin}/python").exists() )
        error "No Python environment at '${params.venv_bin}', which the REGISTER step needs for refgenie. Run run_pipeline.sh once to create it, point --venv_bin at one, or pass --register false to download without cataloguing."

    // checkIfExists turns a wrong --input into an immediate, readable error. Without
    // it fromPath fails in two different unhelpful ways, both verified against
    // Nextflow 26.04:
    //
    //   * a literal path that does not exist is emitted anyway, so the channel holds
    //     one entry pointing at nothing and the run dies later inside splitCsv with a
    //     NoSuchFileException that names neither the parameter nor the cause;
    //   * a pattern that matches nothing yields an EMPTY channel, so no FETCH_GENOME
    //     task is created and the run exits 0 having downloaded not one genome.
    //
    // The second is the dangerous one: success is indistinguishable from doing
    // nothing. With the flag, both become "No files match pattern ..." before any
    // work starts.
    //
    // The smoke test cannot catch either: --dry-run exits before Nextflow is reached,
    // and `nextflow lint` is static analysis that never resolves a path. So the CI
    // guarding the rest of the request-file plumbing would sail straight past this,
    // which is what makes the flag load-bearing rather than decorative.
    channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        // row.release is null when the column is absent, which is the unpinned case
        // and must stay valid: older resolved CSVs predate the column.
        //
        // It is validated here as well as in the process script, because it now
        // forms part of the publish path, and publishDir is resolved independently
        // of the script block that carries the other check.
        //
        // The other three fields go through safe() here for the same reason - they
        // form the rest of that path - and so that a bad value stops the run before
        // any task is created, whether the script or the stub would have run.
        .map { row ->
            def rel = (row.release ?: '').toString().trim()
            if( rel && !(rel ==~ /^[0-9]+$/) )
                error "Unsafe or non-numeric release: '${rel}'"
            def url = safeUrl(row.url, 'url')
            def annotation_url = safeUrl(row.annotation_url, 'annotation_url')
            def provider = safe(row.provider, 'provider')
            // The address is the whole of such a request. Without it genomepy
            // would be asked to fetch a file called after the assembly.
            if( provider.toLowerCase() == 'url' && !url )
                error "Provider 'url' needs an address in the url column (assembly '${row.assembly}')"
            if( provider.toLowerCase() != 'url' && (url || annotation_url) )
                error "The url columns are only used with provider 'url' (assembly '${row.assembly}', provider '${provider}')"
            tuple(safe(row.species, 'species'), safe(row.assembly, 'assembly'),
                  provider, row.annotation, rel, url, annotation_url)
        }
        // Collected so the whole list can be checked for duplicates before any
        // task starts; flatMap then emits the rows one at a time again. The check
        // sits in its own map, not in the flatMap - see rejectDuplicates.
        .toList()
        .map { rows -> rejectDuplicates(rows) }
        .set { checked_ch }

    checked_ch.flatMap { rows -> rows }.set { requests_ch }

    // Which of the requests to download. With the preflight on, only those the
    // vault does not already satisfy; otherwise all of them. The preflight needs
    // the vault, so it is skipped along with registration.
    //
    // The key built here must match source_key() in bin/update_refgenie.py, which
    // is what the preflight writes; keep the two in step.
    if( do_register && isTrue(params.preflight) ) {
        PREFLIGHT(file(params.input, checkIfExists: true), vault,
                  checked_ch.map { rows -> rows.size() })

        // Wrapped in a list so `combine` appends the whole set of keys as ONE
        // element; unwrapped, it would spread them across the tuple.
        PREFLIGHT.out.keys
            .map { keys_file -> [ keys_file.readLines().findAll { line -> line } ] }
            .set { wanted_ch }

        requests_ch
            .combine(wanted_ch)
            .filter { _species, assembly, provider, _annotation, release, _url, _annotation_url, wanted ->
                wanted.contains("${provider.toLowerCase()}/${assembly}@${release}".toString())
            }
            .map { species, assembly, provider, annotation, release, url, annotation_url, _wanted ->
                tuple(species, assembly, provider, annotation, release, url, annotation_url)
            }
            .set { fetch_ch }
    }
    else {
        requests_ch.set { fetch_ch }
    }

    FETCH_GENOME(fetch_ch)

    if( do_register ) {
        // One (keys, directories) pair for the whole run, so REGISTER is a single
        // task that starts only when every download has finished. The key format
        // is source_key() in bin/update_refgenie.py; keep the two in step.
        //
        // An empty pair is passed on, not filtered out. When the preflight finds
        // every request already in the vault nothing is downloaded, and REGISTER
        // must still run: it is what applies an alias added to the request since
        // the genome was built, and what removes kept sources under --cleanup.
        FETCH_GENOME.out
            .map { provider, assembly, release, dir ->
                tuple("${provider.toLowerCase()}/${assembly}@${release}".toString(), dir)
            }
            .toList()
            .map { pairs -> tuple(pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] }) }
            .set { downloads_ch }

        REGISTER(file(params.input, checkIfExists: true), vault, downloads_ch)

        // Transcriptomes: part of every run unless `--txome false` turns them off.
        // REGISTER's manifest lists every requested
        // genome that is in the vault, with the annotation release the row means
        // and the releases a transcriptome already exists for. So this covers a
        // genome downloaded in this run and one that has been in the vault for
        // months alike, and derives only what is missing.
        //
        //   * a row with no annotation, or whose GTF is not in the vault, has
        //     nothing to derive from and is left alone;
        //   * two rows can mean the same genome and release (a pinned 116 and an
        //     unpinned row while 116 is current), hence `unique`;
        //   * the values come out of a file, so they go through safe() again
        //     before a script interpolates them.
        if( isTrue(params.txome) ) {
            REGISTER.out.manifest
                .splitCsv(header: true, sep: '\t')
                .filter { row ->
                    row.fasta && row.gtf && row.annotation_tag &&
                        !(row.fasta_txome_tags ?: '').tokenize(';').contains(row.annotation_tag)
                }
                .unique { row -> "${row.vault_name}:${row.annotation_tag}".toString() }
                .map { row ->
                    tuple(safe(row.vault_name, 'vault_name'),
                          safe(row.annotation_asset, 'annotation_asset'),
                          safe(row.annotation_tag, 'annotation_tag'),
                          file(row.fasta), file("${row.fasta}.fai"), file(row.gtf))
                }
                .set { txome_ch }

            DERIVE_TXOME(txome_ch)

            // The annotation release is both the tag and a recorded parent, so the
            // vault says which GTF a transcriptome came from.
            DERIVE_TXOME.out
                .map { genome, annotation_asset, release, transcripts, provenance ->
                    tuple(genome, 'fasta_txome', release, ['transcripts', 'provenance'],
                          [transcripts, provenance],
                          ["annotation=${annotation_asset}:${release}".toString()])
                }
                .set { txome_assets_ch }
        }
        else {
            channel.empty().set { txome_assets_ch }
        }

        // Third-party resources, for the rows that ask for any. One task per
        // genome and resource; `held` is the versions of it the genome already
        // has, which is what lets the fetch step skip a download.
        REGISTER.out.manifest
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.resources && row.chrom_sizes }
            .flatMap { row ->
                row.resources.tokenize(';').collect { request ->
                    def parts = request.tokenize(':')
                    def name = safe(parts[0], 'resource')
                    def variant = safe(parts[1], 'resource build')
                    def have = (row.assets ?: '').tokenize(';').find { item -> item.startsWith("${name}:") }
                    def held = have ? have.substring(name.size() + 1) : ''
                    if( held && !(held ==~ /^[A-Za-z0-9._,-]+$/) )
                        error "Unsafe characters in the versions held of '${name}': '${held}'"
                    tuple(safe(row.vault_name, 'vault_name'), name, variant, held,
                          file(row.chrom_sizes),
                          file("${projectDir}/resources/${name}.json", checkIfExists: true))
                }
            }
            .unique { item -> "${item[0]}/${item[1]}".toString() }
            .set { resource_ch }

        FETCH_RESOURCE(resource_ch)

        // Nothing to prepare where the fetch step found the vault up to date.
        FETCH_RESOURCE.out
            .filter { _genome, _name, _variant, _sizes, _registry, fetched -> !file("${fetched}/SKIP").exists() }
            .set { fetched_ch }

        PREPARE_RESOURCE(fetched_ch)

        // prepared/ holds the data file, its index if it has one, provenance.json
        // and TAG. Each file is handed over under the key its recipe expects.
        PREPARE_RESOURCE.out
            .map { genome, name, prepared ->
                def version = safe(file("${prepared}/TAG").text.trim(), 'resource version')
                def files = prepared.listFiles().findAll { item -> item.name != 'TAG' }.sort { item -> item.name }
                def keys = files.collect { item ->
                    item.name == 'provenance.json' ? 'provenance' : (item.name.endsWith('.tbi') ? 'tbi' : name)
                }
                tuple(genome, name, version, keys, files, [])
            }
            .set { resource_assets_ch }

        REGISTER_ASSET(vault, txome_assets_ch.mix(resource_assets_ch))
    }
}
