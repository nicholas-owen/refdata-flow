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
// annotation, release. Provider is compared case-insensitively because genomepy
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
    tuple val(species), val(assembly), val(provider), val(annotation), val(release)

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

    // download.py is called without an explicit path - Nextflow automatically adds
    // the project bin/ directory to PATH inside every process (fixes issue 3.1).
    """
    download.py \\
        --species "${sp}" \\
        --assembly "${asm}" \\
        --provider "${prv}" \\
        ${annotation_flag} \\
        ${release_flag}
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
    def rel_tag = release ?: 'stub'
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
    tag "${keys.size()} genome(s)"
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
            tuple(safe(row.species, 'species'), safe(row.assembly, 'assembly'),
                  safe(row.provider, 'provider'), row.annotation, rel)
        }
        // Collected so the whole list can be checked for duplicates before any
        // task starts; flatMap then emits the rows one at a time again. The check
        // sits in its own map, not in the flatMap - see rejectDuplicates.
        .toList()
        .map { rows -> rejectDuplicates(rows) }
        .flatMap { rows -> rows }
        .set { requests_ch }

    FETCH_GENOME(requests_ch)

    if( do_register ) {
        // One (keys, directories) pair for the whole run, so REGISTER is a single
        // task that starts only when every download has finished. The key format
        // is source_key() in bin/update_refgenie.py; keep the two in step.
        //
        // The filter drops the empty case - a request file with a header and no
        // rows - where there is nothing to stage and nothing to register.
        FETCH_GENOME.out
            .map { provider, assembly, release, dir ->
                tuple("${provider.toLowerCase()}/${assembly}@${release}".toString(), dir)
            }
            .toList()
            .filter { pairs -> pairs.size() > 0 }
            .map { pairs -> tuple(pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] }) }
            .set { downloads_ch }

        REGISTER(file(params.input, checkIfExists: true), vault, downloads_ch)
    }
}
