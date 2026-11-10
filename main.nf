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
    debug params.debug

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
        if( task.exitStatus == 3 || task.attempt > params.max_retries )
            return 'finish'
        sleep((params.retry_wait as long) * 1000 * task.attempt)
        return 'retry'
    }
    maxRetries params.max_retries

    // Every task is a download from one of four public providers. Unbounded, a long
    // request list opens one connection per row at once, which is unkind to them and
    // tends to get throttled.
    maxForks params.max_downloads

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
    publishDir {
        release ? "${params.outdir}/${provider}/${species}/release-${release}" : "${params.outdir}/${provider}/${species}"
    }, mode: 'copy'

    input:
    tuple val(species), val(assembly), val(provider), val(annotation), val(release)

    output:
    path "${assembly}"

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
    printf '%s\\n' '{"stub": true, "species": "${sp}", "assembly": "${asm}", "provider": "${prv}", "annotation_requested": ${want_gtf}${release_json}}' \\
        > "${asm}/provenance.json"
    """
}

workflow {
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
}
