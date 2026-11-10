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
def rejectDuplicates(List rows) {
    def seen = [:]
    rows.each { row ->
        def provider = (row[2] ?: '').toString().toLowerCase()
        def key = "${provider}/${row[1]}@${row[4] ?: 'current'}".toString()
        if( seen.containsKey(key) )
            error "Duplicate request for ${row[1]} from '${row[2]}' (release: ${row[4] ?: 'current'}). Each provider, assembly and release may appear once; merge the rows."
        seen[key] = true
    }
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
        .map { row ->
            def rel = (row.release ?: '').toString().trim()
            if( rel && !(rel ==~ /^[0-9]+$/) )
                error "Unsafe or non-numeric release: '${rel}'"
            tuple(row.species, row.assembly, row.provider, row.annotation, rel)
        }
        // Collected so the whole list can be checked for duplicates before any
        // task starts; flatMap then emits the rows one at a time again.
        .toList()
        .flatMap { rows -> rejectDuplicates(rows) }
        .set { requests_ch }

    FETCH_GENOME(requests_ch)
}
