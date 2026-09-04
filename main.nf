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
    publishDir { "${params.outdir}/${provider}/${species}" }, mode: 'copy'

    input:
    tuple val(species), val(assembly), val(provider), val(annotation)

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

    // download.py is called without an explicit path - Nextflow automatically adds
    // the project bin/ directory to PATH inside every process (fixes issue 3.1).
    """
    download.py \\
        --species "${sp}" \\
        --assembly "${asm}" \\
        --provider "${prv}" \\
        ${annotation_flag}
    """
}

workflow {
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.species, row.assembly, row.provider, row.annotation) }
        .set { requests_ch }

    FETCH_GENOME(requests_ch)
}
