#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FETCH_GENOME {
    tag "${species}_${assembly}"

    publishDir "${params.outdir}/${provider}/${species}", mode: 'copy'

    input:
    tuple val(species), val(assembly), val(provider), val(annotation)

    output:
    path "${assembly}"

    script:
    // Validate inputs against a strict allowlist before any shell interpolation.
    // Only alphanumeric characters, dots, underscores, and hyphens are permitted.
    // Nextflow's error() terminates the pipeline immediately with a clear message
    // if an unexpected value is encountered, preventing shell injection.
    def safe = { String val, String field ->
        if (!(val ==~ /^[A-Za-z0-9._-]+$/))
            error "Unsafe characters detected in ${field}: '${val}'"
        return val
    }
    def sp  = safe(species,  'species')
    def asm = safe(assembly, 'assembly')
    def prv = safe(provider, 'provider')
    // annotation_flag is safe: derived from a ternary with two hardcoded strings,
    // no user-supplied data is interpolated into it.
    def annotation_flag = annotation.toString().toLowerCase() == 'true' ? '--annotation' : ''

    // download.py is called without an explicit path — Nextflow automatically adds
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
