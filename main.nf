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
    def annotation_flag = annotation.toString().toLowerCase() == 'true' ? '--annotation' : ''
    
    """
    python3 ${baseDir}/bin/download.py \\
        --species "${species}" \\
        --assembly "${assembly}" \\
        --provider "${provider}" \\
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
