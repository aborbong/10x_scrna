#!/usr/bin/env nextflow

https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-aggr

process CELLRANGER_MULTI {

    input:
        tuple path(reads), path(libraries)
    
    output:
        path "${params.outdir}/cellranger_count/${sample_id}"//"cellranger_output_${sample_id}"

    script:
    """
    mkdir -p ${params.outdir}/cellranger_count/${sample_id}

        cellranger multi --id=aggre
                        --csv=libraries
                         

    """

}