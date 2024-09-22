#!/usr/bin/env Rscript

process SEURAT {

    input:
        tuple path(matrix_control), path(matrix_treatment)
    
    output:
        path outdir

    script:
    """
    R ../bin/data_analysis.r matrix_control matrix_treatment

    """

}