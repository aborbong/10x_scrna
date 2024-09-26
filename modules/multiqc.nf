#!/usr/bin/env nextflow

/*
 * Create a multiqc report for all samples
 */

process MULTIQC {

    input:
        path fastqc 
    
    output:
        path "${params.outdir}/multiqc_report.html" 

    script:
    """
        multiqc .
    """
}