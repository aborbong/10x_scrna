#!/usr/bin/env nextflow

/*
 * Create a multiqc report for all samples
 */

process MULTIQC {

    input:
        path fastqc
    
    output:
        path '${params.outdir}/multiqc_final.html'

    script:
    """
        multiqc ${params.outdir}/fastqc -o ${params.outdir}
    """
}