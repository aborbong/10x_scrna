#!/usr/bin/env nextflow
/*
 * Check quality with FASTQC
 */

 process FASTQC {

    input:
        path reads

    output:
        path fastqc
        

    """
    mkdir -p ${params.outdir}/fastqc

     fastqc \
        ${reads} \
        --threads 4 \
        --memory 4G \
        --outdir ${params.outdir}/fastqc \
     

   """
}
