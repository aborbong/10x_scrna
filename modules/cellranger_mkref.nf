#!/usr/bin/env nextflow

process CELLRANGER_MKREF {

    input:
         path genome
    output:
        path genome
    script:
    """
        mkdir -p ${outdir}/cellranger_count

        cd ${outdir}/cellranger_count

        cellranger mkref \
                --nthreads={CPUs} \
                --genome=${genome}
                --fasta=${genome}
                --genes=${genome.fq.gz}.gtf


    """

}

