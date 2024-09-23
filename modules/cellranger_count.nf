#!/usr/bin/env nextflow

process CELLRANGER_COUNT {

    input:
        tuple path(reads),path(libraries),path(genome),val(sample_id)
    output:
        path cellranger_count

    script:
    """
        mkdir -p ${outdir}/cellranger_count

        cd ${outdir}/cellranger_count

        cellranger count --id=${sample_id} \
                         --transcriptome=${genome} \
                         --sample=${sample_id} \
                         --create-bam=true \
                         --fastqs ${reads} \
                         --localcores=8 \
                         --localmem=24

    """

}

