#!/usr/bin/env nextflow

/*-----------------------------
// DEFINE DEFAULT PARAMETERS
-----------------------------*/

params.projectDir = "$PWD/10x_scrna"
params.transcriptome = "${params.projectDir}/data/ref/"
params.reads = "${params.projectDir}/data/reads/_*R{1,2}_.fastq.gz"
params.outdir = "${params.projectDir}/results"
params.fastqc = "${params.outdir}/fastqc"
params.cellranger_count = "${params.outdir}/cellranger_count"
/*------------------
// IMPORT MODULES
-------------------*/

include { FASTQC             } from './modules/fastqc.nf'
include { MULTIQC            } from './modules/multiqc.nf'
include { CELLRANGER_COUNT   } from './modules/cellranger_count.nf'
include { SEURAT             } from './modules/seurat.nf'

// Workflow definition
workflow {

    //Create a channel with pairs of reads for each sample
    reads_ch = Channel.fromFilePairs("${params.reads}")

    //Extract sample IDs from the reads channel
    sample_ids = reads_ch.map {tuple -> tuple[0].basename()}

    //Run FASTQC on each pair of reads
    fastqc_raw = FASTQC(reads_ch)

    //Generate a multiqc report for all samples
    multiqc = MULTIQC(fastqc_raw)

    //Generate count matrix with CellRanger
    cellranger = CELLRANGER_COUNT(reads_ch)

    //Clustering analysis and visualization 
    seurat = SEURAT(cellranger)
}