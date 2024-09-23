#!/usr/bin/env nextflow
/*
 * Check if the input reads are in bcl format, and converts them into fastq if necessary
   //https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-fq
*/

 process BCL_TO_FASTQ {

    input: 
        path reads
        
    output:
        path "*.fastq"

    script:
    """
    if [[ \$(basename ${reads}) == *.bcl || \$(basename ${reads}) == *.BCL ]]; then

        cellranger mkfastq --run=${reads} 

     

    fi
    """
} 
