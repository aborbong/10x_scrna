# 10X Genomics scRNA-seq Analysis Pipeline

This pipeline performs primary and secondary analysis of 10X Genomics single-cell RNA-seq (scRNA-seq data). It includes the following steps:

#### 1. Primary analysis (Cell Ranger):
Preprocessing of raw FASTQ files generated by the 10X Genomics platform
Read quality assessment using FastQC and quality report generated with multiqc
Alignment of reads to a reference genome 
Gene quantification and UMI counting
Generation of count matrices for downstream analysis*

Note: Supports the analysis of multiple samples (running `cellranger multi`)

#### 2. Secondary analysis (Cell Ranger & Seurat) 
Data normalization
Dimensionality reduction & cell clustering
Differential gene expression analysis
Visualization

## Installation

Requirements:

Before running the pipeline, make sure that the following dependencies and tools are installed and included in your path:

1. Nextflow: See more in https://www.nextflow.io/docs/latest/overview.html
2. Cell ranger: See more in https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
3. FastQC: See more in 
4. MultiQC: See more in https://multiqc.info/docs/getting_started/installation/ 

### Installing dependencies

#### 1. Install Nextflow

1. Install Nextflow
```
`curl -s https://get.nextflow.io | bash`
```
2. Make Nextflow executable
```
`chmod +x nextflow`
```
3. Move Nextflow into an executable path:
```
`sudo mv nextflow /usr/local/bin`     
 ```
4. Confirm that Nextflow is installed correctly" 
```
`nextflow info`
```
#### 2. Install Cell Ranger

1. Download and unpack the cellranger `cellranger-x.y.z.tar.gz/.tar.xz` compressed folder from the 10xgenomics website: https://www.10xgenomics.com/support/software/cell-ranger/downloads

2. Add Cell Ranger to your PATH
```
`export PATH=/file-location:$PATH`
```

#### 3. Install FastQC 
Using conda:
```bash
`conda install bioconda::fastqc`
```

Using apt install:

```bash
`sudo apt-get install fastqc`
#or download from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
```

#### 4. Install MultiQC 
Using conda:

```bash
`conda install multiqc`
```

Using apt install:
```bash
`sudo apt-get install multiqc`
```

### Alternative: Create a conda environment to manage dependencies:

```bash
conda create -n scrna-seq-pipeline -c bioconda nextflow fastqc multiqc trimmomatic
conda activate scrna-seq-pipeline
```
