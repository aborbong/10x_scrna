data_analysis_scrnaseq.R

##Setup

#Install required packages
install.packages("Seurat")
install.packages("ggplot2")
install.packages("cowplot")

#Load libraries
library(Seurat)
library(ggplot2)
library(cowplot)

#Import functions
source("./bin/functions_data_analysis.r", chdir = TRUE)

#Get command-line arguments for importing control/treatment count matrices
args = commandArgs(trailingOnly = TRUE)

#Check if the right number of arguments is provided
if(length(args) > 1) {
    stop("This function only supports one sample. Two comape control-treatment datasets, \n
          please use --mode case_control in the nextflow.config file")
}

#Import count matrices containing Genes (Rows) by Sample (Columns) from cellranger count analysis
input.data = read.table(file = args[1],sep = "\t", header=TRUE )

#Set up Seurat object 
input.seurat = create_seurat_object_onesample(input.data)
                            

# --------------------
# INTEGRATED ANALYSIS
# --------------------

DefaultAssay(input.seurat) = "integrated"

# --------------------
# CLUSTERING AND VISUALIZATION
# --------------------

#Clustering
clustering.object = ScaleData(input.seurat,verbose = FALSE)
input.seurat = RunPCA(input.seurat,npcs=30,verbose=FALSE)

#tSNE
input.seurat = RunUMAP(input.seurat,reduction = "pca",dims = 1:20)
input.seurat = FindNeighbors(input.seurat,reduction="pca",dims = 1:20)
input.seurat = FindClusters(input.seurat,resolution = 0.5)

#Visualization
p1 = DimPlot(input.seurat, reduction = "umap", label = TRUE)


# Create an output directory (if it doesn't exist)
output_dir <- "results/stats_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

#Export plots to output directory 

list.plots = list(p1)

for(i in seq_along(list.plots)){
    output_file = file.path(output_dir,paste0("plot_",i,".pdf"))
    ggsave(output_file, list.plots[[i]],width=10,height=5)
    cat("Saved plot to",output_file, "\n")
}


# -----------------------------------------------------
# DIFFERENTIAL ABUNDANCE ANALYSIS OF CELL TYPES/GENES
# -----------------------------------------------------

DefaultAssay(combined.objects)<-"RNA"

markers = FindConservedMarkers(combined.objects,
                                ident.1 = 7,
                                grouping.var = "variable_value",
                                verbose = FALSE)




