data_analysis_scrnaseq.R
###----------------------
##         Setup
###----------------------

#Required packages
packages = c("Seurat","ggplot2","cowplot")

#Function to install packages
install.load.pkg <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}

#Install and load needed packages if not installed
install.load.pkg(packages)


#Import functions
source("./bin/functions_data_analysis.r", chdir = TRUE)

#Get command-line arguments for importing control/treatment count matrices
args = commandArgs(trailingOnly = TRUE)

###----------------------------------------------------------------
###  Check whether to process case-control data or one-sample data
###----------------------------------------------------------------

if (length(args)==2) {
#Import count matrices containing Genes (Rows) by Sample (Columns) from cellranger count analysis
control.data = read.table(args[1],sep = "\t", header=TRUE )
treatment.data = read.table(args[2], sep= "\t", header=TRUE )

# Define variables:
condition = c("EXP_CONTROL","EXP_TREATMENT")
variable_value = c("CONTROL","TREATMENT")

#Set up Seurat object - control
control = create_seurat_object(control.data,
                                condition[1],
                                variable_value[1])
                            

#Set up Seurat object - perturbation treatment
treatment = create_seurat_object(treatment.data,
                                 condition[2],
                                 variable_value[2])



# --------------------
# INTEGRATION OF CONTROL AND TREATMENT DATA
# --------------------

anchors.data = FindIntegrationAnchors(object.list = list(control,treatment), dims = 1:20)
combined.objects = IntegrateData(anchorset = anchors.data , dims = 1:20)

# --------------------
# INTEGRATED ANALYSIS
# --------------------

DefaultAssay(combined.objects) = "integrated"

# --------------------
# CLUSTERING AND VISUALIZATION
# --------------------

#Clustering
combined.objects = ScaleData(combined.objects,verbose = FALSE)
combined.objects = RunPCA(combined.objects,npcs=30,verbose=FALSE)

#tSNE
combined.objects = RunUMAP(combined.objects,reduction = "pca",dims = 1:20)
combined.objects = FindNeighbors(combined.objects,reduction="pca",dims = 1:20)
combined.objects = FindClusters(combined.objects,resolution = 0.5)

#Visualization

p1 = DimPlot(combined.objects, reduction = "umap", group.by = "variable_value")
p2 = DimPlot(combined.objects, reduction = "umap", label = TRUE)
p3 = plot_grid(p1, p2)

#Split by experiment
p4 = DimPlot(combined.objects,reduction = "umap", split.by = "variable_value")


# Create an output directory (if it doesn't exist)
output_dir <- "results/stats_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
###-------------------------------
#Export plots to output directory 
###-------------------------------

list.plots = list(p1,p2,p3,combined.plot)

for(i in seq_along(list.plots)){
    output_file = file.path(output_dir,paste0("plot_",i,".pdf"))
    ggsave(output_file, list.plots[[i]],width=10,height=5)
    cat("Saved plot to",output_file, "\n")
}
} else if (length(args)==1) { #### If only one sample is provided

#Import count matrices containing Genes (Rows) by Sample (Columns) from cellranger count analysis
input.data = read.table(args,sep = "\t", header=TRUE )

#Set up Seurat object 
input.seurat = create_seurat_object_onesample(input.data)
                            

# --------------------
# RNA ANALYSIS
# --------------------

DefaultAssay(input.seurat) = "RNA"

# --------------------
# CLUSTERING AND VISUALIZATION
# --------------------

#Clustering
input.seurat = ScaleData(input.seurat,verbose = TRUE)
input.seurat = RunPCA(input.seurat,npcs=30,verbose=FALSE)

#tSNE
input.seurat = RunUMAP(input.seurat,reduction = "pca",dims = 1:20)
input.seurat = FindNeighbors(input.seurat,reduction="pca",dims = 1:20)
input.seurat = FindClusters(input.seurat,resolution = 0.5)

#Visualization
p1 = DimPlot(input.seurat, reduction = "umap", label = TRUE)
p1

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

} else if (args > 2) {
   stop("Error: Only one-sample and pairwise control-treatment comparisons are supported by this analysis")
}


