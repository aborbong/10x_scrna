
#Function creates Seurat objects from count matrices for control and case samples
#Input:
#Output: Seurat objects
create_seurat_object <- function(counts_data, 
                                    condition, 
                                    variable_value, 
                                    min_cells = 5, 
                                    min_features = 500, 
                                    n_var_features = 2000) {
  
  # Create Seurat object
  seurat_obj = CreateSeuratObject(counts = counts_data, 
                                    project = condition, 
                                    min.cells = min_cells)
  
  # Add stim label
  seurat_obj$variable_value = variable_value
  
  # Subset based on number of features
  seurat_obj = subset(seurat_obj, subset = nFeature_RNA > min_features)
  
  # Normalize data
  seurat_obj = NormalizeData(seurat_obj, verbose = FALSE)
  
  # Find variable features
  seurat_obj = FindVariableFeatures(seurat_obj, 
                                    selection.method = "vst", 
                                    nfeatures = n_var_features)
  
  return(seurat_obj)
}

create_seurat_object_onesample <- function(counts_data, 
                                    min_cells = 5, 
                                    min_features = 500, 
                                    n_var_features = 2000) {
  
  # Create Seurat object
  seurat_obj = CreateSeuratObject(counts = counts_data, 
                                    min.cells = min_cells)
  
  # Subset based on number of features
  seurat_obj = subset(seurat_obj, subset = nFeature_RNA > min_features)
  
  # Normalize data
  seurat_obj = NormalizeData(seurat_obj, verbose = FALSE)
  
  # Find variable features
  seurat_obj = FindVariableFeatures(seurat_obj, 
                                    selection.method = "vst", 
                                    nfeatures = n_var_features)
  
  return(seurat_obj)
}