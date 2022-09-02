library(Seurat)
library(viridis)
library(tidyverse)

hca_path = "~/research_project/markers/hca_human_heart.rds"

# Load HCA annotated object 
hca_obj <- readRDS(hca_path)

print("Loaded HCA object")

# Display metadata fields
print(names(hca_obj@meta.data))

# Object lacks RNA assay key which prevents sub-setting
hca_obj@assays$RNA@key <- "rna_"

# Select LV Nuclei used in the Atlas and filter for genes present in at least 
# 10 cells and cells with at least 300 genes
hca_obj <- hca_obj[rowSums(hca_obj[["RNA"]]@counts >0) > 10, 
                   hca_obj$source == "Nuclei" &
                     hca_obj$tissue == "heart left ventricle" &
                     hca_obj$Used == "Yes" &
                     hca_obj$cell_states != "doublets" &
                     hca_obj$nFeaturess_RNA > 300]

print("Successfully subset object")

gc()

# Set identity class for the object to cell_type
Idents(hca_obj) <- "cell_type"
# Find DE gene markers for all identity classes in hca_obj
print("Finding Markers")
ct_markers <- FindAllMarkers(hca_obj, assay = "RNA", 
                             verbose = T, test.use = "wilcox")

write.table(ct_markers, 
            file = "~/research_project/markers/hca_ct_markers.csv",
            sep = ",")



