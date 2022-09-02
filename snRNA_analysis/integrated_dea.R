library(Seurat)
library(optparse)
library(viridis)
library(tidyverse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--marker_path"), action = "store", type = "character"))

sample_path = "~/research_project/dataset/snRNA/integrated_object/cca_integrated_obj_DS.rds"

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load integrated object and set identity class to "opt_clust"
integrated_obj <- readRDS(sample_path)
Idents(integrated_obj) <- "opt_clust"

# Set default assay to SCT data
DefaultAssay(integrated_obj) <- "SCT"

# Find  definitional marker genes for each cluster
integrated_obj <- PrepSCTFindMarkers(integrated_obj, 
                   assay = "SCT",
                   verbose = T)

cluster_markers <- FindAllMarkers(integrated_obj,
                                  test.use = "wilcox",
                                  assay = DefaultAssay(integrated_obj),
                                  slot = "data")

write.table(cluster_markers, 
            file = paste0(marker_path, sample_id, "_cluster_markers.csv"),
            sep = ",")



