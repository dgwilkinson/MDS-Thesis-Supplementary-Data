library(Seurat)
library(zellkonverter)
library(optparse)
library(tidyverse)

options <- list(
  make_option(c("--seurat_path"), action = "store", type = "character"),
  make_option(c("--seurat_name"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

seurat_obj <- readRDS(paste0(seurat_path, seurat_name))

# Select the RNA or Spatial assay
assay <- grep("RNA|Spatial", names(seurat_obj), value = T)

# h5ad only requires RNA as c2l acts on raw counts and will account for technical variation in NB regression
writeH5AD(as.SingleCellExperiment(seurat_obj, assay = assay),
          file = paste0(seurat_path, gsub(".rds", ".h5ad", seurat_name)))

