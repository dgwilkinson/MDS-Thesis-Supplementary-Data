library(Seurat)
library(tidyverse)
library(optparse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load in Processed Object
seurat_obj <- readRDS(sample_path)

# Downsample clusters to 10%
valid_barcodes <- seurat_obj@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  group_by("opt_clust") %>%
  slice_sample(prop = 0.1) %>%
  unnest(cols = everything()) %>%
  pull("barcodes")

# Subset into downsampled object
DS_obj <- seurat_obj[,valid_barcodes]

# Save downsampled object
saveRDS(DS_obj, file=gsub(".rds", "_DS.rds", sample_path))