library(Seurat)
library(tidyverse)
library(optparse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load in Processed Object
seurat_obj <- readRDS(sample_path)

# Perform SCTransform and append normalised result to "SCT" assay
seurat_obj <- SCTransform(seurat_obj,
                          assay = "RNA",
                          verbose = F)

# Verify addition of "SCT" assay
names(seurat_obj)

# Perform PCA and UMAP, scaling data not advisable for SCT workflow
seurat_obj <- seurat_obj %>%
  RunPCA(verbose = F,
         reduction.name = "SCT_pca",
         npcs = 30) %>%
  RunUMAP(reduction = "SCT_pca",
          dims = 1:30, verbose = F,
          reduction.name = "SCT_umap")

sct_umap <- DimPlot(seurat_obj, group.by = "opt_clust",
                    reduction = "SCT_umap")+
  labs(title = paste(sample_id, "SCT UMAP"),
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

saveRDS(seurat_obj, file=sample_path)

pdf(file = paste0(fig_path, sample_id, "_SCT_UMAP.pdf"),
    width = 10, height = 10)
sct_umap
dev.off()




