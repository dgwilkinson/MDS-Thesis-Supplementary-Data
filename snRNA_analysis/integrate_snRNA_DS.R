library(Seurat)
library(cluster)
library(optparse)
library(tidyverse)
library(viridis)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

options(future.globals.maxSize = 8000 * 1024^2)

options <- list(
  make_option(c("--objects_dir"), action = "store", type = "character"),
  make_option(c("--out_path"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--metadata_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

#Locate objects
objects <- list.files(objects_dir, pattern = "^.*_DS.rds", full.names = T)
names(objects) <- str_extract(objects, "[CKL]{2}[0-9]+")

# Remove KL objects lacking metadata info and samples that fail qc
exclusion <- c("CK161", "CK369")
obj_filter <- grepl("KL", objects) | names(objects) %in% exclusion
objects <- objects[!obj_filter]

# Load in valid objects
loaded_objects <- map(objects, readRDS)

print("Successfully loaded all objects")

# Overwrite pca reduction with SCT_PCA to work with rpca integration
loaded_objects <- lapply(loaded_objects, function(obj){
  obj@reductions$pca <- obj@reductions$SCT_pca
  obj@reductions$SCT_pca <- NULL
  return(obj)
})

gc()

# Find integration anchors for the loaded objects

integration_features <- SelectIntegrationFeatures(loaded_objects, 
                                                  nfeatures = 3000)

loaded_objects <- PrepSCTIntegration(object.list = loaded_objects,
                   assay = "SCT",
                   anchor.features = integration_features,
                   verbose = F)

gc()

integration_anchors <- FindIntegrationAnchors(object.list = loaded_objects,
                                              normalization.method = "SCT",
                                              reduction = "rpca",
                                              anchor.features = integration_features)

print("Found Integration Anchors")

integrated_obj <- IntegrateData(integration_anchors,
                                new.assay.name = "integrated",
                                normalization.method = "SCT")

rm(loaded_objects)
rm(integration_anchors)

print("Run Dimensionality Reduction")

integrated_obj <- integrated_obj %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F,
         npcs = 30)  %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          verbose = F, reduction.name = "umap_original")

print("Plotting DR")

original_umap <- integrated_obj %>% DimPlot(., reduction = "umap_original",
                                               group.by = "orig.ident",
                                               pt.size = .1,
                                               cols = viridis(length(unique(.$orig.ident)))) +
  labs(title = "Downsampled RPCA UMAP",
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

saveRDS(integrated_obj, file = paste(out_path, "rpca_integrated_obj_DS.rds"))

pdf(file = paste(fig_path, "integrated_snRNA_rpca_plots_DS.pdf", sep = "/"),
    width = 10, height = 10)
original_umap
dev.off()




