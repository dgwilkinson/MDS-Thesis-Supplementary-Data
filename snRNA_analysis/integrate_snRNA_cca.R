library(Seurat)
library(harmony)
library(cluster)
library(optparse)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)
library(future)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

# Enable multiprocessing
# plan("multisession", workers = 16)
# # Allow up to 8 GB per worker (core)

# Allow up to 64 GiB for SCT Prep
options(future.globals.maxSize = 64 * 1024^3)

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
objects <- list.files(objects_dir, pattern = "^[^_]*.rds", full.names = T)
names(objects) <- str_extract(objects, "[CKL]{2}[0-9]+")

# Remove KL objects lacking metadata info and samples that fail qc
exclusion <- c("CK161", "CK369")
obj_filter <- grepl("KL", objects) | names(objects) %in% exclusion
objects <- objects[!obj_filter]

# Select only samples from the first four patients
selected_samples <- read_xlsx(metadata_path, skip = 1) %>%
  select(c("sample_id", "patient", "patient_group")) %>%
  filter(patient %in% c("P1", "P2", "P3", "P4")) %>%
  pull(sample_id)

# Load in valid objects
loaded_objects <- map(objects[names(objects) %in% selected_samples], readRDS)

print("Successfully loaded all objects")

# Overwrite pca reduction with SCT_PCA to work with rpca integration
# loaded_objects <- lapply(loaded_objects, function(obj){
#   obj@reductions$pca <- obj@reductions$SCT_pca
#   obj@reductions$SCT_pca <- NULL
#   return(obj)
# })

gc()

# Find integration anchors for the loaded objects
integration_features <- SelectIntegrationFeatures(loaded_objects, 
                                                  nfeatures = 3000)

loaded_objects <- PrepSCTIntegration(object.list = loaded_objects,
                   assay = "SCT",
                   anchor.features = integration_features,
                   verbose = F)

gc()

# Uncomment reference to use Control "CK158" as reference
integration_anchors <- FindIntegrationAnchors(object.list = loaded_objects,
                                              #reference = c(1),
                                              normalization.method = "SCT",
                                              reduction = "cca",
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
  labs(title = "P1-4 CCA UMAP",
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

saveRDS(integrated_obj, file = paste(out_path, "cca_integrated_p14.rds", sep = "/"))

pdf(file = paste(fig_path, "integrated_snRNA_cca_plots_p14.pdf", sep = "/"),
    width = 10, height = 10)
original_umap
dev.off()




