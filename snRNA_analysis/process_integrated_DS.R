library(Seurat)
library(cluster)
library(optparse)
library(tidyverse)
library(viridis)
library(readxl)
library(harmony)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--metadata_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load integrated object
integrated_obj <- readRDS(sample_path)

# Read metadata for patient, sample and batch information
batch_metadata <- read_xlsx(metadata_path, skip = 1) %>%
  select(c("sample_id", "patient", "patient_group"))

# Add batch metadata to the integrated object
print("Adding metadata")
integrated_obj@meta.data <- integrated_obj@meta.data %>%
  rownames_to_column() %>%
  left_join(batch_metadata, by = c("orig.ident" = "sample_id")) %>%
  column_to_rownames()


integrated_PCA <- integrated_obj %>% DimPlot(., reduction = "pca",
                                              group.by = "orig.ident",
                                              pt.size = .1,
                                              cols = viridis(length(unique(.$orig.ident)))) +
  labs(title = paste(sample_id, "PCA"),
       x = expression(1^st~" Principal Component"),
       y = expression(2^nd~" Principal Component"))

# Plot correlation of PCs with confounding variables 
conf_variables <- c("orig.ident", "patient", "patient_group")

corr_conf <- function(integrated_obj, assay = "PCA"){
  library(scater)
  
  print("Plotting confounding variables for PCA")
  
  conf_plot <- as.SingleCellExperiment(integrated_obj) %>% 
    getExplanatoryPCs(variables = conf_variables,
                      dimred = assay, n_dimred =30) %>%
    as.data.frame() %>%
    rownames_to_column("components") %>%
    pivot_longer(-c(components), names_to = "variable") %>%
    ggplot(mapping = aes(x = factor(components, unique(components)),
                         y = log(value,base = 10), group = variable, col = variable)) +
    geom_line(stat = "identity")+
    geom_point(shape = 1, size = 2) + 
    scale_color_viridis(discrete=T, labels = c("Sample ID", "Patient", "Batch")) + 
    theme_classic() + 
    scale_y_continuous(breaks = seq(-3, 2),
                       labels = sprintf("%g", 10**seq(-3, 2))) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    labs(x= "Principal Component", y = "Explained Variance (%)", color = "Confounding Variables")
  
  detach("package:scater", unload = T)
  
  return(conf_plot)
}

conf_plot_original <- corr_conf(integrated_obj)

# Run Harmony Batch effect correction
print("Running Harmony")

integrated_obj <- RunHarmony(integrated_obj,
                                group.by.vars = conf_variables,
                                reduction = "pca",
                                assay.use = DefaultAssay(integrated_obj),
                                plot_convergence = T,
                                max.iter.harmony = 20)

harmony_pca <- DimPlot(integrated_obj, 
                       reduction = "harmony", 
                       group.by = "orig.ident",
                       pt.size = .1) +
  labs(title = "Harmony PCA by Sample", 
       x = expression(1^st~" Principal Component"),
       y = expression(2^nd~" Principal Component")) +
  scale_color_viridis(discrete = T)

# Run UMAP on the pca harmony cell embedding
print("Rerunning UMAP")
integrated_obj <- integrated_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,
          verbose = F, reduction.name = "umap_harmony")

conf_plot_harmony <- corr_conf(integrated_obj, assay = "HARMONY")

## Optimise clustering for the integrated object (Very intensive)
print("Made it to clustering optimisation")
gc()
integrated_obj$opt_clust <- NULL
integrated_obj <- optimiseClusters(integrated_obj,
                                      seq(0.5, 1.5, 0.2),
                                      reduction = "harmony")
gc()
print("Completed clustering optimisation")

umap_clust <- DimPlot(integrated_obj, reduction = "umap_original",
                      group.by = "opt_clust", pt.size = .1) +
  labs(title = "Original UMAP Clusters", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

harmony_clust <- DimPlot(integrated_obj, reduction = "umap_harmony",
        group.by = "opt_clust", pt.size = .1) +
  labs(title = "Harmony UMAP Clusters", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

# Save back into the same object
saveRDS(integrated_obj, file = sample_path)

pdf(file = paste0(fig_path, sample_id, "_snRNA_clustering.pdf"),
    width = 10, height = 10)
integrated_PCA
conf_plot_original
conf_plot_harmony
harmony_pca
umap_clust
harmony_clust
dev.off()


