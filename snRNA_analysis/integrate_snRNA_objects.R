library(Seurat)
library(harmony)
library(cluster)
library(optparse)
library(tidyverse)
library(readxl)
library(viridis)
library(clustree)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

options <- list(
  make_option(c("--objects_dir"), action = "store", type = "character"),
  make_option(c("--out_path"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--metadata_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Read metadata for patient, sample and batch information
batch_metadata <- read_xlsx(metadata_path, skip = 1) %>%
  select(c("sample_id", "patient", "patient_group"))

# Locate all rds object paths 
objects <- list.files(objects_dir, pattern = "^[^_]*.rds", full.names = T)
names(objects) <- str_extract(objects, "[CKL]{2}[0-9]+")

# Remove KLXXX object paths for which there is no metadata
objects <- objects[names(objects) %in% batch_metadata$sample_id]

# Load in rds objects (memory intensive ~38GB total disk usage in folder)
loaded_objects <- map(objects, readRDS)

print("Successfully loaded all objects")

# Find top 3000 highly variable genes per sample then sort by decreasing occurrence
feature_list <- lapply(loaded_objects, FUN = function(seurat_obj){
  assay <- DefaultAssay(seurat_obj)
  FindVariableFeatures(seurat_obj,
                       selection.method = "vst",
                       nfeatures = 3000,
                       assay = assay)
  seurat_obj[[assay]]@var.features
}) %>% 
  unlist() %>% 
  table() %>%
  sort(decreasing=T)

hvg_dist <- feature_list %>% enframe(name = "gene", value = "count") %>%
  group_by(count) %>%
  mutate(count = as.factor(count)) %>%
  summarise(ngenes_at_count = n()) %>%
  ggplot(mapping = aes(x = count, y = ngenes_at_count)) +
  geom_bar(aes(fill = count), stat = "identity") + 
  scale_fill_viridis(discrete=T, guide = "none") + 
  theme_bw() +
  scale_y_continuous(n.breaks = 10) +
  labs (x = "Common Samples per Highly Variable Gene", 
        y = "Gene Count")

# Filter for top 3000 occurring genes
top_features <- feature_list %>%
  as.data.frame() %>%
  slice(1:3000) %>%
  pull(".")
 
# Calculate original cell counts 
orig_counts <- sapply(loaded_objects, FUN= function(obj){
  nrow(obj@meta.data) 
}) %>% enframe(name = "orig.ident", value = "num_cells") %>%
  arrange(orig.ident)

# Merge the list of loaded objects into a single object
print("Merging objects...")
#integrated_object <- loaded_objects %>% reduce(merge, merge.data = T)

# Modified merging procedure to limit memory usage, merges objects individually 
# then removes from the loaded object list until list empties.
integrated_object <- loaded_objects[[1]]
loaded_objects[[1]] <- NULL
gc()
while (length(loaded_objects)){
  print(paste("Merging sample object:", names(loaded_objects)[1]))
  integrated_object <- merge(integrated_object, loaded_objects[[1]], merge.data = T)
  loaded_objects[[1]] <- NULL
  gc()
}
rm(loaded_objects)

# Check that cell counts have been preserved through the merge
int_counts <- integrated_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(num_cells = n()) %>%
  arrange(orig.ident)

if (identical(int_counts, orig_counts)){
  print("Sample objects successfully merged")
} else {
  print("Integrated object has different number of cells")
  print("Integrated:")
  int_counts
  print("Original:")
  orig_counts
}

# Add batch metadata to the integrated object
print("Adding metadata")
integrated_object@meta.data <- integrated_object@meta.data %>% 
  rownames_to_column() %>%
  left_join(batch_metadata, by = c("orig.ident" = "sample_id")) %>%
  column_to_rownames()

# Run Basic Processing on the integrated object
print("Run Dimensionality Reduction")
integrated_object <- integrated_object %>%
  ScaleData(verbose = F) %>%
  RunPCA(features = top_features,
         verbose = F,
         npcs = 30)  %>%
  RunUMAP(reduction = "pca", dims = 1:30, 
          verbose = F, reduction.name = "umap_original")

# Plot original UMAP and PCAs 
print("Plotting DR")
original_umap <- integrated_object %>% DimPlot(., reduction = "umap_original",
                                               group.by = "orig.ident",
                                               pt.size = .1,
                                               cols = viridis(length(unique(.$orig.ident)))) +
  labs(title = "Uncorrected UMAP", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

original_pca <- integrated_object %>% DimPlot(., reduction = "pca",
                                              group.by = "orig.ident",
                                              pt.size = .1,
                                              cols = viridis(length(unique(.$orig.ident)))) +
  labs(title = "Uncorrected PCA", 
       x = expression(1^st~" Principal Component"),
       y = expression(2^nd~" Principal Component"))

# PCA scree plot 
print("Plotting Scree Plot")
tot_var <- integrated_object[["pca"]]@misc$total.variance
eig_vals <- integrated_object[["pca"]]@stdev^2
pcs <- names(integrated_object[["pca"]])

scree_plot <- data.frame(component = factor(pcs, levels = pcs),
                         explained_var = eig_vals/tot_var * 100) %>%
  ggplot(aes(x = component, y = explained_var)) +
  geom_point(shape = 1, size = 2,  col = viridis(1)) +
  geom_line(group=1, col = viridis(1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  labs(y = "Explained Variance (%)", x = "Principal Component")

# Plot correlation of PCs with confounding variables 
library(scater)

print("Plotting confounding variables for PCA")

conf_variables <- c("orig.ident", "patient")

conf_plot <- as.SingleCellExperiment(integrated_object) %>% 
  getExplanatoryPCs(variables = conf_variables,
                    dimred = "PCA", n_dimred =30) %>%
  as.data.frame() %>%
  rownames_to_column("components") %>%
  pivot_longer(-c(components), names_to = "variable") %>%
  ggplot(mapping = aes(x = factor(components, unique(components)),
                       y = log(value,base = 10), group = variable, col = variable)) +
  geom_line(stat = "identity")+
  geom_point(shape = 1, size = 2) + 
  scale_color_viridis(discrete=T, labels = c("Sample ID", "Patient")) + 
  theme_classic() + 
  scale_y_continuous(breaks = seq(-3, 2),
                     labels = sprintf("%g", 10**seq(-3, 2))) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  labs(x= "Principal Component", y = "Explained Variance (%)", color = "Confounding Variables")

detach("package:scater", unload = T)

# Run Harmony Batch effect correction
print("Running Harmony")
print(DefaultAssay(integrated_object))
integrated_object <- RunHarmony(integrated_object,
                                group.by.vars = conf_variables,
                                reduction = "pca",
                                assay.use = DefaultAssay(integrated_object),
                                plot_convergence = T,
                                max.iter.harmony = 20)

harmony_pca <- DimPlot(integrated_object, 
                       reduction = "harmony", 
                       group.by = "orig.ident",
                       pt.size = .1) +
  labs(title = "Harmony PCA by Sample", 
       x = expression(1^st~" Principal Component"),
       y = expression(2^nd~" Principal Component"))

# Run UMAP on the pca harmony cell embedding
print("Rerunning UMAP")
integrated_object <- integrated_object %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,
          verbose = F, reduction.name = "umap_harmony")

## Optimise clustering for the integrated object (Very intensive)
# print("Made it to clustering optimisation")
# gc()
# integrated_object$opt_clust <- NULL
# integrated_object <- optimiseClusters(integrated_object,
#                                       seq(0.5, 1.5, 0.2),
#                                       reduction = "harmony")
# gc()
# print("Completed clustering optimisation")

# Cluster optimisation requires a lot of time and memory
print("Performing Clustering")
integrated_object$opt_clust <- NULL

integrated_object <- integrated_object %>% 
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = seq(0.5, 1.5, 0.2), verbose = F)

cluster_tree <- clustree(integrated_object,
                         prefix = paste0(DefaultAssay(integrated_object),
                                         "_snn_res."))

integrated_object$opt_clust <- integrated_object$SCT_snn_res.0.9

print("Completed Clustering")

harmony_clust <- DimPlot(integrated_object, reduction = "umap_harmony",
        group.by = "opt_clust", pt.size = .1) +
  labs(title = "Harmony UMAP Clusters", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

umap_clust <- DimPlot(integrated_object, reduction = "umap_original",
        group.by = "opt_clust", pt.size = .1) +
  labs(title = "Uncorrected UMAP Clusters", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

saveRDS(integrated_object, file = paste(out_path, "integrated_SCT_harmony.rds", sep = "/"))

pdf(file = paste(fig_path, "integrated_snRNA_plots_nopg.pdf", sep = "/"),
    width = 10, height = 10)
hvg_dist
original_pca
original_umap
scree_plot
conf_plot
cluster_tree
harmony_pca
umap_clust
harmony_clust
dev.off()


