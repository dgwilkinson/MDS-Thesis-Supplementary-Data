# This script is a modification on https://github.com/saezlab/visium_heart/blob/b0d9c1e92c20ad680d410bd08d22ed2085324ec1/st_snRNAseq/01.1_spatial_QC_reading/run_singleprocessing_spatial.R

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(optparse)
library(viridis)
library(ggsci)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

# Enforce reproducible clustering 
set.seed(42)

# Create options list
options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--obj_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load Spatial sample from Space Ranger output
visium_obj <- Load10X_Spatial(data.dir = sample_path,
                filename = paste0("spatial_",sample_id,"_raw_feature_bc_matrix.h5"),
                filter.matrix = F)

# Record if barcodes derive from on or off tissue spots
visium_obj$orig.ident <- sample_id
visium_obj$tissue <- ifelse(GetTissueCoordinates(visium_obj,
                                          cols = c("tissue"),
                                          scale = NULL) == 1, 
                            "on_tissue", "off_tissue")


# Compare QC metrics between on and off tissue barcodes

tissue_qc <- visium_obj@meta.data %>%
  select(-c(orig.ident)) %>%
  pivot_longer(-tissue, names_to = "metric") %>%
  ggplot(mapping = aes(x = factor(tissue, unique(tissue)), y = value)) +
  geom_violin(aes(fill = tissue)) +
  scale_fill_viridis(discrete = T, guide = "none") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_log10() +
  facet_wrap(~metric, ncol = 2, scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(nCount_Spatial = "UMI Count",
                                      nFeature_Spatial = "Gene Count"))) +
  labs(x = "Tissue Presence", y = NULL)

gc()

# Remove off-tissue spots 
visium_obj <- subset(visium_obj, subset = tissue == "on_tissue")

# Calculate percentage MT genes
visium_obj$percent.mt <- PercentageFeatureSet(visium_obj,
                                              pattern = "^MT-")
# Filter mitochondrial and ribosomal genes 
mt_rb_bool <- grepl("(^MT-)|(^RPS)|(^MRP)|(^RPL)", row.names(visium_obj))
visium_obj <- visium_obj[!mt_rb_bool,]

# Filter genes present in less than 10 cells + recalculate UMI/gene counts
visium_obj <- visium_obj[rowSums(visium_obj[["Spatial"]]@counts>0) > 10,]
visium_obj$nFeature_filt <- colSums(visium_obj[["Spatial"]]@counts>0)
visium_obj$nCount_filt <- colSums(visium_obj[["Spatial"]]@counts)

# Scatter plot filter metrics 

visium_obj@meta.data <- visium_obj@meta.data %>%
  mutate(passed = nFeature_filt>300 &
           nCount_filt>500)

count_feat_plot <- visium_obj@meta.data %>%
  ggplot(aes(x = nCount_filt, y = nFeature_filt)) +
  theme_classic() +
  geom_point(aes(col = passed)) + 
  scale_color_viridis(discrete = T) +
  geom_vline(xintercept = 500, col = "red", linetype = "dashed") +
  geom_hline(yintercept = 300, col = "red", linetype = "dashed") +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "UMI Count per Spot", y = "Gene Count per Spot", color = "Passed Filters")

visium_obj@meta.data <- visium_obj@meta.data %>%
  mutate(passed = nFeature_filt>300)

feat_MT_plot <- visium_obj@meta.data %>%
  ggplot(aes(x = nFeature_filt, y = percent.mt)) +
  theme_classic() +
  geom_point(aes(col = passed)) + 
  scale_color_viridis(discrete = T) +
  geom_vline(xintercept = 300, col = "red", linetype = "dashed") +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "Gene Count per Spot", y = "Percent MT per Spot", color = "Passed Filters")

visium_obj$passed <- NULL

# Just a way to rename the colour bars in the Spatial Feature Plot
visium_obj@meta.data <- visium_obj@meta.data %>%
  mutate(`UMI Count` = nCount_Spatial,
         `Gene Count` = nFeature_Spatial,
         `MT Percentage (%)` = percent.mt)

# Output formatting uses patchwork syntax
spatial_plot <- SpatialFeaturePlot(visium_obj,
                   features = c("UMI Count",
                                "Gene Count",
                                "MT Percentage (%)"),
                   ncol = 3,  combine = T) & 
  scale_fill_viridis() &
  guides(fill = guide_colourbar(label.theme = element_text(hjust = 0.5, 
                                                           vjust = 0.5,
                                                           size = 8,
                                                           angle = -90)))

# Filter Spatial data
visium_obj <- subset(visium_obj, subset = nFeature_filt > 300 & 
                       nCount_filt > 500)

# Perform SCT normalisation + DR (Shouldn't scale the SCT assay)
visium_obj <- visium_obj %>%
  SCTransform(assay = "Spatial",
              verbose = F) %>%
  RunPCA(npcs = 30,
         verbose = F,
         assay = "SCT") %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          verbose = F)
  
# Optimise Clustering on the visium object
visium_obj <- optimiseClusters(visium_obj, resolutions = seq(0.5, 1.5, 0.1),
                assay = "SCT", reduction = "pca")

spatial_clusters <- DimPlot(visium_obj, group.by = "opt_clust",
                            reduction = "umap") +
  scale_color_npg() +
  labs(title = paste(sample_id, "Spatial Clusters"),
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"),
       col = "Cluster")

spatial_dimplot <- SpatialDimPlot(visium_obj, group.by = "opt_clust") +
  scale_fill_npg(guide = "none") + 
  labs(fill = "Cluster")

feature_plot <- FeaturePlot(visium_obj, 
                            features = c("UMI Count", 
                                         "Gene Count",
                                         "MT Percentage (%)"),
                            ncol = 3,
                            cols = c("lightgrey", viridis(1))) &
  labs(x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension")) 

# Remove dummy variable names
visium_obj$`UMI Count` <- NULL
visium_obj$`Gene Count` <- NULL
visium_obj$`MT Percentage (%)` <- NULL

# Save RDS to out path
Idents(visium_obj) <- "opt_clust"
saveRDS(visium_obj, file = paste0(obj_path,sample_id,".rds"))

# Print Plots to PDF

pdf(file = paste0(fig_path, sample_id, "_initial_qc.pdf"),
    width = 15, height = 15)
tissue_qc/(count_feat_plot+feat_MT_plot)
spatial_plot/feature_plot
plot_grid(spatial_clusters, spatial_dimplot, ncol = 2)
dev.off()
