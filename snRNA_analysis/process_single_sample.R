# This script is a modification on https://github.com/saezlab/visium_heart/blob/b0d9c1e92c20ad680d410bd08d22ed2085324ec1/st_snRNAseq/01_snuc_QC_doublets_majorannotation/run_singleprocessing.R

library(Seurat)
library(tidyverse)
library(cowplot)
library(optparse)
library(viridis)
library(scDblFinder)
source("/home/dwzj28/research_project/src/optimise_clustering.R")

# Enforce reproducible clustering 
set.seed(42)

# Create options list
options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--obj_path"), action = "store", type = "character"),
  make_option(c("--dis_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

dis_path <-"/home/dwzj28/research_project/markers/dissociation_markers.csv"
dis_gset <- read.table(dis_path, header = TRUE, sep=",")

# Geneset natively arranged by increasing PValue
# Select top 100 markers as gene_symbols by PValue
top_markers <- dis_gset %>% 
  select(gene_symbol, PValue) %>%
  slice(1:100) %>%
  pull(gene_symbol)

# Process single sample snRNA-seq data
print(paste("Processing sample:", sample_id))

# sample_path <- "/home/dwzj28/research_project/dataset/snRNA/sample_data/CK158/outs/unzipped"
sample_obj <- CreateSeuratObject(counts = Read10X(sample_path),
                                 project = sample_id,
                                 min.cells = 10)

gc()

#Plot Distribution of MT Percentage with 5% cutoff
sample_obj[["percent.mt"]] <- PercentageFeatureSet(sample_obj, pattern = "^MT-")

mtp_hist <- ggplot(sample_obj[["percent.mt"]]) +
  theme_classic() +
  geom_histogram(aes(log(percent.mt)),fill = viridis(1),  bins = 30) +
  geom_vline(aes(xintercept = log(5)), 
             linetype = "dashed", col = "red") +
  labs(x = expression(log[10]*" of MT percentage"), y = "Frequency")

# Barcodes above the 99th quantile of Feature counts are likely to be doublets
upper_feature_threshold <- quantile(sample_obj$nFeature_RNA, 0.99)

sample_obj@meta.data <- sample_obj@meta.data %>%
  mutate(passed = nFeature_RNA<upper_feature_threshold &
           nFeature_RNA>300 &
           nCount_RNA>500)

feat_count_plot <- sample_obj@meta.data %>%
  ggplot(aes(x=nCount_RNA, y= nFeature_RNA)) +
  theme_classic() +
  geom_point(aes(col=passed)) +
  scale_color_viridis(discrete=TRUE) +
  geom_hline(yintercept = 300, linetype = "dashed") + 
  geom_hline(yintercept = upper_feature_threshold, linetype = "dashed") +
  geom_vline(xintercept = 500, linetype = "dashed")+
  labs(x = "Read Count", y="Feature Count", colour = "Passed Filters")

sample_obj@meta.data <- sample_obj@meta.data %>%
  mutate(passed = nCount_RNA>500 & percent.mt<5 )

count_MT_plot <- sample_obj@meta.data %>%
  ggplot(aes(x=nCount_RNA, y = percent.mt)) +
  theme_classic() +
  geom_point(aes(col=passed)) +
  scale_color_viridis(discrete=TRUE) +
  geom_hline(yintercept= 5, linetype = "dashed") +
  geom_vline(xintercept = 500, linetype = "dashed") +
  labs(x = "Read Count", y = "Mitochondrial Percentage (%)")

sample_obj@meta.data$passed <- NULL

# Perform basic filtering
sample_obj <- subset(sample_obj,
                     subset = 300<nFeature_RNA &
                       nFeature_RNA<upper_feature_threshold &
                       percent.mt<5 & nCount_RNA>500)

# Evaluate likely doublets
doublets <- scDblFinder(sce = as.matrix(GetAssayData(sample_obj, "counts")))
sample_obj$doublet_score <- doublets$scDblFinder.score
sample_obj$doublet <- doublets$scDblFinder.class

# Normalise + DR
sample_obj <- sample_obj %>%
  NormalizeData(normalization.method = 'LogNormalize',
                scale.factor = 1e4,
                verbose = F) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 3e3,
                       verbose = F) %>%
  ScaleData(.,
            features = rownames(.), 
            verbose = F) %>%
  RunPCA(verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)

gc()

# DimPlot doublet class assignments 
doublet_plot <- DimPlot(sample_obj, group.by= "doublet", 
        cols = viridis(2), reduction = "umap") +
  labs(title = "Doublet Assignments", 
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

# Score gene expression against dissociation marker geneset
sample_obj <- AddModuleScore(sample_obj,
                             assay = "RNA",
                             features = list("dissociation_score" = top_markers),
                             name = "dis_score")

# Visualise Dissociation Score distribution 
sample_obj@meta.data <- sample_obj@meta.data %>%
  rename(c("dis_score" = "dis_score1"))

upper_dis_threshold <- quantile(sample_obj$dis_score, 0.99)

dis_hist <- sample_obj@meta.data %>%
  ggplot(mapping = aes(x = dis_score)) +
  theme_classic() +
  geom_histogram(fill = viridis(1), bins = 30) +
  geom_vline(xintercept = upper_dis_threshold, 
             linetype = "dashed", col = "red") + 
  labs(x = "Dissociation Score", y = "Frequency")

# Feature Plots on latent dimension (UMAP)
feature_plots <- FeaturePlot(sample_obj, features = c("nCount_RNA", 
                                                      "nFeature_RNA", 
                                                      "percent.mt",
                                                      "doublet_score",
                                                      "dis_score"),
                             cols = c("lightgrey", viridis(1)))

# Filter Out likely doublets and high dissociation score barcodes
sample_obj <- subset(sample_obj,
                      subset = dis_score < upper_dis_threshold &
                        doublet == "singlet")

# Repeat dimensionality reduction on sample data
sample_obj <- sample_obj %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 3e3,
                       verbose = F) %>%
  ScaleData(.,
            features = rownames(.), 
            verbose = F) %>%
  RunPCA(verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)

#Call cluster optimisation function from optimise_clustering script
resolutions <- seq(0.5, 1.5, 0.1)
sample_obj <- optimiseClusters(sample_obj, resolutions = resolutions, assay = "RNA")

# Visualise clusters in the UMAP embedding
clust_embed <- DimPlot(sample_obj, group.by = "opt_clust") +
  labs(title = paste(sample_id, "Optimised Clusters"),
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"))

# Save the Seurat object to the out path
Idents(sample_obj) <- "opt_clust"
saveRDS(sample_obj, paste0(obj_path, sample_id, ".rds"))

#Save plots as pdf
pdf(file = paste0(fig_path, sample_id, "_initial_QC.pdf"), 
    width = 15, height = 15)

plot_grid(nrow = 2, ncol = 2, mtp_hist, dis_hist, count_MT_plot, feat_count_plot)
plot(doublet_plot)
plot(feature_plots)
plot(clust_embed)

dev.off()


