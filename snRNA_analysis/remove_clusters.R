library(Seurat)
library(tidyverse)
library(viridis)
library(scales)
library(optparse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

seurat_obj <- readRDS(sample_path)
# Calculate Features by cluster 
clust_metrics <- seurat_obj@meta.data %>% 
  mutate(complexity = log10(nFeature_RNA)/log10(nCount_RNA)) %>%
  group_by(opt_clust) %>%
  summarise(mean_complexity = mean(complexity),
            mean_feat = mean(nFeature_RNA),
            mean_count = mean(nCount_RNA))%>% 
  arrange(mean_feat, decreasing = F)

complexity_threshold = 0.8
count_threshold = quantile(seurat_obj$nCount_RNA, 0.05)
feature_threshold = quantile(seurat_obj$nFeature_RNA, 0.05)

# Get Cluster Colours
clust_pal = hue_pal()(length(unique(seurat_obj$opt_clust)))

# Plot bar plots of cluster metrics
metric_bar <- clust_metrics %>% 
  pivot_longer(-opt_clust, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = factor(opt_clust, unique(opt_clust)),
             y = value)) +
  geom_bar(aes(fill = opt_clust), stat = "identity") +
  geom_hline(data = data.frame(metric = c("mean_complexity",
                                          "mean_count",
                                          "mean_feat"),
                               yint = c(complexity_threshold,
                                        count_threshold,
                                        feature_threshold)),
             aes(yintercept = yint),
             col = "red", linetype = "dashed") +
  scale_fill_manual(values = clust_pal, guide = "none") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous(n.breaks = 6) +
  facet_wrap(.~metric, scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(mean_complexity = "Mean Complexity",
                                      mean_count = "Mean UMI Count",
                                      mean_feat = "Mean Gene Count")),
             ncol = 1) +
  labs(x = "Optimised Cluster", y = NULL)

# Filter Seurat Object for Poor Quality Clusters 
# (Complexity <0.8 or nUMI & nFeature below 5% quantile)
poor_clusters <- clust_metrics %>%
  filter(mean_complexity < complexity_threshold | 
           (mean_count < count_threshold & mean_feat < feature_threshold)) %>%
  pull(opt_clust)

seurat_obj <- seurat_obj[,!(seurat_obj$opt_clust %in% poor_clusters)]
saveRDS(sample_path)

pdf(paste0(fig_path, sample_id, "_cluster_filters.pdf"), 
    width = 8, height = 10)
metric_bar
dev.off()
