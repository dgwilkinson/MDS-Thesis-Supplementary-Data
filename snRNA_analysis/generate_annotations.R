library(optparse)
library(scales)
library(tidyverse)

options <- list(
  make_option(c("--sample_marker_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--anno_marker_path"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Read in annotation & sample markers
anno_markers <- read.table(anno_marker_path, sep = ",")
sample_markers <- read.table(sample_marker_path, sep=",")

# Select the top 150 annotation markers in each cluster (cell type) by log FC
# Also translate markers from ENSEMBL ID to SYMBOL using org.Hs.eg.db
group_gsets <- anno_markers %>%
  mutate(gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                     keys = gene, 
                                     keytype = "ENSEMBL", 
                                     column = "SYMBOL",
                                     multiVals = "first")) %>%
  filter(p_val_adj < 0.005 & !is.na(gene)) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:150) %>%
  summarise(gene_list = list(gene)) %>%
  deframe()

# Filter for the top 300 sample markers in each cluster by log FC
sample_markers <- sample_markers %>%
  filter(p_val_adj < 0.005) %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
  group_by(cluster) %>%
  slice(1:300) %>%
  summarise(gene_list = list(gene)) %>%
  deframe()

# Perform Gene Set Enrichment analysis
gse_res <- map(sample_markers, function(clust_markers){
  anno_genes <- unique(unlist(group_gsets))
  # Restrict cluster markers to only those present in annotation markers
  clust_markers <- intersect(clust_markers, anno_genes)
  # Hypergeometric test for marker enrichment from each gene set
  lapply(group_gsets, function(gset){
    # q - number of cluster markers shared with gene set
    n_shared = length(intersect(clust_markers, gset))
    # m - number of annotation markers in the gene set
    n_gset = length(gset)
    # n - number of annotation markers not in the gene set
    n_xgset = length(anno_genes) - n_gset
    # k - number of cluster markers
    n_clust = length(clust_markers)
    pval = phyper(q = n_shared, 
                  m = n_gset,
                  n = n_xgset,
                  k = n_clust,
                  lower.tail = F,
                  log.p = F)
    return(data.frame(n_shared = n_shared,
             n_gset = n_gset,
             n_xgset = n_xgset,
             n_clust = n_clust,
             pval = pval)) 
  }) %>% 
    enframe(name = "geneset") %>% 
    unnest(cols = value) %>%
    mutate(adj_pval = p.adjust(pval, method = "BH")) %>%
    arrange(adj_pval, desc(n_shared))
}) %>% 
  enframe(name = "cluster") %>% 
  unnest(cols = value)

dot_plot <- gse_res %>%
  group_by(cluster) %>%
  mutate(best = adj_pval== min(adj_pval)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(as.numeric(cluster)), 
             y = geneset, 
             size = -log10(adj_pval),
             fill = as.factor(as.numeric(cluster)),
             col = best)) +
  geom_point(pch = 21, stroke = 1.25) +
  scale_fill_manual(values = hue_pal()(length(unique(gse_res$cluster)))) +
  scale_color_manual(values = c("white", "black")) +
  theme_classic() +
  guides(col = "none", fill = "none") +
  labs(x = "Optimised Cluster", y = "Annotation Cell Type", 
       size = expression(log[10]~"p-value"))

proposed_annos <- gse_res %>% 
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>%
  select(cluster, geneset) %>%
  rename("opt_clust" = cluster,
         "cell_type" = geneset)

write.table(proposed_annos,
            file = paste0(fig_path, sample_id, "_cluster_annotations.csv"),
            quote = F,
            sep = ",")

pdf(file = paste0(fig_path, sample_id, "_annotation_dot_plot.pdf"),
    width = 10, height = 10)
dot_plot
dev.off()

