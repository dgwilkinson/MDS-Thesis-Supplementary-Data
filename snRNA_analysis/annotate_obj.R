library(Seurat)
library(optparse)
library(tidyverse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--anno_path"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load in integrated object and annotations 
integrated_obj <- readRDS(sample_path)
annos <- read.table(anno_path, 
                    sep = ",",
                    header = T) %>% # Had to remove a comma in one of the cell type names
  remove_rownames() %>%
  column_to_rownames(var = "opt_clust")

# Evaluate the non-definitive annotations
if (sample_id == "integrated_DS"){
  # Select Fibroblast and CM marker genesets from hca
  fib_and_cm <- read.table("~/research_project/markers/hca_ct_markers.csv", sep = ",") %>%
    mutate(gene = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                        keys = gene, 
                                        keytype = "ENSEMBL", 
                                        column = "SYMBOL",
                                        multiVals = "first")) %>%
    filter(p_val_adj < 0.005 & !is.na(gene), 
           cluster %in% c("regular ventricular cardiac myocyte",
                          "fibroblast")) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    slice(1:150) %>%
    summarise(gene_list = list(gene)) %>%
    deframe()
  
  # Subset out optimised clusters 5 and 6 and calculate Module score for each geneset
  clust_56 <- subset(integrated_obj, opt_clust %in% c(5,6)) %>%
    AddModuleScore(features = fib_and_cm,
                   assay = "integrated",
                   search = F) 
  
  # Average cell type score across clusters 
  ct_results <- clust_56@meta.data %>% 
    group_by(opt_clust) %>%
    summarise(`Fibroblast` = mean(Cluster1), 
              `Cardiac Myocycte`= mean(Cluster2)) 
  
  # Plot values in bar plot
  cm_fib_plot <- ct_results %>%
    pivot_longer(-opt_clust, names_to = "cell_type_score") %>%
    ggplot(aes(x=opt_clust, y = value, fill = cell_type_score)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_bw() +
    labs(x = "Cluster", y = "Score", fill = "Cell Type")
  
  # Write results to table
  write.table(ct_results,
              file = paste0(fig_path, "nondefinitive_ct_scores.csv"),
              sep = ",")
}

# Transfer annotations to the integrated object metadata
integrated_obj@meta.data <- integrated_obj@meta.data %>%
  mutate(cell_type = annos[as.character(opt_clust),]) %>%
  mutate_at("cell_type", ~factor(.x, levels = c("regular ventricular cardiac myocyte",
                                                "macrophage",
                                                "fibroblast",
                                                "endothelial cell of artery",
                                                "capillary endothelial cell",
                                                "vein endothelial cell",
                                                "neural cell",
                                                "CD8-positive alpha-beta cytotoxic T cell",
                                                "pericyte cell",
                                                "smooth muscle cell",
                                                "dendritic cell")))

# Set cell type as the new identity class and DimPlot
Idents(integrated_obj) <- "cell_type"
ct_umap <- DimPlot(integrated_obj, group.by = "cell_type",
        reduction = "umap_harmony") +
  labs(title = paste(sample_id, "HCA Cell Type Annotations"),
       x = expression(1^st~" UMAP Dimension"),
       y = expression(2^nd~" UMAP Dimension"),
       col = "Cell Type") 

# Save annotated object back to the sample path
saveRDS(integrated_obj, file = sample_path)

pdf(file = paste0(fig_path, sample_id, "_annotated_dim_plot.pdf"),
    width = 15, height = 15)
if (sample_id == "integrated_DS"){
  cm_fib_plot
}
ct_umap
dev.off()



