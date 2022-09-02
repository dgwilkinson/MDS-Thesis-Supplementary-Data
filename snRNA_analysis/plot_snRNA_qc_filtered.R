library(Seurat)
library(tidyverse)
library(viridis)
library(cowplot)

processed_path <- "/home/dwzj28/research_project/dataset/snRNA/processed_objects"
media_path <- "/home/dwzj28/research_project/media"

# Find processed objects and extract sample names
processed_objects <- paste(processed_path, list.files(processed_path, pattern = "^[^_]+$"), sep = "/")
names(processed_objects) <- str_extract(processed_objects,
                                        pattern = "[CKL]{2}[0-9]+")

get.metadata <- function(rds_path){
  seurat_obj <- readRDS(rds_path)
  metadata <- seurat_obj@meta.data %>%
    select(c("orig.ident", "nCount_RNA",
             "nFeature_RNA", "percent.mt",))
  return(metadata)
}

# Bind the metadata dataframes on rows
all_metadata <- map(processed_objects, get.metadata) %>%
  bind_rows()

spatial_metrics <- read.table(file = "~/research_project/media/all_spatial_metrics.csv",
                              sep = ",",check.names = F)

# Calculate qc metrics by Sample
qc_metrics <- all_metadata %>%
  group_by(orig.ident) %>%
  summarise(n_cells = n(),
            median_umi = median(nCount_RNA),
            median_features = median(nFeature_RNA),
            mean_percent_MT = sprintf("%.2f%%",mean(percent.mt))) %>%
  rename("sample_name" = orig.ident)

write.table(qc_metrics, row.names = F, col.name = T, sep = ",",
            file = paste(media_path, "snRNA_filtered_qcs.csv", sep = "/"))

# Violin Plot for qc metrics across samples

metric_names = set_names(c("UMI Count",
                           "Feature Count",
                           "Percentage MT (/%)"),
                           names(all_metadata)[-1])

violin_plots <- lapply(names(all_metadata)[-1], function(metric){
  ggplot(all_metadata, aes(x = orig.ident, y = .data[[metric]])) +
    geom_violin(aes(fill = orig.ident), scale = "width",
                draw_quantiles = seq(0.25, 0.75, by = 0.25),
                col = "black") +
    scale_fill_viridis(discrete = T, guide = "none") +
    theme_classic() +
    scale_y_log10(n.breaks = 10) +
    labs(x = "Sample", y = metric_names[metric]) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
})

# Print to pdf in media out path

pdf(file = paste(media_path, "snRNA_filtered_qcs.pdf", sep = "/"),
    width = 15, height = 15)

plot_grid(plotlist = violin_plots, nrow = 3, ncol = 1)

dev.off()


