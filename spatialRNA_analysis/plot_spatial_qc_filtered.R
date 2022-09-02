library(Seurat)
library(tidyverse)
library(viridis)
library(cowplot)

processed_path <- "/home/dwzj28/research_project/dataset/spatialRNA/processed_objects/"
media_path <- "/home/dwzj28/research_project/media"

# Find processed objects and extract sample names
processed_objects <- paste0(processed_path, list.files(processed_path,
			    pattern = "^.*\\.rds$"))
names(processed_objects) <- str_extract(processed_objects,
                                        pattern = "[\\w\\d]*(?=\\.rds)")

get.metadata <- function(rds_path){
  seurat_obj <- readRDS(rds_path)
  metadata <- seurat_obj@meta.data %>%
    select(c("orig.ident", "nCount_filt",
             "nFeature_filt", "percent.mt",))
  gc()
  return(metadata)
}

# Bind the metadata dataframes on rows
all_metadata <- map(processed_objects, get.metadata) %>%
  bind_rows()

# Extract Number of Spots Under Tissue from the spatial metrics collated 
# in initial_spatial_qc.R (used later to calculate spot loss)
tissue_spots <- read.table(file = "~/research_project/media/all_spatial_metrics.csv",
                              sep = ",", check.names = F) %>%
  select(c(sample_names,`Number of Spots Under Tissue`)) %>%
  rename(nspots_under_tissue = `Number of Spots Under Tissue`)

# Calculate qc metrics by Sample
qc_metrics <- all_metadata %>%
  group_by(orig.ident) %>%
  summarise(n_spots = n(),
            median_umi_spot = median(nCount_filt),
            median_features_spot = median(nFeature_filt),
            mean_percent_MT_spot = sprintf("%.2f%%",mean(percent.mt))) %>%
  rename("sample_names" = orig.ident)

# Append the tissue spots column to qc_metrics
qc_metrics <- left_join(qc_metrics, tissue_spots, by = "sample_names")

# Save qc_metrics to csv
write.table(qc_metrics, row.names = F, col.name = T, sep = ",",
            file = paste(media_path, "spatial_filtered_qcs.csv", sep = "/"))

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
    scale_y_log10() +
    labs(x = "Sample", y = metric_names[metric]) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
})

violin_plots <- c(violin_plots[1:2],
		  list(violin_plots[[3]] + scale_y_continuous(n.breaks = 10)))

# Print to pdf in media out path
pdf(file = paste(media_path, "spatial_filtered_qcs.pdf", sep = "/"),
    width = 15, height = 15)

plot_grid(plotlist = violin_plots, nrow = 3, ncol = 1)

dev.off()


