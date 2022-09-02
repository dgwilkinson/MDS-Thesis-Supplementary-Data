library(Seurat)
library(tidyverse)
library(optparse)
library(patchwork)
library(cowplot)
library(viridis)

# Create options list
options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--abd_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Accidentally deleted results for this sample 
if (sample_id == "ACH0010"){
  quit()
}

# Load sample and cell abundances per spot
ct_mean_abd <- read.csv(abd_path, row.names = 1, check.names = F)
cell_types <- colnames(ct_mean_abd)
visium_obj <- readRDS(sample_path)

# Calculate cell type proportions
ct_props <- ct_mean_abd / rowSums(ct_mean_abd)

# Add cell type proportions to the spot metadata 
visium_obj@meta.data <- visium_obj@meta.data %>%
  bind_cols(ct_props)

# Calculate maximum proportion across all values for use in colorbar limits
max_prop <- visium_obj@meta.data %>%
  select(cell_types) %>%
  max() %>%
  plyr::round_any(., 0.1, f = ceiling)

# Spatial Feature Plot
plot_list = SpatialFeaturePlot(visium_obj, features = cell_types, combine = F)

plot_list = map(plot_list, function(plot){
  plot +
    scale_fill_viridis(limits = c(0,max_prop)) +
    guides(fill = guide_colourbar(label.theme = element_text(hjust = 0.5, 
                                                             vjust = 0.5,
                                                             size = 8,
                                                             angle = -90)))
})

plot_lists <- split(plot_list, cut(1:12, 3))
grouped_plots <- map(plot_lists, function(plot_list){
  plot_grid(plotlist = plot_list, ncol = 2)
})

# Print Plots to PDF
pdf(file = paste0(fig_path, sample_id, "_cell_abundances.pdf"),
    width = 10, height = 10)
grouped_plots
dev.off()
