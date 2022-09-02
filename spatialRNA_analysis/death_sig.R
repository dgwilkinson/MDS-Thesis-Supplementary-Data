library(Seurat)
library(tidyverse)
library(optparse)
library(patchwork)
library(viridis)
library(decoupleR)

# Create options list
options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--marker_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load Unfiltered spatial sample from Space Ranger output
visium_obj <- Load10X_Spatial(data.dir = sample_path,
                filename = paste0("spatial_",sample_id,"_raw_feature_bc_matrix.h5"),
                filter.matrix = T)

# Load all death signature genesets into a list
gset_paths <- list.files(marker_path, full.names = T)
gset_list <- map(gset_paths, ~ pull(read_tsv(.x, skip=1),1))
names(gset_list) <- str_extract(string = gset_paths, 
                                 pattern = "[^/]*(?=.tsv)")

# Calculate the Module score for all genesets on the unfiltered sample
# visium_obj <- AddModuleScore(visium_obj,
#                              features = gset_list,
#                              name = "DS")
# visium_obj@meta.data <- visium_obj@meta.data %>%
#   rename(setNames(paste0("DS", seq_along(gset_list)), names(gset_list)))

# Calculate GSE using run_wmean from decoupleR instead

# Transform geneset list into acceptable input
wmean_in <- gset_list %>% 
  enframe(name = "source", value = "target") %>%
  unnest(everything()) %>%
  mutate(mor = 1)

# Run Weighted Mean GSEA
wmean_out <- run_wmean(network = wmean_in, 
                       mat = GetAssayData(visium_obj, "data"),
                       times = 100)

# Transform output to select normal distribution z-scores 
wmean_out <- wmean_out %>%
  filter(statistic == "norm_wmean") %>%
  select(-c(statistic, p_value)) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  t() %>%
  as.data.frame()

# Bind GSE scores to the slide metadata
visium_obj@meta.data <- visium_obj@meta.data %>%
  bind_cols(wmean_out)

# Density plot of Module scores across spots
gse_density <- visium_obj@meta.data %>%
  select(c(names(gset_list))) %>%
  pivot_longer(everything(), names_to = "gset", values_to = "score") %>%
  ggplot(aes(x = score, col = gset, fill = gset)) +
  geom_density(alpha = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  labs(y = "Density", x = "Module Score", col = "Geneset", fill = "Geneset")

# Spatial feature plot of module scores
gse_spatial <- SpatialFeaturePlot(visium_obj,
                                   features = names(gset_list),
                                   ncol = 3,  combine = T) & 
  scale_fill_viridis() &
  guides(fill = guide_colourbar(label.theme = element_text(hjust = 0.5, 
                                                           vjust = 0.5,
                                                           size = 8,
                                                           angle = -90)))

# Write Enrichment Scores to csv
write.table(wmean_out, 
            file = paste0(fig_path, sample_id, "_DS_GSE.csv"),
            sep = ",")

# Print Plots to PDF
pdf(file = paste0(fig_path, sample_id, "_death_signature.pdf"),
    width = 15, height = 15)
gse_density
gse_spatial
dev.off()
