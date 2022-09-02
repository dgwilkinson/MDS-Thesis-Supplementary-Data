library(tidyverse)
library(optparse)
library(readxl)
library(viridis)

# Create options list
options <- list(
  make_option(c("--sample_ids"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load Death Signature GSEA results generated in death_sig.R (barcode x)
sample_names <- as.vector(t(read.table(file = sample_ids)))
gse_scores <- map(sample_names, ~read.csv(paste0(fig_path, .x, "/", .x, "_DS_GSE.csv")))

# Load spatial metadata for patient region information 
spatial_metadata <- read.csv(file = "~/research_project/dataset/spatialRNA/metadata-Visium.csv") %>%
  select(hca_sample_id, major_labl)

# Join region information to nested gse score matrices
names(gse_scores) <- sample_names
gse_scores <- gse_scores %>% 
  enframe(name = "sample_name", value = "gse_mat") %>%
  left_join(spatial_metadata, by = c("sample_name" = "hca_sample_id"))

# Plot enrichment score by region 
enrich_region <- gse_scores %>% 
  unnest(cols = gse_mat) %>% 
  pivot_longer(-c(sample_name, major_labl), 
               names_to = "gset", 
               values_to = "score") %>%
  ggplot(aes(x = factor(major_labl, levels = c("CTRL", "RZ", "BZ", "IZ", "FZ")),
             y = score, col = gset, fill = gset)) +
  geom_boxplot(position = "dodge", alpha = 0.2, outlier.size = 0.1) +
  scale_y_continuous(n.breaks = 15, limits = c(-1, 20)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  labs(y = "Enrichment Z-Score", x = "Region", col = "Geneset", fill = "Geneset")

# # Alternatively plot mean enrichment score across sample by region 
# gse_scores %>%
#   unnest(cols = gse_mat) %>%
#   group_by(sample_name) %>%
#   filter(!if_any(everything(), is.infinite)) %>%
#   summarise(across(where(is.numeric), mean), major_labl = major_labl[1]) %>%
#   pivot_longer(-c(sample_name, major_labl), 
#                names_to = "gset", 
#                values_to = "score") %>%
#   ggplot(aes(x = factor(major_labl, levels = c("CTRL", "RZ", "BZ", "IZ", "FZ")),
#              y = score, col = gset, fill = gset)) +
#   geom_boxplot(position = "dodge", alpha = 0.2, outlier.size = 0.1) +
#   scale_color_viridis(discrete = T) +
#   scale_fill_viridis(discrete = T) +
#   theme_bw() +
#   labs(y = "Enrichment Z-Score", x = "Region", col = "Geneset", fill = "Geneset")

# Print QC plot to PDF
pdf(file = paste0(fig_path, "region_death_signature.pdf"),
    width = 10, height = 10)
enrich_region
dev.off()
