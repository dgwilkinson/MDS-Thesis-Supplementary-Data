library(tidyverse)
library(optparse)
library(readxl)
library(viridis)

# Create options list
options <- list(
  make_option(c("--fig_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load filtered qc metrics and metadata for spatial and snRNA experiments
snRNA_metadata <- read_xlsx(path = "~/research_project/dataset/snRNA/metadata-RNA-seq.xlsx",
                            skip = 1)
spatial_metadata <- read.csv(file = "~/research_project/dataset/spatialRNA/metadata-Visium.csv")
snRNA_qc <- read.csv(paste0(fig_path, "snRNA_filtered_qcs.csv"))
spatial_qc <- read.csv(paste0(fig_path, "spatial_filtered_qcs.csv"))

# Attach region id from respective metadata files to qc metrics and group n_cells/spots upon it
n_cells <- left_join(snRNA_qc, snRNA_metadata, 
                     by = c("sample_name" = "sample_id")) %>%
  filter(!is.na(patient_region_id)) %>%
  group_by(patient_region_id) %>%
  summarise(n_cells = sum(n_cells.x),
            patient_group = patient_group[1])
spot_loss <- left_join(spatial_qc, spatial_metadata, 
                     by = c("sample_names" = "hca_sample_id")) %>%
  filter(!is.na(patient_region_id)) %>%
  group_by(patient_region_id) %>%
  summarise(n_spots = sum(n_spots),
            nspots_under_tissue = sum(nspots_under_tissue)) %>%
  mutate(loss_percent = 1 - n_spots/nspots_under_tissue) 

# Plot number of cells against spot loss coloured by patient group

loss_plot <- merge(n_cells, spot_loss) %>%
  ggplot(aes(x = 100 * loss_percent, y = n_cells, col = patient_group)) +
  theme_classic() +
  geom_point(pch = 4, size = 3, stroke = 1.25) +
  geom_smooth(aes(col = NULL), method = "lm", formula = y~x, 
              se = F, linetype = "dashed", col = "red") +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_viridis(discrete = T, labels = c("Control/RZ", "IZ/BZ", "FZ")) +
  labs(x = "Percentage Spot Loss (%)", y = "Number of Cells", col = "Region")

# Print QC plot to PDF
pdf(file = paste0(fig_path, "loss_plot.pdf"),
    width = 10, height = 10)
loss_plot
dev.off()
