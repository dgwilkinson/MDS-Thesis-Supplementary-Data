library(tidyverse)
library(optparse)
library(viridis)

# Create options list
options <- list(
  make_option(c("--sample_dir"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

# Load metrics_summary.csv file from each sample SpaceRanger output 
samples <- as.vector(t(read.table(paste0(sample_dir,"sample_names.txt"))))
metrics_files <- paste0(sample_dir, samples,
                        "/outs/spatial_", samples, 
                        "_metrics_summary.csv")

spatial_metrics <- map(metrics_files, read.csv) %>%
  enframe(name = "sample_names") %>%
  unnest(cols = everything()) %>%
  mutate(sample_names = samples) %>%
  rename_with(~gsub("[.]", " ", .x))

spatial_qc <- spatial_metrics %>%
  select(c("sample_names",
           "Number of Spots Under Tissue",
           "Median Genes per Spot",
           "Mean Reads per Spot",
           "Median UMI Counts per Spot",
           "Fraction of Spots Under Tissue",
           "Fraction Reads in Spots Under Tissue")) %>%
  arrange(`Median Genes per Spot`) %>%
  pivot_longer(-sample_names, names_to = "metric") %>%
  ggplot(aes(x = factor(sample_names, unique(sample_names)),
             y = value,
             fill = factor(sample_names, unique(sample_names)))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = T, guide = "none") +
  facet_wrap(~metric, ncol = 3, scales = "free_y",
             strip.position = "left") +
  theme_bw() +
  scale_y_continuous(n.breaks = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  labs(x = "Sample ID", y = NULL)

write.table(spatial_metrics,
            file = paste0(fig_path, "all_spatial_metrics.csv"),
            sep = ",")

# Print QC plot to PDF
pdf(file = paste0(fig_path, "spatial_all_qcs.pdf"),
    width = 16, height = 9)
spatial_qc
dev.off()
