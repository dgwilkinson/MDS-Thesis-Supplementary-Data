library(optparse)
library(tidyverse)
library(viridis)
library(cowplot)

options <- list(
  make_option("--data_path", action = "store", type = "character"),
  make_option("--fig_path", action = "store", type = 'character')
)

#Parse Command Line Arguments into equivalently named variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir = globalenv()))

# Find all sample directories on the data path
sample_dirs <- list.dirs(data_path) %>%
  as.vector() %>%
  .[grepl(paste0(data_path, "/[CKL]{2}[0-9]{3}$"), .)]

# Locate the metadata metrics for each sample 
sample_metadata <- sapply(sample_dirs, function(dir){
  paste(dir, list.files(dir, pattern="metrics_summary.csv", recursive = T), sep = "/")
})

# Extract sample names 
names(sample_metadata) <- str_extract(names(sample_metadata), "[CKL]{2}[0-9]{3}")

# KL named samples have an extra column (20 vs 19)
CK_samples <- sample_metadata %>%
  enframe(name = "sample_name", value = "metadata_path") %>%
  filter(str_detect(sample_name, "KL", negate = T)) %>%
  mutate(read_csv(metadata_path))

KL_samples <- sample_metadata %>%
  enframe(name = "sample_name", value = "metadata_path") %>%
  filter(!sample_name %in% CK_samples$sample_name) %>%
  mutate(read_csv(metadata_path))

# Extra column refers to the Phred Quality score of second read
#setdiff(names(KL_samples), names(CK_samples))

# Remove extra column and combine the samples 
KL_samples <- KL_samples %>%
  select(-c("Q30 Bases in RNA Read 2"))

all_samples <- bind_rows(CK_samples, KL_samples)

# Select and transform qc metrics
all_samples_qc <- all_samples %>% 
  select(c("sample_name",
           "Estimated Number of Cells",
           "Mean Reads per Cell",
           "Median Genes per Cell",
           "Fraction Reads in Cells",
           "Total Genes Detected",
           "Median UMI Counts per Cell")) %>%
  mutate_at("Fraction Reads in Cells", ~as.double(str_replace(.x, "[%]", "")))

qc_plots <- lapply(names(all_samples_qc)[-1], function(metric){
  ggplot(data = all_samples_qc, mapping = aes(x = sample_name, y = .data[[metric]])) +
    theme_bw() +
    geom_bar(aes(fill = sample_name), stat = "identity") +
    scale_fill_viridis(discrete = T, guide = "none") +
    scale_y_continuous(n.breaks = 10) +
    labs(x = "Sample", y = metric) +
    theme(axis.text.x  = element_text(angle = -90, hjust = 0, vjust = 0.5),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 13))
})

pdf(height = 10, width = 15, file = paste0(fig_path,"snRNA_all_qcs.pdf"))
plot_grid(plotlist = qc_plots, nrow = 2, ncol = 3)
dev.off()

write.table(all_samples, file = paste0(fig_path, "snRNA_all_qc_metrics.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")
