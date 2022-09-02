library(Seurat)
library(optparse)
library(viridis)
library(tidyverse)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"))

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

integrated_metadata <- readRDS(sample_path)@meta.data
gc()

cells_by_sample <- integrated_metadata %>%
  group_by(orig.ident) %>%
  summarise(ncells = n(),
            patient_group = unique(patient_group)) %>%
  arrange(desc(patient_group), ncells) %>%
  ggplot(aes(x = factor(orig.ident, orig.ident),
             fill = patient_group,
             y = ncells)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = T, labels = c("Control/RZ", "IZ/BZ", "FZ")) +
  scale_y_continuous(n.breaks = 10) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, 
                                   vjust = 0.5,
                                   hjust = 0)) +
  labs(x = "Sample ID", y = "Contributed Cells", fill = "Region")

viol_plot <- integrated_metadata %>%
  select(c(orig.ident, nFeature_RNA, nCount_RNA, percent.mt)) %>%
  arrange(nFeature_RNA) %>%
  pivot_longer(-orig.ident, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = factor(orig.ident, unique(orig.ident)),
             y = value,
             fill = factor(orig.ident, unique(orig.ident)))) +
  geom_violin(scale = "width") +
  theme_classic() +
  facet_wrap(~metric, scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(c(nFeature_RNA = "Gene Count",
                                      nCount_RNA = "UMI Count",
                                      percent.mt = "Percentage MT (%)")),
             ncol = 1, nrow = 3) +
  scale_fill_viridis(discrete = T, guide = "none") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -90,
                                   vjust = 0.5,
                                   hjust = 0)) +
  labs(x = "Sample ID", y = NULL)

pdf(file = paste0(fig_path, sample_id, "_qc_plots.pdf"),
    width = 10, height = 10)
cells_by_sample
viol_plot
dev.off()



