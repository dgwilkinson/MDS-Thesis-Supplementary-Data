# Default markers taken from https://github.com/saezlab/visium_heart/blob/b0d9c1e92c20ad680d410bd08d22ed2085324ec1/st_snRNAseq/01_snuc_QC_doublets_majorannotation/plot_knownmarkers.R

library(Seurat)
library(optparse)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(scales)

options <- list(
  make_option(c("--sample_path"), action = "store", type = "character"),
  make_option(c("--sample_id"), action = "store", type = "character"),
  make_option(c("--fig_path"), action = "store", type = "character"),
  make_option(c("--marker_path"), action = "store", type = "character", default = NA),
  make_option(c("--marker_id"), action = "store", type = "character", default = "default"),
  make_option(c("--assay"), action="store", type = "character", default = "RNA")
)

# Read command-line arguments given to the Rscript and assign to variables
args <- parse_args(OptionParser(option_list = options))
walk(names(args), ~ assign(.x, args[[.x]], envir=globalenv()))

if (is.na(marker_path)){
  markers <- list(
    fibroblasts = c("PDGFRA","C7","COL15A1"),
    pericytes = c("RGS5","ABCC9","KCNJ8"),
    adipocytes = c("GPAM","LEP","FASN"),
    smooth_muscle = c("MYH11", "TAGLN", "ACTA2"),
    neuronal = c("PLP1","NRXN1","NRXN3"),
    endothelial = c("VWF","PECAM1","CDH5"),
    atrial_cardio = c("NPPA","MYL7","MYL4"),
    ventricular_cardio = c("MYH7","MYL2","PRDM16"),
    myeloid = c("CD14","C1QA","FOLR2"),
    lymphoid = c("CD8A","NKG7","LCK")) %>%
    enframe(name = "cell_type", value = "marker_gene") %>%
    unnest(cols = everything()) %>% 
    select("marker_gene", "cell_type")
}else{
  # May use this later to load a predefined list of markers
  markers <- readRDS(marker_path) %>%
    enframe(name="cell_type", value = "marker_gene")%>%
    unnest(cols = everything()) %>%
    select("marker_gene", "cell_type")
}

# Load Gene Count Matrix (GCM) and subset markers to features present within it
seurat_obj <- readRDS(sample_path)
print("Successfully Loaded object")
labels <- seurat_obj$opt_clust
GCM <- seurat_obj[[assay]]@scale.data
print("Loaded object into GCM")
rm(seurat_obj)
absent_markers <- setdiff(markers$marker_gene, rownames(GCM))
#markers <- markers %>% filter(marker_gene %in% rownames(GCM))

GCM <- as.data.frame(GCM)
GCM[absent_markers,] <- NA

# Reorder cell observations according to cluster assignments & select markers

ord <- order(labels, decreasing = F)
GCM <- GCM[markers$marker_gene,ord]
labels <- labels[ord]

#clust_cols <- setNames(viridis(length(levels(labels))), unique(sort(labels)))

aGCM <- GCM %>% t() %>% 
  as.data.frame() %>% 
  mutate(clusters = labels) %>% 
  group_by(clusters)%>% 
  summarise(across(everything(), mean)) %>% 
  select(-clusters) %>% t()
print("aggregated GCM")

colnames(aGCM) <- paste0("cluster_", 0:(ncol(aGCM)-1))
cluster_pal <- hue_pal()(ncol(aGCM))
clust_cols <- setNames(cluster_pal, colnames(aGCM))

ann <- HeatmapAnnotation(df = data.frame("cluster" = colnames(aGCM)),#labels), 
                         col = list("cluster"=clust_cols), 
                         show_annotation_name = F, show_legend = F,
                         annotation_legend_param = list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp = gpar(fontsize = 10),
                                                        labels_gp = gpar(fontsize = 10)))


heatmap <- Heatmap(aGCM, cluster_rows = F, cluster_columns = F,
                   col = c("white", viridis(1)), top_annotation = ann, bottom_annotation = ann,
                   row_split = factor(markers$cell_type, levels = unique(markers$cell_type)),
                   row_title_rot = 0,
                   column_split = factor(colnames(aGCM), colnames(aGCM)),
                   column_title_rot = -90, cluster_column_slices = FALSE,
                   show_column_names = F, row_names_side = "left",
                   show_heatmap_legend = F, row_names_gp = gpar(fontsize=7))

pdf(paste0(fig_path, sample_id, "_marker_heatmap_", marker_id, ".pdf"), 
    width = 10, height = 10)
heatmap
dev.off()

