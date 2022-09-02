library(Seurat)
library(dplyr)
library(cluster)

optimiseClusters <-function(seurat_obj, resolutions, assay=character(1), reduction="pca"){
  # Resolutions should be a vector sequence of clustering resolutions e.g. seq(0.2, 1, 0.1)
  # Assay should be a character reference to an existing assay on the Seurat object
  
  if (!assay %in% names(seurat_obj)){
    assay <- DefaultAssay(seurat_obj)
  }
  
  seurat_obj <- FindNeighbors(seurat_obj,
                              reduction= reduction, 
                              dims=1:30)
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = resolutions)
  
  # Optimise upon cluster silhouette score
  dist_mat <- dist(seurat_obj@reductions$pca@cell.embeddings,
                   method = "euclidean")
  
  cluster_info <- seurat_obj@meta.data[, grepl(paste0(assay,"_snn_res"),
                                               colnames(seurat_obj@meta.data))] %>%
    mutate_all(as.numeric)
  
  res_scores <- apply(cluster_info, MARGIN = 2, FUN = function(res){
    sil_val <- silhouette(res, dist_mat)
    if (!any(is.na(sil_val))){
      mean(sil_val[, "sil_width"])
    }
    else {
      NA
    }
  })
  
  seurat_obj[["opt_clust"]] <- seurat_obj[[names(which.max(res_scores))]]
  
  # Remove sub-optimal cluster columns from metadata
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    select(!contains("_snn_res"))
  
  return(seurat_obj)
}