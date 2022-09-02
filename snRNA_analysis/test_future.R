library(Seurat)
pbmc <- readRDS(file = "/home/dwzj28/research_project/src/pbmc3k_final.rds")

timing.comparisons <- data.frame(fxn = character(), time = numeric(), strategy = character())
plan("sequential")
start <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "ScaleData", time = as.numeric(end -
                                                                                                  start, units = "secs"), strategy = "sequential"))

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "FindMarkers", time = as.numeric(end -
                                                                                                    start, units = "secs"), strategy = "sequential"))

plan("multisession", workers = 4)
start <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "ScaleData", time = as.numeric(end -
                                                                                                  start, units = "secs"), strategy = "multiprocess"))

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "FindMarkers", time = as.numeric(end -
                                                                                                    start, units = "secs"), strategy = "multiprocess"))
print(timing.comparisons)
