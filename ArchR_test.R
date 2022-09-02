library(parallel)
library(ArchR)
path="/home/dwzj28/research_project/dataset/snATAC/10X_ATAC_CK166.tsv.gz"
base=sub("[^/]*$", "", path)
input_file=c(CK166=path)
addArchRThreads(threads=parallel::detectCores()-2)
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(input_file,
			       sampleNames = names(input_file),
			       outputNames = names(input_file),
			       minTSS = 4,
			       minFrags = 3000,
			       QCDir = paste0(base,"QualityControl"),
			       addTileMat = TRUE,
			       addGeneScoreMat = TRUE)
