library(rhdf5)

path = "/home/dwzj28/research_project/dataset/snRNA/sample_data/CK158/outs/snRNA_CK158_molecule_info.h5"

#Print the group structure of the H5DF file 
h5dump(path, load=FALSE)

#Open the H5DF file at path
obj <- H5Fopen(path)
paste("Alignment genome:", obj$"barcode_info/genomes")
#GEM group is the number appended after the barcode e.g. AG...TC-1
paste("Library information:", obj$library_info)
#Confirms whitelist as 3M-february-2018
paste("JSON metrics:", obj$metrics_json)
#print if barcode_info$pass_filter contains values other than zero
any(obj$barcode_info$pass_filter[2:3,]!=0)
#Max barcode value: 2588190
max(obj$barcode_info$pass_filter[1,])
#Even distribution of barcode labels across entire index -> no bias 
hist(obj$barcode_info$pass_filter[1,], breaks=30)
#Total number of present barcodes
length(unique(obj$barcode_idx))
#max count at 910
max(obj$count)
#zero-based indexing corresponding to TNNT2 - cardiac troponin (makes sense)
obj$"features/name"[obj$feature_idx[which.max(obj$count)] + 1]
#Only Gene Expression features counted
unique(obj$"features/feature_type")

h5closeAll()
gc()

