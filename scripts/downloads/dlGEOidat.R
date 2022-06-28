library(minfi)
library(GEOquery)

#increase file download timeout
options(timeout = 600)

#download GEO object
gse <- getGEO("GSE93646", GSEMatrix = TRUE)
#get phenotype data - sample sheet
pd = pData(gse[[1]])

#get raw data - idats, processed beta matrix, etc.
getGEOSuppFiles("GSE93646")

#decompress idats
untar("GSE93646/GSE93646_RAW.tar", exdir = "GSE93646/idat")
#list files
head(list.files("GSE93646/idat", pattern = "idat"))
idatFiles <- list.files("GSE93646/idat", pattern = "idat.gz$", full = TRUE)
#decompress individual idat files
sapply(idatFiles, gunzip, overwrite = TRUE)
#read idats and create RGSet
RGSet <- read.metharray.exp("GSE93646/idat")

saveRDS(RGSet, "RGSet_GSE93646.RDS")