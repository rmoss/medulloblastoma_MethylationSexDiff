#!/usr/bin/Rscript

RGSet <- readRDS("RGset.rds")
colnames(RGSet) <- gsub("_.*","",colnames(RGSet))
dir <- "/home/spectorl/moss0134/projects/research/MBSexDiff/data/Cavalli/"
pheno <- read_csv(paste0(dir, "phenoCav.csv"))
RGSet_filt <- RGSet[,pheno$geo_accession]

# Create methylset
Mset <- preprocessRaw(RGSet_filt)
saveRDS(Mset, paste0(dir,"/MSet_Cavalli_filt.RDS"))