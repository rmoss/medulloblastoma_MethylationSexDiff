#!/usr/bin/Rscript

### Load packages
library(minfi)
library(conumee)
library(GEOquery)
library(tidyverse)
library(readr)
library(filesstrings)
library(gridExtra)
library(optparse)

# Details of list of options for the script
option_list = list(
make_option(c("-a", "--mset"), type = "character", default = NULL, help = "path/file of methylset", metavar = "character"),
make_option(c("-b", "--phenofile"),type = "character", default = NULL, help = "path/pheno_file"),
make_option(c("-c", "--ctrlmset"),type = "character", default = NULL, help = "path/ctrl_methylset"),
make_option(c("-d", "--ctrlphenofile"),type = "character", default = NULL, help = "path/ctrl_phenofile"),
make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "name for output file", metavar = "character"),
make_option(c("-s", "--subgroup"), type="character", default = NULL, help = "subgroup to work with"),
make_option(c("-o", "--outDir"), type = "character", default = NULL, help ="write final plots and tables to outDir ", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$outDir)

### Specify files and paths
### Define functions

### Reading data

# We'll first read in the data for samples we are processing and control data as methylsets. I'll also read in the pheno data for the samples that contains subgroups

MSet <- readRDS(opt$mset)
pheno <- read_csv(opt$phenofile)
Mset_Ctrl <- readRDS(opt$ctrlmset)
phenoCtrl <- read_csv(opt$ctrlphenofile)
phenoCtrl <- phenoCtrl[match(colnames(Mset_Ctrl),phenoCtrl$Basename),]

## Use Conumee to call CNVs for Newcastle
# Now that we have done all of our setup, we can move on to the CNV analysis using conumee

### Quick stats
# Describe the data we have.

cat("Sample data contains:", ncol(MSet),"medulloblastoma samples")
cat("\n") # print return between
cat("Counts of samples by subgroup: ")
table(pheno$subgroup)
cat("Control data contains:", ncol(Mset_Ctrl), "samples")


### Create annotation object

data(exclude_regions)
data(detail_regions)
head(detail_regions, n = 2)
anno <- CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions)
anno

# Break into subgroup
pheno_sub <- filter(pheno, subgroup==opt$subgroup)
MSet_sub <- Mset[,pheno_sub$geo_accession]

### Combine intensity values
sample.data <- CNV.load(MSet_sub)
controls.data <- CNV.load(Mset_Ctrl)
names(sample.data)
sample.data



### Now loop through all the samples and save plot and table

datalist <- list()

for (samp in names(sample.data)){
  x <- CNV.fit(sample.data[samp], controls.data, anno)
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
  CNV.genomeplot(x)
  plist[[samp]] <- CNV.genomeplot(x)
  pdf(file=paste0(opt$outDir,"/",opt$subgroup,"_" sprintf("p%s.pdf", samp)),
      width = 6, height = 4, onefile = T)
  CNV.genomeplot(x)
  dev.off()
  dat <- CNV.write(x, what = "segments")
  dat$samp <- samp
  datalist[[samp]] <- data
}

big_data <- do.call(rbind, datalist)
write_csv(big_data, file=paste0(opt$outDir,"/",opt$prefix,"_",opt$subgroup,"_CNVsegments.csv"))


## Session information
sessionInfo()