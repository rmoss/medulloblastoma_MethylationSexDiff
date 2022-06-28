#!/usr/bin/Rscript
# Run R script to load GRSet object and pheno file and then run bumphunter to find differentially methylated regions
# Rscript --vanilla $scriptsDir/DMRbySex.R -g $grset.rds -f $filename -p $prefix -n $number -o $outDir

# Load required libraries
library(optparse)
library(limma)
library(minfi)
library(tidyverse)
library(doParallel)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Details of list of options for the script
option_list = list(
make_option(c("-g", "--grset"), type = "character", default = NULL, help = "path/file", metavar = "character"),
make_option(c("-f", "--filename"),type = "character", default = NULL, help = "path/pheno_file"),
make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "name for output file", metavar = "character"),
make_option(c("-n", "--number"), type="integer", default = NULL, help = "number of iterations"),
make_option(c("-o", "--outDir"), type = "character", default = NULL, help ="write final SNP density profiles to outDir ", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$outDir)

# cutoff for bumphunter
cut <- 0.05
# number of permutations
num <- opt$number

#' Title
#'
#' @param GR a genomic ratio set methylation object created by minfi
#' @param design design or model matrix for the analysis
#' @param cutoff a numeric value that sets the cutoff of % difference in beta values 
#' @param perm integer denoting the number of permutations or resamples to use to compute null distributions
#'
#' @return
#' @export
#'
#' @examples
findDMRs <- function(GR, design, cutoff, perm, dir, pre, pheno){
  dmrs <- bumphunter(GR, design = design, coef=2,
                    cutoff = cutoff, B=perm, type="Beta", nullMethod=c("bootstrap"))
  saveRDS(dmrs, paste0(dir, "/", pre,"_DMRs_sexdiff.rds", sep=""))

  dmrs_table = dmrs$table
  dmrs_table <- as.data.frame(dmrs_table)
  write.table(dmrs_table, paste0(dir, "/", pre,"_sexdmrs_raw.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  dmrs_table$strand <- "."
  dmrs_table <- unite(dmrs_table,name, chr, start, end, sep = "_", remove = F)

  # Write to bed file for processing with bedtools closest
  write.table(filter(dmrs_table, fwer < 0.05)[,c("chr","start","end","name","fwer","strand")], file = paste0(dir, "/", pre, "_sexdmrs_sigRegions.bed"), sep = "\t",quote = F, col.names = F, row.names = F)

  # beta <- getBeta(GR)
  # # filtering DMP’s to significant ones (qvalue < 0.05)
  # dmrSig <- dmrs_table[dmrs_table$fwer < .05,]
  # sigBeta <- beta[rownames(dmrSig), ]
  # sigBeta <- as.data.frame(sigBeta)
  # sigBeta <- rownames_to_column(sigBeta, var = "probeID")
  # dim(dmrSig)

  # # finding delta beta values (Male - Female)
  # betaLong2 <- pivot_longer(sigBeta, -probeID, names_to = "sampleID", values_to = "betaVal")
  # # pheno2 <- as.data.frame(pData(GRSet)) 
  # # pheno2 <- rownames_to_column(pheno2, var = "sampleID")
  # betaLong2 <- left_join(betaLong2, pheno) 
  # betaLong2 <- separate(betaLong2, "sampleID", "GSMid", sep = "_", remove = F) ### "sampleID" has messy characters that need to be removed

  # # creating a table (called “summary2”) showing the average beta value for males and females on each probe
  # summary2 <- betaLong2 %>% 
  #   group_by(probeID, sex) %>%   
  #   dplyr::summarize(avgBeta = mean(betaVal))
  # summaryWide <- pivot_wider(summary2[,c("probeID","sex","avgBeta")], names_from = sex, values_from = avgBeta)
  # diffValue <- summaryWide
  # diffValue$diff <- summaryWide$`M`-summaryWide$`F`
  # head(diffValue)

  # write_csv(diffValue, paste0(dir,"/",pre,"_DeltaBetapredSex.csv"))

  # # finding negative or positive differences (Male - Female)
  # neg_diff <- diffValue %>%
  #   filter(diff < 0)
  # nrow(neg_diff)
  # pos_diff <- diffValue %>%
  #   filter(diff > 0)
  # nrow(pos_diff)

  # # Creating Full Results table (with probeID, gene name, p-value, delta beta):
  # dmpSig <- dmpSig %>%
  #   rownames_to_column(var = "probeID")
  # full_predSex <- full_join(dmpSig, diffValue, by = "probeID")
  # ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # ann450k_df <- as.data.frame(ann450k)
  # ann450kNew <- filter(ann450k_df, Name %in% full_predSex$probeID)
  # ann450kNew <- ann450kNew %>%
  #   rownames_to_column(var = "probeID")
  # gene_probeID <- ann450kNew %>%
  #   dplyr::select(probeID, UCSC_RefGene_Name)
  # full_predSex <- full_predSex %>%
  #   left_join(dplyr::select(gene_probeID, probeID, UCSC_RefGene_Name), by = "probeID")
  # full_predSex <- full_predSex %>%
  #   filter(!UCSC_RefGene_Name == "")

  # write_csv(full_predSex, paste0(dir,"/",pre,"_MethFullResults.csv"))
  # # Return the DMRs with full results including annotation, p-value, and delta beta
  return(dmrs)
}

GRSet <- readRDS(opt$grset)
pheno <- read_csv(opt$filename)

keep <- !is.na(pheno$age)
pheno_filt <- pheno[keep,]
GRSet_filt <- GRSet[,pheno_filt$geo_accession]

pheno_sex <- factor(pheno_filt$sex)
pheno_age <- pheno_filt$age
pheno_agecat <- factor(pheno_filt$agecat)

registerDoParallel(cores = 8)

# Run bumphunter to find DMRs by sex with age adjustment
designMatrix1 <- model.matrix(~ pheno_sex + pheno_age)
dmrs_all_age <- findDMRs(GRSet_filt, designMatrix1, cut, num, opt$outDir, paste0(opt$prefix,"_ageadj"), pheno_filt)

# Do the same but with age as a category
designMatrix2 <- model.matrix(~ pheno_sex + pheno_agecat)
dmrs_all_agecat <- findDMRs(GRSet_filt, designMatrix2, cut, num, opt$outDir, paste0(opt$prefix,"_agecatadj"), pheno_filt)

# Break dataset into subgroups and analyze for DMRs within each group

subgroups <- c("SHH","WNT","Group3","Group4")

lsDMRsSubAge <- list()
lsDMRsSubAgeCat <- list()
for (sub in subgroups){
  pheno_sub <- filter(pheno_filt, subgroup==sub)
  GRSet_sub <- GRSet_filt[,pheno_sub$geo_accession]
  prefix <- paste0(opt$prefix, "_",sub)
  pheno_sex <- pheno_sub$sex
  pheno_age <- pheno_sub$age
  pheno_agecat <- pheno_sub$agecat
  designMatrix1 <- model.matrix(~ pheno_sex + pheno_age)
  lsDMRsSubAge[[sub]] <- findDMRs(GRSet_sub, designMatrix1, cut, num, opt$outDir, paste0(prefix,"_ageadj"), pheno_sub)
  designMatrix2 <- model.matrix(~ pheno_sex + pheno_agecat)
  lsDMRsSubAgeCat[[sub]] <- findDMRs(GRSet_sub, designMatrix2, cut, num, opt$outDir, paste0(prefix,"_agecatadj"), pheno_sub)
}




