library(ggplot2)
library(read.xlsx)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(tidyverse)
library(readxl)
library(Biobase)
library(GEOquery)
library(sva)
library(devtools)
library(ComplexHeatmap)
library(doParallel)
library(gtools)
library(Polychrome)
library(openxlsx)

# Run detP function on RGSet & normalize using preprocessFunnorm
detP <- detectionP(RGSet_filt)

# Check predicted sex and remove samples that that don't match reported sex     
predictedSex <- getSex(GRSet, cutoff = -2)$predictedSex
head(predictedSex)
GRSet <- addSex(GRSet) ### adding predictedSex to GRSet
#png(file="./Results_Newcastle/MBMeth_Newcastle_predictedsexPlot.png")
plotSex(GRSet)
#dev.off()

phenoData <- as.data.frame(pData(GRSet))
keep <- filter(phenoData, ifelse(gender.ch1 == predictedSex,T,F))$sampleName
keep
GRSet.flt <- GRSet[,keep]

# removing SNP’s from GRSet
snps <- getSnpInfo(GRSet.flt)
GRSet.flt <- addSnpInfo(GRSet.flt)
GRSet.flt <- dropLociWithSnps(GRSet.flt, snps=c("SBE","CpG"), maf=0)

# removing cross reactive probes
xReactiveProbes <- read.xlsx("../../HB_Methylation_Project/HBMeth_CombinedAnalysis/Data3/48639-non-specific-probes-Illumina450k.xlsx", 1)
xReactiveProbes <- xReactiveProbes[match(featureNames(GRSet.flt),xReactiveProbes$TargetID),]
keep_nonxreact <- !(featureNames(GRSet.flt) %in% xReactiveProbes$TargetID)
table(keep_nonxreact)
GRSet.flt <- GRSet.flt[keep_nonxreact,]
GRSet.flt

# Removing Sex chromosomes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(GRSet.flt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
GRSet.flt <- GRSet.flt[keep, ]

# Removing probes with high p-values
detP <- detP[match(featureNames(GRSet.flt),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(GRSet)
table(keep)
# keep
# FALSE   TRUE 
# 67510 362434 

GRSet.flt <- GRSet.flt[keep,]

# saveRDS(GRSet.flt, "GRSet_GSE93646_filt.RDS")

#### Start here if loading filtered GRSet ####
setwd("/Volumes/GoogleDrive/My Drive/Research/SexDiffMedullo/Analysis_RMM_new")
# GRSet.flt <- readRDS("GRSet_GSE93646_filt.RDS")

# phenoData <- as.data.frame(pData(GRSet.flt))
# samples <- phenoData[,c("sampleName","predictedSex")]
# samples <- samples %>%
#   separate(sampleName, c("accession","Sentrix_ID","Sentrix_Position"), sep="_") %>%
#   rownames_to_column(var="Sample_Name") %>%
#   add_column(Sample_Plate=NA) %>%
#   add_column(Sample_Well=NA) %>%
#   add_column(Pool_ID=NA) %>%
#   rename(Sample_Group=predictedSex) %>%
#   select(Sample_Name,Sample_Plate,Sample_Group,Sample_Well,Pool_ID,Sentrix_ID,Sentrix_Position)

# write_csv(phenoData, file="./data/Newcastle_pheno_filt.csv")
# write_csv(samples,file="./data/sampleSheet.csv")


# getting beta Values from new filtered GRSet
# beta_matrix <- getBeta(GRSet.flt)
# beta <- as.data.frame(beta_matrix)
#write_csv(beta, file="./data/Newcastle_beta_filt.csv")

# Read in beta values from filtered GRSet.flt
beta <-read_csv("./data/Newcastle_beta_filt.csv")
phenoData <- read_csv("./data/Newcastle_pheno_filt.csv")

# getting DMP’s from beta values
sex <- phenoData$predictedSex
dmp <- dmpFinder(beta, pheno = sex, type = "categorical")

# filtering DMP’s to significant ones (qvalue < 0.05)
dmpSig <- dmp[as.numeric(as.character(dmp$qval))<0.05,]
sigBeta <- beta[rownames(dmpSig), ]
sigBeta <- as.data.frame(sigBeta)
sigBeta <- rownames_to_column(sigBeta, var = "probeID")
dim(dmpSig)

# finding delta beta values (Male - Female)
betaLong2 <- pivot_longer(sigBeta, -probeID, names_to = "sampleID", values_to = "betaVal")
pheno2 <- as.data.frame(pData(GRSet.flt)) 
pheno2 <- rownames_to_column(pheno2, var = "sampleID")
betaLong2 <- left_join(betaLong2, pheno2) 
betaLong2 <- separate(betaLong2, "sampleID", "GSMid", sep = "_", remove = F) ### "sampleID" has messy characters that need to be removed

# creating a table (called “summary2”) showing the average beta value for males and females on each probe
summary2 <- betaLong2 %>% 
  group_by(probeID, predictedSex) %>%   
  dplyr::summarize(avgBeta = mean(betaVal))
summaryWide <- pivot_wider(summary2[,c("probeID","predictedSex","avgBeta")], names_from = predictedSex, values_from = avgBeta)
diffValue <- summaryWide
diffValue$diff <- summaryWide$`M`-summaryWide$`F`
head(diffValue)

# write_csv(diffValue, "./Results_Newcastle/DeltaBetapredSex_v2.csv")

# finding negative or positive differences (Male - Female)
neg_diff <- diffValue %>%
  filter(diff < 0)

nrow(neg_diff)

pos_diff <- diffValue %>%
  filter(diff > 0)

nrow(pos_diff)

# Creating Full Results table (with probeID, gene name, p-value, delta beta):
dmpSig <- dmpSig %>%
  rownames_to_column(var = "probeID")
full_predSex <- full_join(dmpSig, diffValue, by = "probeID")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k_df <- as.data.frame(ann450k)
ann450kNew <- filter(ann450k_df, Name %in% full_predSex$probeID)
ann450kNew <- ann450kNew %>%
  rownames_to_column(var = "probeID")
gene_probeID <- ann450kNew %>%
  dplyr::select(probeID, UCSC_RefGene_Name)
full_predSex <- full_predSex %>%
  left_join(dplyr::select(gene_probeID, probeID, UCSC_RefGene_Name), by = "probeID")
full_predSex <- full_predSex %>%
  filter(!UCSC_RefGene_Name == "")

# write_csv(full_predSex, "./Results_Newcastle/Newcastle_MethFullResults_v2.csv")

#### I stopped here ####

# Making Heatmaps
# creating column annotation of sex
sex <- as.data.frame(sex)
col_fun = list(sex = c("M" = "black", "F"="grey"))

colAnn <- HeatmapAnnotation(df = sex, 
                                which = "col",
                                col = col_fun,
                                annotation_height = 0.6,
                                annotation_width = unit(1, 'cm'),
                                show_legend = FALSE,
                                gap = unit(1, 'mm'))
# creating row annotation of chromosome info
# First, getting chromosome information:
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k_df <- as.data.frame(ann450k)
ann450k <- filter(ann450k_df, Name %in% sigBeta$probeID)
chr2 <- unique(ann450k$chr)
armytage <- as.vector(unlist(green.armytage.colors(20)))
armytage_test2 <- list(armytage = c(chr2 = armytage))
chr2 <- as.data.frame(chr2)
rowAnn <- HeatmapAnnotation(df = chr2, 
                                which = "row",
                                col = armytage_test2,
                                annotation_width = unit(1.5, "cm"),
                                show_legend = FALSE,
                                show_annotation_name = FALSE)

# Need to re-run these two lines right before the heatmap so that sigBeta can properly be put back into a matrix:
dmpSig <- dmp[as.numeric(as.character(dmp$qval))<0.05,]
sigBeta <- beta[rownames(dmpSig), ]

# Running Heatmap
sigBeta_matrix <- as.matrix(sigBeta)

Heatmap(sigBeta_matrix, 
        name = "Newcastle w/o chr",
        cluster_rows = TRUE,
        show_row_dend = TRUE, 
        row_title = "DMP's (q-value < 0.05)",
        row_title_side = 'left',
        row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
        row_title_rot = 90,
        show_row_names = FALSE,
        row_dend_width = unit(25,'mm'),
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        column_title = '',
        column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 12, fontface
                               ='bold'),
        column_title_rot = 0,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE,
        top_annotation = colAnn,
        left_annotation = rowAnn
)

