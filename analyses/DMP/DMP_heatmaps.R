#!/usr/bin/Rscript
# Run R script to load GRSet object and pheno file and then run bumphunter to find differentially methylated regions
# Rscript --vanilla $scriptsDir/DMP_heatmaps.R -g $datadir/grsetfile.rdfs -f $datadir/$phenofilename -p $cohortprefix -d $datadir -o $outDir

# Load required libraries
library(optparse)
library(limma)
library(minfi)
library(tidyverse)
library(doParallel)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(siggenes)
library(ComplexHeatmap)
require(RColorBrewer)
require(circlize)

# Details of list of options for the script
option_list = list(
make_option(c("-g", "--grset"), type = "character", default = NULL, help = "path/file", metavar = "character"),
make_option(c("-f", "--filename"),type = "character", default = NULL, help = "path/pheno_file"),
make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "name for output file", metavar = "character"),
make_option(c("-d", "--dir"), type = "character", default = NULL, help = "dir for dmp data file",),
make_option(c("-o", "--outDir"), type = "character", default = NULL, help ="write output to outDir ", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd(opt$outDir)

GRSet <- readRDS(opt$grset)
pheno <- read_csv(opt$filename)
keep <- !is.na(pheno$age)
pheno_filt <- pheno[keep,]
GRSet_filt <- GRSet[,pheno_filt$geo_accession]

if (opt$prefix=="Cavalli") {
  pheno_filt$sex <- factor(pheno_filt$sex, levels=c("M","F"))
  pheno_filt$Histology <- ifelse(is.na(pheno_filt$Histology),"Unknown",pheno_filt$Histology)
  pheno_filt$agecat <- factor(pheno_filt$agecat)
  pheno_filt$agecat <- recode_factor(pheno_filt$agecat, "1"="<5", "2"="5-9","3"="10-14","4"="15-19","5"=">=20")
  pheno_filt$`Metastatic Disease` <- ifelse(is.na(pheno_filt$`Metastatic Disease`), "Unknown", pheno_filt$`Metastatic Disease`)
  pheno_filt$`Metastatic Disease` <- recode_factor(pheno_filt$`Metastatic Disease`, "1"="Yes", "0"="No", "NA"="Unknown")
  # pheno_filt$`Vital Status` <- factor(pheno_filt$`Vital Status`)
  pheno_filt$`Vital Status` <- ifelse(is.na(pheno_filt$`Vital Status`), "Unknown", pheno_filt$`Vital Status`)
  pheno_filt$`Vital Status` <- recode_factor(pheno_filt$`Vital Status`, "1"="Deceased", "0"="Alive")
} else {
  pheno_filt$sex <- factor(pheno_filt$sex, levels=c("M","F"))
  pheno_filt$pathology <- ifelse(is.na(pheno_filt$pathology),"Unknown",pheno_filt$pathology)
  colnames(pheno_filt)[8] <- "Histology"
  pheno_filt$Histology <- recode_factor(pheno_filt$Histology, "CLA"="Classic","DN"="Desmoplastic","LCA"="LCA","MBEN"="MBEN","NOS"="MB Not Otherwise Specified")
  pheno_filt$agecat <- factor(pheno_filt$agecat)
  pheno_filt$agecat <- recode_factor(pheno_filt$agecat, "1"="<5", "2"="5-9","3"="10-14","4"="15-19","5"=">=20")
}

subgroups <- c("SHH","WNT","Group3","Group4")

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
ann450k_df <- as.data.frame(ann450k)
# load(paste0(opt$dir,"/../","chromColors.RDS"))

# myCol <- colorRampPalette(c('dodgerblue', 'white', 'red'))(100)
# myBreaks <- seq(0, 1, length.out = 100)

# lsDMPsSubAgeCat <- list()
for (sub in subgroups){
  prefix <- paste0(opt$prefix, "_",sub)
  dmpSig <- read_csv(paste0(opt$dir,"/",prefix,"_ageacatadj_dmp_MethFullResults.csv"))
  dmpSig <- column_to_rownames(dmpSig, var="probeID")
  pheno_sub <- filter(pheno_filt, subgroup==sub)
  GRSet_sub <- GRSet_filt[,pheno_sub$geo_accession]
  beta <- getBeta(GRSet_sub)
  sigBeta <- beta[rownames(dmpSig), ] ### only including significant beta values of significant DMP's 
  sigBeta <- as.data.frame(sigBeta) ### to make a "probeID" column (helpful in the future), changing our beta value matrix into a dataframe
  sigBeta <- rownames_to_column(sigBeta, var = "probeID")
  # write_csv(sigBeta, paste0(opt$outDir,"/tables/",prefix,"_betavalues_20221013.csv"))

  betaLong2 <- pivot_longer(sigBeta, -probeID, names_to = "sampleID", values_to = "betaVal") ###
  pheno2 <- as.data.frame(pData(GRSet_sub)) ### creating pheno data from Group GRset
  pheno2 <- rownames_to_column(pheno2, var = "sampleID")
  betaLong2 <- left_join(betaLong2, pheno2) ### joining pData to beta Data by "sampleID"
  betaLong2 <- separate(betaLong2, "sampleID", "GSMid", sep = "_", remove = F) ### "sampleID" has messy characters 
  summary2 <- betaLong2 %>% 
      group_by(probeID, predictedSex) %>%   
      dplyr::summarize(avgBeta = mean(betaVal))
  ### converting the summary table to a wider format
  summaryWide <- pivot_wider(summary2[,c("probeID","predictedSex","avgBeta")], names_from = predictedSex, values_from = avgBeta)
  diffValueG <- summaryWide
  diffValueG$diff <- summaryWide$`M`-summaryWide$`F` ### calculating delta beta values M-F
  head(diffValueG)
  neg_diff_G <- diffValueG %>%
      filter(diff < 0)
  nrow(neg_diff_G)
  pos_diff_G <- diffValueG %>%
      filter(diff > 0)
  nrow(pos_diff_G)
  dmpSig <- dmpSig %>%
    rownames_to_column(var = "probeID")
  full_GpredSex <- full_join(diffValueG, dmpSig, by = "probeID")
  # full_GpredSex_genes <- full_GpredSex %>% rowwise() %>%
  #   mutate(UCSC_RefGene_Name = toString(unique(unlist(strsplit(UCSC_RefGene_Name,";")))))
  full_GpredSex_genes <- full_GpredSex %>% separate_rows(UCSC_RefGene_Name) %>% distinct()
  genes <- as.data.frame(full_GpredSex_genes$UCSC_RefGene_Name)
  colnames(genes) <- c(paste0(opt$prefix,"_",sub,"_DMP_genes"))
  colnames(full_GpredSex) <- c("probeID","Mean Female beta","Mean Male beta","diff Beta", "logFC","AveExpr","t","p.value","adj.p.value", "B","UCSC_RefGene_Name")
  full_GpredSex <- full_GpredSex[,c(1,5:10,2,3,4,11)]
  write_csv(full_GpredSex, paste0(opt$outDir,"/tables/",prefix,"_DMPFullResults_v2.csv"))
  # write_csv(genes, paste0(opt$outDir,"/tables/",prefix,"_DMP_genelist.csv"))

#   ann450kG <- filter(ann450k_df, Name %in% full_GpredSex$probeID)
#   ann450kG <- ann450kG %>%
#     rownames_to_column(var = "probeID")
  
#   sex <- pheno_sub$sex
#   #sex <- as.data.frame(sex)
#   # age <- pheno_sub$age
#   agecat <- pheno_sub$agecat
#   #agecat <- as.data.frame(agecat)
#   histology <- pheno_sub$Histology
#   #histology <- as.data.frame(histology)
#   if (opt$prefix=="Cavalli") {
#     met <- pheno_sub$`Metastatic Disease`
#     #met <- as.data.frame(met)
#     dead <- pheno_sub$`Vital Status`
#     #dead <- as.data.frame(dead)
    
#     ann1 <- data.frame(
#       histology = histology,
#       dead = dead,
#       met_status = met,
#       age = agecat,
#       sex = sex)
#     colours <- list(
#       histology = c("Classic" = "deepskyblue", "Desmoplastic" = "royalblue1", "LCA" = "blue3", "MBEN" = "cyan", "Unknown" = "grey91"),
#       dead = c("Alive" = "palegreen3", "Deceased" = "darkgreen", "Unknown" = "grey91"),
#       met_status = c("No" = "royalblue1", "Yes" = "navyblue", "Unknown" = "grey91") ,
#       age = c("<5" = "lightyellow1", "5-9" = "yellow2", "10-14" = "gold2", "15-19" = "goldenrod2", ">=20" = "darkgoldenrod3", "Unknown" = "grey91"),
#       sex = c("M" = "black", "F"="grey"))
#     label <- c("Histology","Vital Status","Metastatic Disease","Age at Diagnosis","Sex")
#   } else {
#     ann1 <- data.frame(
#       histology = histology,
#       age = agecat,
#       sex = sex)
#     colours <- list(
#       histology = c("Classic" = "deepskyblue", "Desmoplastic" = "royalblue1", "LCA" = "blue3", "MBEN" = "cyan", "Unknown" = "grey91", "MB Not Otherwise Specified" = "grey91"),
#       age = c("<5" = "lightyellow1", "5-9" = "yellow2", "10-14" = "gold2", "15-19" = "goldenrod2", ">=20" = "darkgoldenrod3", "Unknown" = "grey91"),
#       sex = c("M" = "black", "F"="grey"))
#     label <- c("Histology","Age at Diagnosis","Sex")
#   }
#   colAnn_G <- HeatmapAnnotation(df = ann1,
#       which = "col",
#       col = colours,
#       annotation_height = 0.6,
#       annotation_width = unit(1, 'cm'),
#       show_legend = TRUE,
#       annotation_label = label,
#       gap = unit(0, 'mm'))    
#   # ann450kG <- filter(ann450k_df, Name %in% full_GpredSex$probeID)
#   # ann450kG <- ann450kG %>%
#   #   rownames_to_column(var = "probeID")
#   # chr <- ann450kG$chr
#   # armytage_test2 <- list(chr = armytage)
#   # chr2 <- as.data.frame(chr)
#   DMP_diff <- full_GpredSex$diff

#   # ann2 <- data.frame(
#   #   chr = chr,
#   #   DMP_diff = DMP_diff   
#   # )  
#   # colours2 <- list(
#   #   chr = armytage,
#   #   DMP_diff = colorRamp2(c(-0.2, 0, 0.2), c("purple4", "white", "mediumorchid1"))
#   # )
#   #col_diff <- colorRampPalette(c('#7fbf7b', '#f7f7f7', '#af8dc3'))(100)
#   # myBreaks <- seq(-0.2, 0, 0.2, length.out = 100)
#   myBreaks <- c(-0.2,0,0.2)
#   row_fun <- colorRamp2(myBreaks, c("mediumorchid1", "white", "purple4"))
#   rowAnn_G_diff <- HeatmapAnnotation(DMP_diff = DMP_diff, 
#                                        which = "row",
#                                        simple_anno_size = unit(1,"cm"),
#                                        col = list(DMP_diff = row_fun),
#                                        show_legend = TRUE,
#                                        show_annotation_name = FALSE,
#                                        annotation_label="Meth. Diff.")
  
#   # rowAnn_G2 <- HeatmapAnnotation(df = chr2, 
#   #                               which = "row",
#   #                               col = armytage_test2,
#   #                               annotation_width = unit(1.5, "cm"),
#   #                               show_legend = TRUE,
#   #                               show_annotation_name = FALSE)
  
#   # rowAnn_G <- HeatmapAnnotation(df = ann2, 
#   #                               which = "row",
#   #                               col = colours2,
#   #                               annotation_width = unit(1.5, "cm"),
#   #                               show_legend = TRUE,
#   #                               show_annotation_name = FALSE,
#   #                               annotation_label=c("","Meth. Diff."))

#   sigBeta <- sigBeta %>% column_to_rownames(var="probeID")
#   sigBeta <- as.matrix(sigBeta)
#   sigBeta_scaled_mat = t(scale(t(sigBeta)))
#   # ht <- Heatmap(sigBeta, 
#   #     name = sub,
#   #     cluster_rows = TRUE,
#   #     show_row_dend = TRUE, 
#   #     row_title = "DMP's (q-value < 0.05)",
#   #     row_title_side = 'left',
#   #     row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
#   #     row_title_rot = 90,
#   #     show_row_names = FALSE,
#   #     row_dend_width = unit(25,'mm'),
#   #     cluster_columns = TRUE,
#   #     show_column_dend = TRUE,
#   #     column_title = '',
#   #     column_title_side = 'bottom',
#   #     column_title_gp = gpar(fontsize = 12, fontface
#   #                                             ='bold'),
#   #     column_title_rot = 0,
#   #     show_column_names = FALSE,
#   #     show_heatmap_legend = TRUE,
#   #     top_annotation = colAnn_G,
#   #     left_annotation = c(rowAnn_G_diff,rowAnn_G2)
#   # )

#   ht1 <- Heatmap(sigBeta, 
#     name = sub,
#     width = 10,
#     height = 8,
#     cluster_rows = TRUE,
#     show_row_dend = TRUE, 
#     row_title = "DMP's (q-value < 0.05)",
#     row_title_side = 'left',
#     row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
#     row_title_rot = 90,
#     show_row_names = FALSE,
#     # row_dend_width = unit(25,'mm'),
#     cluster_columns = TRUE,
#     show_column_dend = TRUE,
#     column_title = '',
#     column_title_side = 'bottom',
#     column_title_gp = gpar(fontsize = 12, fontface
#                                             ='bold'),
#     column_title_rot = 0,
#     show_column_names = FALSE,
#     show_heatmap_legend = TRUE,
#     top_annotation = colAnn_G #,
#     #left_annotation = rowAnn_G_diff
#   )

#   ht2 <- Heatmap(sigBeta_scaled_mat, 
#     name = sub,
#     width = 10,
#     height = 8,
#     cluster_rows = TRUE,
#     show_row_dend = TRUE, 
#     row_title = "DMP's (q-value < 0.05)",
#     row_title_side = 'left',
#     row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
#     row_title_rot = 90,
#     show_row_names = FALSE,
#     # row_dend_width = unit(25,'mm'),
#     cluster_columns = TRUE,
#     show_column_dend = TRUE,
#     column_title = '',
#     column_title_side = 'bottom',
#     column_title_gp = gpar(fontsize = 12, fontface
#                                             ='bold'),
#     column_title_rot = 0,
#     show_column_names = FALSE,
#     show_heatmap_legend = TRUE,
#     top_annotation = colAnn_G #,
#     #left_annotation = rowAnn_G_diff
#   )
#   # draw the heatmap
#   pdf(file = paste0(opt$outDir,"/heatmaps/",prefix,"_DMPHeatmap_v2.pdf"), width = 7, height = 7)
#   heatDraw <- draw(ht1,
#     heatmap_legend_side = 'right',
#     annotation_legend_side = 'right',
#     row_sub_title_side = 'left')
#   dev.off()
  
#   pdf(file = paste0(opt$outDir,"/heatmaps/",prefix,"_DMPHeatmap_scaled_v2.pdf"), width = 7, height = 7)
#   heatDraw <- draw(ht2,
#     heatmap_legend_side = 'right',
#     annotation_legend_side = 'right',
#     row_sub_title_side = 'left')
#   dev.off()
}

