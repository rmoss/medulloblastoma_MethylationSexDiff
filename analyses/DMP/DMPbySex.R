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
library(siggenes)

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

dmpFinder_custom <- function(dat, pheno, cov, pre, dir) {
    design <- model.matrix(~0 + pheno + cov)
    fit <- lmFit(dat, design)
    contrasts <- makeContrasts(phenoM-phenoF, levels=colnames(design))
    fitC <- contrasts.fit(fit,contrasts)
    fitC <- eBayes(fitC)
    dt <- decideTests(fitC)
    tab <- topTable(fitC, n=500000, p.value=1.0)
    png(file=paste0(dir,"/",pre,"_sexDiffDMP.png"))
    plotMD(fitC, column=1, status=dt[,1], main=colnames(fitC)[1], xlim=c(-0.5,2))
    dev.off()
    # Creating Full Results table (with probeID, gene name, p-value, delta beta):
    dmpSig <- tab[tab$adj.P.Val < .05,] %>% 
        rownames_to_column(var="probeID")
    ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann450k_df <- as.data.frame(ann450k)
    ann450kNew <- filter(ann450k_df, Name %in% dmpSig$probeID)
    ann450kNew <- ann450kNew %>%
        rownames_to_column(var = "probeID")
    gene_probeID <- ann450kNew %>%
        dplyr::select(probeID, UCSC_RefGene_Name)
    full_predSex <- dmpSig %>%
        left_join(dplyr::select(gene_probeID, probeID, UCSC_RefGene_Name), by = "probeID")
    full_predSex <- full_predSex %>%
        filter(!UCSC_RefGene_Name == "")
    dim(full_predSex)
    print(dir)
    print(pre)
    write_csv(full_predSex, file=paste0(dir,"/",pre,"_dmp_MethFullResults.csv"))
    # Return the DMRs with full results including annotation, p-value, and delta beta
    return(full_predSex)
}
    

# dmpFinderAlt <- function(dat, pheno, cov, type = c("categorical", "continuous"),
#                       qCutoff = 1, shrinkVar = FALSE) {
#     # Check inputs
#     if (is(dat, "MethylSet")) {
#         .isMatrixBackedOrStop(dat, "dmpFinder")
#         M <- getM(dat)
#     } else if (is(dat, "DelayedMatrix")) {
#         stop("This function does not yet support DelayedArray objects")
#     } else {
#         stopifnot(is.numeric(dat))
#         M <- dat
#         if (is.vector(M)) M <- matrix(M, nrow = 1)
#     }
#     n <- length(pheno)
#     if (n != ncol(M)) stop("length of pheno does not equal number of samples")

#     if (type == "categorical") {
#         pheno <- factor(as.character(pheno))
#         design <- model.matrix(~0 + pheno + cov)
#         print(design)
#         fit <- lmFit(M, design)
#         if (shrinkVar) {
#             fit <- contrasts.fit(fit, contrasts(pheno))
#             fit <- eBayes(fit)
#             tab <- data.frame(
#                 intercept = fit$coefficients[, 1],
#                 f = fit[["F"]],
#                 pval = fit[["F.p.value"]])
#         } else {
#             fit1 <- lmFit(M)
#             RSS1 <- rowSums2((M - fitted(fit1))^2)
#             RSS <- rowSums2((M - fitted(fit))^2)
#             df1 <- length(levels(pheno)) - 1
#             df2 <- n - length(levels(pheno))
#             Fstat <- ((RSS1 - RSS)/df1)/(RSS/df2)
#             if (df2 > 1e+06) {
#                 F.p.value <- pchisq(df1 * Fstat, df1, lower.tail = FALSE)
#             }
#             else {
#                 F.p.value <- pf(Fstat, df1, df2, lower.tail = FALSE)
#             }
#             tab <- data.frame(
#                 intercept = fit$coefficients[, 1],
#                 f = Fstat,
#                 pval = F.p.value)
#         }
#     }
#     else if (type == "continuous") {
#         design <- model.matrix(~0 + pheno + cov)
#         fit <- lmFit(M, design)
#         if (shrinkVar) {
#             fit <- eBayes(fit)
#             sigma <- sqrt(fit$s2.post)
#             df <- fit$df.prior + fit$df.residual
#         } else {
#             sigma <- fit$sigma
#             df <- fit$df.residual
#         }
#         coef <- fit$coefficients
#         if (is.vector(coef)) coef <- matrix(coef, ncol = 2)
#         stdev <- fit$stdev.unscaled
#         if (is.vector(stdev)) stdev <- matrix(stdev, ncol = 2)
#         t <- coef[, 2] / (stdev[, 2] * sigma)
#         pval <- 2 * pt(abs(t), df = df, lower.tail = FALSE)
#         tab <- data.frame(
#             intercept = coef[, 1],
#             beta = coef[, 2],
#             t = t,
#             pval = pval)
#     }
#     p0 <- pi0.est(tab$pval[!is.na(tab$pval)])$p0
#     tab$qval <- qvalue.cal(tab$pval, p0)
#     if (qCutoff < 1) tab <- tab[tab$qval <= qCutoff,]
#     o <- order(tab$pval)
#     tab <- tab[o, ]
#     ## tab <- cbind(tab, annotate(rownames(tab)))
#     if (nrow(tab) == 0) {
#         message(sprintf("No significant DMPs at FDR cutoff of %s.", qCutoff))
#     }
#     tab
# }

GRSet <- readRDS(opt$grset)
pheno <- read_csv(opt$filename)

keep <- !is.na(pheno$age)
pheno_filt <- pheno[keep,]
GRSet_filt <- GRSet[,pheno_filt$geo_accession]

pheno_sex <- factor(pheno_filt$sex, levels=c("M","F"))
pheno_age <- pheno_filt$age
pheno_agecat <- factor(pheno_filt$agecat)

beta <- getBeta(GRSet_filt)
# dmp <- dmpFinder_custom(beta, pheno=pheno_sex, type="categorical")
# saveRDS(dmp, paste0(opt$outDir, "/", opt$prefix,"_DMPs_sexdiff.rds", sep=""))
dmpAgeAdj <- dmpFinder_custom(beta, pheno=pheno_sex, cov=pheno_age, pre=paste0(opt$prefix,"_ageadj"),dir=opt$outDir)
saveRDS(dmpAgeAdj, paste0(opt$outDir, "/", opt$prefix,"_DMPs_AgeAdj_sexdiff.rds", sep=""))
dmpAgeCatAdj <- dmpFinder_custom(beta, pheno=pheno_sex, cov=pheno_agecat, pre=paste0(opt$prefix,"_agecatadj"), dir=opt$outDir)
saveRDS(dmpAgeCatAdj, paste0(opt$outDir, "/", opt$prefix,"_DMPs_AgeCatAdj_sexdiff.rds", sep=""))


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
# findDMRs <- function(GR, design, cutoff, perm, dir, pre, pheno){
#   dmp <- dmpFinder(GR, pheno = sex, type = "categorical") 
  


#   dmrs <- bumphunter(GR, design = design, coef=2,
#                     cutoff = cutoff, B=perm, type="Beta", nullMethod=c("bootstrap"))
#   saveRDS(dmrs, paste0(dir, "/", pre,"_DMRs_sexdiff.rds", sep=""))

#   dmrs_table = dmrs$table
#   dmrs_table <- as.data.frame(dmrs_table)
#   write.table(dmrs_table, paste0(dir, "/", pre,"_sexdmrs_raw.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
#   dmrs_table$strand <- "."
#   dmrs_table <- unite(dmrs_table,name, chr, start, end, sep = "_", remove = F)

  # Write to bed file for processing with bedtools closest
# write.table(filter(dmrs_table, fwer < 0.05)[,c("chr","start","end","name","fwer","strand")], file = paste0(dir, "/", pre, "_sexdmrs_sigRegions.bed"), sep = "\t",quote = F, col.names = F, row.names = F)

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
#   return(dmrs)
# }

# GRSet <- readRDS(opt$grset)
# pheno <- read_csv(opt$filename)

# keep <- !is.na(pheno$age)
# pheno_filt <- pheno[keep,]
# GRSet_filt <- GRSet[,pheno_filt$geo_accession]

# pheno_sex <- pheno_filt$sex
# pheno_age <- pheno_filt$age
# pheno_agecat <- pheno_filt$agecat

# beta <- getBeta(GRSet_filt)


# registerDoParallel(cores = 4)

# # Run bumphunter to find DMRs by sex with age adjustment
# designMatrix1 <- model.matrix(~ pheno_sex + pheno_age)
# dmrs_all_age <- findDMRs(GRSet_filt, designMatrix1, cut, num, opt$outDir, paste0(opt$prefix,"_ageadj"), pheno_filt)

# # Do the same but with age as a category
# designMatrix2 <- model.matrix(~ pheno_sex + pheno_agecat)
# dmrs_all_agecat <- findDMRs(GRSet_filt, designMatrix2, cut, num, opt$outDir, paste0(opt$prefix,"_agecatadj"), pheno_filt)

# Break dataset into subgroups and analyze for DMRs within each group

subgroups <- c("SHH","WNT","Group3","Group4")

lsDMPsSubAge <- list()
lsDMPsSubAgeCat <- list()
for (sub in subgroups){
  pheno_sub <- filter(pheno_filt, subgroup==sub)
  GRSet_sub <- GRSet_filt[,pheno_sub$geo_accession]
  beta <- getBeta(GRSet_sub)
  prefix <- paste0(opt$prefix, "_",sub)
  pheno_sex <- factor(pheno_sub$sex, levels=c("M","F"))
  pheno_age <- pheno_sub$age
  pheno_agecat <- factor(pheno_sub$agecat)
  lsDMPsSubAge[[sub]] <- dmpFinder_custom(beta, pheno=pheno_sex, cov=pheno_age, pre=paste0(prefix,"_ageadj"), dir=opt$outDir)
  lsDMPsSubAgeCat[[sub]] <- dmpFinder_custom(beta, pheno=pheno_sex, cov=pheno_agecat, pre=paste0(prefix,"_ageacatadj"),dir=opt$outDir)
}