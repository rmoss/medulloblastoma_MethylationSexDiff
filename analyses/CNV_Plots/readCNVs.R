library(tidyverse)

options(scipen = 20)

setwd("/home/spectorl/moss0134/projects/research/MBSexDiff/results/tables")
CNVdir <- "./"
d <- list.files(path=CNVdir, pattern = "_CNVsegments.csv")
readCounts <- function(name){
  tmp <- read_csv(paste0(CNVdir,name, sep=""))
  #colnames(tmp)[7] <- "counts"
  #colnames(tmp)[1] <- "locID"
  #tmp <- filter(tmp, Chr %in% c(1:22, "X","Y"))
  cohort <- str_remove(name, "_CNVsegments.csv")
  cohort <- gsub("_.*","",cohort)
  tmp$cohort <- cohort
  #tmp <- tmp[,c("locID","counts","sampleID")]
  tmp
}
all <- map_dfr(d, readCounts)
all_bed <- all[,c("chrom","loc.start","loc.end","ID","seg.mean","cohort","subgroup","pval","seg.median")]
write.table(all_bed, file = "CNVcombined.bed", col.names = T, row.names = F, quote = F, sep = "\t")

# for (sub in unique(all_bed$subgroup)) {
#   print(sub)
#   cnv <- filter(all_bed, subgroup==sub)
#   write.table(cnv, file =paste0(sub,"_CNV.bed"), col.names=T, row.names=F, quote=F, sep = "\t")
# }

for (co in unique(all_bed$cohort)) {
  for (sub in unique(all_bed$subgroup)) {
    print(co)
    print(sub)
    cnv <- filter(all_bed, cohort==co & subgroup==sub)
    write.table(cnv, file =paste0(co,"_",sub,"_CNV.bed"), col.names=T, row.names=F, quote=F, sep = "\t")
  }
}




#all_wide <- pivot_wider(all, names_from=sampleID, values_from=counts)
#write.table(all_wide, file = "ES_R21_allmatrix.counts", col.names = T, row.names = F, quote = F, sep = "\t")
