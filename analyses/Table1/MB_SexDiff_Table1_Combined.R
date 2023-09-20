library(readr)
library(tableone)
library(dplyr)
library(table1)
library(boot)

#### Data directory
dir <- "/Volumes/GoogleDrive/My\ Drive/Research/SexDiffMedullo/Analysis_RMM_new/data/"

#### Read in data ####
phenoHov <- read_csv(paste0(dir,"pheno_Hovestadt_subgrouped.csv"))
phenoCav <- read_csv(paste0(dir,"pheno_CavalliCombined_subgrouped.csv"))
phenoNew <- read_csv(paste0(dir,"Newcastle_pheno_filt_MethPedsubgroup.csv"))

# Standardize columns, column names, and add column for source
phenoHovF <- phenoHov[,c(1:2,4:8)]
colnames(phenoHovF) <- c("sample","geo_accession","description","age","sex","subgroup","type")
phenoHovF$source <- c("Hovstadt")
# phenoCavF <- phenoCav[,c(2,15,24,3,5,12,21,8:10)]
# colnames(phenoCavF) <- c("sample","geo_accession","description","age","sex","subgroup","type","Histology","Metastatic Disease","Vital Status")
phenoCavF <- phenoCav[,c(2,15,24,3,5,12,21)]
colnames(phenoCavF) <- c("sample","geo_accession","description","age","sex","subgroup","type")
phenoCavF$source <- c("Cavalli")
phenoCavF$sex <- 
  factor(phenoCavF$sex, 
         levels=c("Female","Male"),
         labels=c("F", # Reference
                  "M"))
phenoNewF <- phenoNew[,c(10,1,2,5,7,14,8)]
#phenoNewF$type <- c("NA")
colnames(phenoNewF) <- c("sample","geo_accession","description","age","sex","subgroup","type")
phenoNewF$source <- c("Newcastle")

# combine the 3 datasets into one
pheno <- rbind(phenoHovF,phenoCavF,phenoNewF)

# create categories for age
pheno <- pheno %>%
  mutate(
    agecat = case_when(
      age < 5 ~ 1,
      age <10 ~ 2,
      age <15 ~ 3,
      age <20 ~ 4,
      age >=20 ~ 5,
      TRUE ~ NA_real_
    )
  )

# standardize subgroup names
pheno$subgroup <- sub("MB_","",pheno$subgroup)
pheno$subgroup <- sub("Gr3","Group3",pheno$subgroup)
pheno$subgroup <- sub("Gr4","Group4",pheno$subgroup)

# filter out non-MB samples
MBsubs <- c("SHH","WNT","Group3","Group4")
phenoF <- filter(pheno, subgroup %in% MBsubs)

# Factor the basic variables that we're interested in
phenoF$sex <- 
  factor(phenoF$sex, 
         levels=c("F","M"),
         labels=c("Females", # Reference
                  "Males"))
phenoF$subgroup <-
  factor(phenoF$subgroup,
         levels=c("SHH","WNT","Group3","Group4"),
         labels=c("SHH","WNT","Group 3","Group 4"))

phenoF$agecat <-  
  factor(phenoF$agecat,
         levels=c(1:5),
         labels=c("0-<5","5-<10","10-<15","15-<20",">=20"))

#### Create table ####
table1::table1( ~ agecat + source | subgroup*sex, data=phenoF, droplevels=F)

# Create non-sex stratified table to get numbers for top row
table1::table1( ~ sex | subgroup, data=phenoF, droplevels=F)

#### Calculate stats ####
chisq.test(phenoF$subgroup, phenoF$sex)
chisq.test(table(phenoF$subgroup, phenoF$sex))


Group3 <- phenoF[which(phenoF$subgroup == "Group 3"),]
Group4 <- phenoF[which(phenoF$subgroup== "Group 4"),]
SHH <- phenoF[which(phenoF$subgroup == "SHH"),]
WNT <- phenoF[which(phenoF$subgroup == "WNT"),]

fisher.test(SHH$agecat, SHH$sex)
fisher.test(WNT$agecat, WNT$sex)
fisher.test(Group3$agecat, Group3$sex)
fisher.test(Group4$agecat, Group4$sex)
fisher.test(phenoF$agecat, phenoF$sex)

fisher.test(SHH$source,SHH$sex)
fisher.test(WNT$source,WNT$sex)
fisher.test(Group3$source,Group3$sex)
fisher.test(Group4$source,Group4$sex)
fisher.test(phenoF$source,phenoF$sex)



#### Experimentation ####

tab1 <- CreateTableOne(data = pheno_filt,
                       vars = c("predictedSex", "agecat","pathology.ch1"),
                       factorVars = c("agecat","pathology.ch1"),
                       strata = c("predictedSex"), # Which variable should be used for stratification
                       addOverall = T, includeNA=TRUE # This allows us to include an Overall column in the summary table
)

tab1 <- CreateTableOne(data = pheno_filt,
               vars = c("predictedSex", "agecat","pathology.ch1"),
               factorVars = c("agecat","pathology.ch1"),
               strata = c("predictedSex","subgroup"), # Which variable should be used for stratification
               addOverall = T, includeNA=TRUE # This allows us to include an Overall column in the summary table
)
print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
              
tab2 <- CreateTableOne(data = pheno,
                       vars = c("predictedSex"),#, "agecat","histology","Met.status..1.Met..0.M0", "Dead", "", "dibpat", "height"),
                       #factorVars = c("factor(predictedSex)","histology","Met.status..1.Met..0.M0","Dead"),
                       strata = c("subgroup","predictedSex"), # Which variable should be used for stratification
                       addOverall = T, includeNA=TRUE # This allows us to include an Overall column in the summary table
)
print(tab2, showAllLevels = TRUE, formatOptions = list(big.mark = ","))


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- lalonde[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ lalonde$treat)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(lalonde$treat)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

table1(~ predictedSex + factor(agecat[!is.na(agecat)]) + factor(histology) + factor(Met.status..1.Met..0.M0.) + factor(Dead)  | Subgroup*predictedSex, data=pheno, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)

table1(~ geo_accession) | subgroup*predictedSex, data=pheno_filt)
table1(~ factor(agecat) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$agecat)),])
table1(~ factor(histology) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$histology)),])
table1(~ factor(Met.status..1.Met..0.M0.) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$Met.status..1.Met..0.M0.)),])
table1(~ factor(Dead) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$Dead)),])






