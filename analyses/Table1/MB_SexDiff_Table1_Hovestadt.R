library(readr)
library(tableone)
library(dplyr)
library(table1)
library(boot)
library(xlsx)

#### Read in Hovestadt data ####
phenoHov <- read_csv("./data/pheno_Hovestadt_subgrouped.csv")

# create categories for age
phenoHov <- phenoHov %>%
  mutate(
    agecat = case_when(
      `age (yrs):ch1` < 5 ~ 1,
      `age (yrs):ch1` <10 ~ 2,
      `age (yrs):ch1` <15 ~ 3,
      `age (yrs):ch1` <20 ~ 4,
      `age (yrs):ch1` >=20 ~ 5,
      TRUE ~ NA_real_
    )
  )

# Factor the basic variables that we're interested in
phenoHov$`gender:ch1` <- 
  factor(phenoHov$`gender:ch1`, 
         levels=c("F","M"),
         labels=c("Females", # Reference
                  "Males"))
phenoHov$`molecular subgroup:ch1` <-
  factor(phenoHov$`molecular subgroup:ch1`,
         levels=c("SHH","WNT","Group3","Group4"),
         labels=c("SHH","WNT","Group 3","Group 4"))

phenoHov$agecat <-  
  factor(phenoHov$agecat,
         levels=c(1:5),
         labels=c("0-<5","5-<10","10-<15","15-<20",">=20"))

#### Create table ####
table1::table1( ~ agecat | `molecular subgroup:ch1`*`gender:ch1`, data=phenoHov, droplevels=F)

# Create non-sex stratified table to get numbers for top row
table1::table1( ~ `gender:ch1` | `molecular subgroup:ch1`, data=phenoHov, droplevels=F)

#### Calculate stats ####
chisq.test(phenoHov$`molecular subgroup:ch1`, phenoHov$`gender:ch1`)
chisq.test(table(phenoHov$`molecular subgroup:ch1`, phenoHov$`gender:ch1`))


Group3 <- phenoHov[which(phenoHov$`molecular subgroup:ch1` == "Group 3"),]
Group4 <- phenoHov[which(phenoHov$`molecular subgroup:ch1`== "Group 4"),]
SHH <- phenoHov[which(phenoHov$`molecular subgroup:ch1` == "SHH"),]
WNT <- phenoHov[which(phenoHov$`molecular subgroup:ch1` == "WNT"),]

fisher.test(SHH$agecat, SHH$`gender:ch1`)
fisher.test(WNT$agecat, WNT$`gender:ch1`)
fisher.test(Group3$agecat, Group3$`gender:ch1`)
fisher.test(Group4$agecat, Group4$`gender:ch1`)
fisher.test(phenoHov$agecat, phenoHov$`gender:ch1`)

fisher.test(SHH$pathology.ch1, SHH$predictedSex)
fisher.test(WNT$pathology.ch1, WNT$predictedSex)
fisher.test(Group3$pathology.ch1, Group3$predictedSex)
fisher.test(Group4$pathology.ch1, Group4$predictedSex)
fisher.test(pheno_filt$pathology.ch1,pheno_filt$predictedSex)

#write_csv(phenoHov, file="./data/Heidelberg/Heidelberg_pheno.csv")

#### Lichter data integration ####
# Read in Lichter supplemental tables and ses which match to Heidelberg cohort: https://www.nature.com/articles/nature22973#Sec44
phenoLichter <- read.xlsx("./data/Heidelberg/41586_2017_BFnature22973_MOESM1_ESM.xlsx", 1)
pheno <- left_join(phenoHov,phenoLichter,by=c("title"="PID"))

LichterS1 <- read.xlsx("./data/Heidelberg/Lichter_TableS1_final.xls",1)
LichterS1 <- LichterS1[,1:14]
pheno2012 <- 

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






