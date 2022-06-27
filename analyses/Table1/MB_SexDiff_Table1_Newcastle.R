library(readr)
library(tableone)
library(dplyr)
library(table1)
library(boot)

#### Read in Newcastle data ####
pheno <- read_csv("./data/Newcastle_pheno_filt_MethPedsubgroup.csv")

# create categories for age
pheno <- pheno %>%
  mutate(
    agecat = case_when(
      age..yrs..ch1 < 5 ~ 1,
      age..yrs..ch1 <10 ~ 2,
      age..yrs..ch1 <15 ~ 3,
      age..yrs..ch1 <20 ~ 4,
      age..yrs..ch1 >=20 ~ 5,
      TRUE ~ NA_real_
    )
  )

# standardize subgroup names
pheno$subgroup <- sub("MB_","",pheno$subgroup)
pheno$subgroup <- sub("Gr","Group",pheno$subgroup)
# filter out non-MB samples
MBsubs <- c("SHH","WNT","Group3","Group4")
pheno_filt <- filter(pheno, subgroup %in% MBsubs)

# Factor the basic variables that we're interested in
pheno_filt$predictedSex <- 
  factor(pheno_filt$predictedSex, 
         levels=c("F","M"),
         labels=c("Females", # Reference
                  "Males"))
pheno_filt$subgroup <-
  factor(pheno_filt$subgroup,
         levels=MBsubs,
         labels=c("SHH","WNT","Group 3","Group 4"))

pheno_filt$agecat <-  
  factor(pheno_filt$agecat,
         levels=c(1:5),
         labels=c("0-<5","5-<10","10-<15","15-<20",">=20"))


write_csv(pheno_filt, file="./data/Newcastle/Newcastle_pheno.csv")

#### Create table ####
table1::table1( ~ agecat + factor(pathology.ch1) | subgroup*predictedSex, data=pheno_filt, droplevels=F)

# Create non-sex stratified table to get numbers for top row
table1::table1( ~ predictedSex | subgroup, data=pheno_filt, droplevels=F)

#### Calculate stats ####
chisq.test(pheno_filt$subgroup, pheno_filt$predictedSex)
chisq.test(table(pheno_filt$subgroup,pheno_filt$predictedSex))


Group3 <- pheno_filt[which(pheno_filt$subgroup == "Group 3"),]
Group4 <- pheno_filt[which(pheno_filt$subgroup == "Group 4"),]
SHH <- pheno_filt[which(pheno_filt$subgroup == "SHH"),]
WNT <- pheno_filt[which(pheno_filt$subgroup == "WNT"),]

fisher.test(SHH$agecat, SHH$predictedSex)
fisher.test(WNT$agecat, WNT$predictedSex)
fisher.test(Group3$agecat, Group3$predictedSex)
fisher.test(Group4$agecat, Group4$predictedSex)
fisher.test(pheno_filt$agecat,pheno_filt$predictedSex)

fisher.test(SHH$pathology.ch1, SHH$predictedSex)
fisher.test(WNT$pathology.ch1, WNT$predictedSex)
fisher.test(Group3$pathology.ch1, Group3$predictedSex)
fisher.test(Group4$pathology.ch1, Group4$predictedSex)
fisher.test(pheno_filt$pathology.ch1,pheno_filt$predictedSex)


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






