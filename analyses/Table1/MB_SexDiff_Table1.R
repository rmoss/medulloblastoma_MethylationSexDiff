library(readr)
library(tableone)
library(dplyr)
library(table1)
library(boot)

pheno <- read_csv("./data/pheno_combined.csv")

pheno <- pheno %>%
  select(1:51) %>%
  mutate(
    agecat = case_when(
      Age < 5 ~ 1,
      Age <10 ~ 2,
      Age <15 ~ 3,
      Age <20 ~ 4,
      Age >=20 ~ 5,
      TRUE ~ NA_real_
    )
  )

tab1 <- CreateTableOne(data = pheno,
               vars = c("predictedSex", "agecat","histology","Met.status..1.Met..0.M0", "Dead", "", "dibpat", "height"),
               factorVars = c("factor(predictedSex)","histology","Met.status..1.Met..0.M0","Dead"),
               strata = c("Subgroup","Gender"), # Which variable should be used for stratification
               addOverall = T, includeNA=TRUE # This allows us to include an Overall column in the summary table
)
print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
              
tab2 <- CreateTableOne(data = pheno,
                       vars = c("predictedSex"),#, "agecat","histology","Met.status..1.Met..0.M0", "Dead", "", "dibpat", "height"),
                       #factorVars = c("factor(predictedSex)","histology","Met.status..1.Met..0.M0","Dead"),
                       strata = c("Subgroup","Gender"), # Which variable should be used for stratification
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

table1(~ factor(agecat) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$agecat)),])
table1(~ factor(histology) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$histology)),])
table1(~ factor(Met.status..1.Met..0.M0.) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$Met.status..1.Met..0.M0.)),])
table1(~ factor(Dead) | Subgroup*predictedSex, data=pheno[!(is.na(pheno$Dead)),])



chisq.test(pheno$Subgroup, pheno$predictedSex)
chisq.test(table(pheno$Subgroup,pheno$predictedSex))


Group3 <- pheno[which(pheno$Subgroup == "Group3"),]
Group4 <- pheno[which(pheno$Subgroup == "Group4"),]
SHH <- pheno[which(pheno$Subgroup == "SHH"),]
WNT <- pheno[which(pheno$Subgroup == "WNT"),]

fisher.test(SHH$agecat, SHH$predictedSex)
fisher.test(WNT$agecat, WNT$predictedSex)
fisher.test(Group3$agecat, Group3$predictedSex)
fisher.test(Group4$agecat, Group4$predictedSex)

fisher.test(SHH$histology, SHH$predictedSex)
fisher.test(WNT$histology, WNT$predictedSex)
fisher.test(Group3$histology, Group3$predictedSex)
fisher.test(Group4$histology, Group4$predictedSex)

fisher.test(SHH$Met.status..1.Met..0.M0., SHH$predictedSex)
fisher.test(WNT$Met.status..1.Met..0.M0., WNT$predictedSex)
fisher.test(Group3$Met.status..1.Met..0.M0., Group3$predictedSex)
fisher.test(Group4$Met.status..1.Met..0.M0., Group4$predictedSex)

fisher.test(SHH$Dead, SHH$predictedSex)
fisher.test(WNT$Dead, WNT$predictedSex)
fisher.test(Group3$Dead, Group3$predictedSex)
fisher.test(Group4$Dead, Group4$predictedSex)

