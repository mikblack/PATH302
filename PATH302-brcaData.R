## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=F-------------------------------------------------------------
## dir.create('Rlibs')
## newPath = c(.libPaths(), paste0(getwd(), "/Rlibs"))
## .libPaths(newPath)
## 
## install.packages(c('magrittr','gplots','httr'), lib='Rlibs')
## 
## install.packages('BiocInstaller', repos='http://www.bioconductor.org/packages/3.6/bioc', lib='Rlibs')
## library(BiocInstaller, lib.loc="Rlibs")
## biocLite("Biobase", lib="Rlibs")
## library("Biobase", lib.loc="Rlibs")

## ---- warning=FALSE, message=FALSE---------------------------------------
library(ggplot2)
library(magrittr)
library(gplots)
library(survival)
library(httr)

## ---- echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE----------------
## load(url("https://github.com/mikblack/PATH302/raw/master/uppsalaCohort.RData"))

## ---- echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE----------------
library(Biobase)
load("uppsalaCohort.RData")

## ------------------------------------------------------------------------
uppExp = exprs(upp)

## ---- eval=FALSE---------------------------------------------------------
## View(uppExp)

## ------------------------------------------------------------------------
uppAnnot = featureData(upp)@data

## ---- eval=FALSE---------------------------------------------------------
## View(uppAnnot)

## ------------------------------------------------------------------------
uppClin  = phenoData(upp)@data

## ---- eval=FALSE---------------------------------------------------------
## View(uppClin)

## ------------------------------------------------------------------------
uppClinSmall = uppClin[,c("size", "age", "er", "grade", "pgr", 
                          "node", "t.rfs", "e.rfs","treatment")]

## ---- eval=FALSE---------------------------------------------------------
## View(uppClinSmall)

## ------------------------------------------------------------------------
uppClinSmall$t.rfs = uppClinSmall$t.rfs / 365

## ------------------------------------------------------------------------
attach(uppClinSmall)

## ------------------------------------------------------------------------
grade

## ---- echo=FALSE---------------------------------------------------------
desc = c("Tumour Size (cm)",
         "Patient age at diagnosis (years)",
         "Estrogen receptor status (0=negative, 1=positive)",
         "Tumour Grade (1, 2, or 3)",
         "Progesterone receptor status (0=negative, 1=positive)",
         "Lymph node status (0=negative, 1=positive)",
         "Recurrence free survival time (years)",
         "Recurrence free survival event (0=censored, 1=recurrence)",
         "Treatment (0=none, 2=endocrine treatment or chemotherapy")
datTab = data.frame(Variable = names(uppClinSmall),
                    Description = desc)
knitr::kable(datTab)

## ---- warning=FALSE, message=FALSE---------------------------------------
ggplot(data=uppClinSmall, aes(x=size)) + 
  geom_histogram(fill="white", colour="black") +
  xlab("Tumour Size") +
  ylab("Frequency") + 
  ggtitle("Uppsala cohort: tumour size")

## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
ggplot(data=uppClinSmall, aes(x=age)) + 
  geom_histogram(fill="white", colour="black") +
  xlab("Age") +
  ylab("Frequency") + 
  ggtitle("Uppsala cohort: patient age")

## ------------------------------------------------------------------------
ggplot(data=uppClinSmall, aes(x=as.factor(grade),y=size)) + 
  geom_boxplot()

## ---- warning=FALSE, message=FALSE---------------------------------------
ggplot(data=subset(uppClinSmall, !is.na(grade)), aes(x=as.factor(grade),y=size)) + 
  geom_boxplot() +
  ylim(0,7) +
  ggtitle("Uppsala cohort: tumour size versus grade")

## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
ggplot(data=subset(uppClinSmall, !is.na(er)), aes(x=as.factor(er),y=size)) + 
  geom_boxplot() + 
  ylim(0,7) +
  ggtitle("Uppsala cohort: tumour size versus ER status")

## ------------------------------------------------------------------------
table(er)

## ------------------------------------------------------------------------
table(grade)
prop.table( table(grade) )

## ------------------------------------------------------------------------
table(er, grade)

## ------------------------------------------------------------------------
## Overall proportions:
prop.table(table(er, grade))

## Proportions by row (rows sum to 1)
prop.table(table(er, grade),1)

## Proportions by columns (columns sum to 1)
prop.table(table(er, grade),2)

## ------------------------------------------------------------------------
fisher.test(table(er, grade))

## ------------------------------------------------------------------------
chisq.test(table(er, grade))

## ---- eval=FALSE---------------------------------------------------------
## ## DON'T RUN THIS CODE:
## library(genefu)
## PAM50Preds = molecular.subtyping(sbt.model = "pam50", data=t(uppExp),
##                                 annot=uppAnnot, do.mapping=TRUE)

## ------------------------------------------------------------------------
## DO RUN THIS CODE:
table(PAM50Preds$subtype)

## ---- warning=FALSE, message=FALSE---------------------------------------
uppClinSmall$subtype = PAM50Preds$subtype
attach(uppClinSmall)

## ------------------------------------------------------------------------
table(subtype, er)

table(subtype, er) %>% 
  prop.table(., 1) %>% 
  round(., 3)

## ------------------------------------------------------------------------
round( prop.table( table(subtype, er), 1) ,3)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
## table(subtype, grade)
## table(subtype, grade) %>%
##   prop.table(., 1) %>%
##   round(., 3)

## ------------------------------------------------------------------------
plot( survfit(Surv(t.rfs, e.rfs) ~1 ), 
      xlab="Time (years)", 
      ylab = "Proportion recurrence free")

## ------------------------------------------------------------------------
plot( survfit(Surv(t.rfs, e.rfs) ~ er ), col=1:2, 
      xlab="Time (years)", 
      ylab = "Proportion recurrence free")
legend('bottomleft', c("ER-", "ER+"), fill=1:2)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
## plot( survfit(Surv(t.rfs, e.rfs) ~ grade ), col=1:3,
##       xlab="Time (years)",
##       ylab = "Proportion recurrence free")
## legend('bottomleft', c("G1", "G2", "G3"), fill=1:3)

## ----]-------------------------------------------------------------------
plot( survfit(Surv(t.rfs, e.rfs) ~ subtype ), col=1:5, 
      xlab="Time (years)", 
      ylab = "Proportion recurrence free")
groups <- names(table(subtype))
legend('bottomleft', groups, fill=1:5)

## ------------------------------------------------------------------------
survdiff(Surv(t.rfs, e.rfs) ~ subtype)

## ---- echo=FALSE, eval=FALSE---------------------------------------------
## survdiff(Surv(t.rfs, e.rfs) ~ grade)

## ------------------------------------------------------------------------
esr1Probes = uppAnnot$probe[ na.omit(uppAnnot$Gene.symbol == 'ESR1') ]
esr1Probes

## ---- eval=FALSE, echo=FALSE---------------------------------------------
## par(mfrow=c(2,5))
## for(i in esr1Probes){
##   x = uppExp[match(i, rownames(uppExp)),]
##   boxplot(x ~ er)
## }

## ------------------------------------------------------------------------
match("205221_at", rownames(uppExp))

## ------------------------------------------------------------------------
esr1probe = match("205221_at", rownames(uppExp))

## ------------------------------------------------------------------------
esr1Dat = uppExp[esr1probe, ]

## ---- evel=FALSE, echo=FALSE---------------------------------------------
## Switched to the above to make the code more readable
esr1Dat = uppExp[match("205221_at", rownames(uppExp)), ]

## ---- warning=FALSE, message=FALSE---------------------------------------
uppClinSmall$esr1Dat = esr1Dat
attach(uppClinSmall)

## ------------------------------------------------------------------------
ggplot(data=subset(uppClinSmall, !is.na(er)), aes(x=as.factor(er), y=esr1Dat)) + geom_boxplot()

## ---- echo=FALSE---------------------------------------------------------
#prolifGenes = c("OSTC","MCM6","RPA3","MCM7","PCNA","XRCC6","KPNA2","ANLN","RNASEH2A",
#                "PBK","GMNN","RRM1","CDC45","MAD2L1","RAN","DUT","RRM2","CDK7",
#                "MLH3","SMC4","SMC3","POLD2","POLE2","BCCIP","GINS2","TREX1",
#                "BUB3","FEN1","DBF4B","MOB4","CCNE1","RPA1","POLE3","RFC4","MCM3",
#                "CHEK1","CCND1","CDC37")
prolifGenes = c("MAD2L1", "RRM2", "ANLN", "MCM6", "PBK", "GINS2", "KPNA2", "PCNA")
knitr::kable(data.frame(Gene=prolifGenes, 
                        Probe=uppAnnot$probe[match(prolifGenes, uppAnnot$Gene.symbol)]))

## ------------------------------------------------------------------------
prolifGenes = c("MAD2L1", "RRM2", "ANLN", "MCM6", "PBK", "GINS2", "KPNA2", "PCNA")

## ------------------------------------------------------------------------
prolifRows = na.omit(match(prolifGenes, uppAnnot$Gene.symbol))

## ------------------------------------------------------------------------
prolifRows

## ------------------------------------------------------------------------
prolifDat = uppExp[prolifRows, ]

## ------------------------------------------------------------------------
prolifDatScale = t(scale(t(prolifDat)))
prolifDatScale[prolifDatScale < -3] = -3
prolifDatScale[prolifDatScale >  3] =  3

## ------------------------------------------------------------------------
heatmap.2(prolifDatScale, trace='none', scale='none', col='bluered', 
          labRow = uppAnnot$Gene.symbol[prolifRows])

## ------------------------------------------------------------------------
prolifMean = colMeans(prolifDat)

## ------------------------------------------------------------------------
ord = order(prolifMean)
prolifCol = bluered(length(prolifMean))[rank(prolifMean)]

## ------------------------------------------------------------------------
heatmap.2(prolifDatScale[,ord], trace='none', scale='none', col='bluered', 
          labRow = uppAnnot$Gene.symbol[prolifRows], Colv=FALSE,
          ColSideColors=prolifCol[ord])

## ---- warning=FALSE, message=FALSE---------------------------------------
uppClinSmall$prolifMean = prolifMean
attach(uppClinSmall)

## ------------------------------------------------------------------------
ggplot(data=subset(uppClinSmall, !is.na(grade)), aes(x=as.factor(grade), y=prolifMean)) + geom_boxplot()

## ---- echo=FALSE, eval=FALSE---------------------------------------------
## ggplot(data=uppClinSmall, aes(x=as.factor(subtype), y=prolifMean)) + geom_boxplot()

## ------------------------------------------------------------------------
prolifHilo = ifelse(prolifMean > median(prolifMean), "HighProilif", "LowProlif")

## ------------------------------------------------------------------------
table(prolifHilo)

## ------------------------------------------------------------------------
table(grade, prolifHilo)

## ------------------------------------------------------------------------
plot( survfit(Surv(t.rfs, e.rfs) ~ prolifHilo ), col=1:2, 
      xlab="Time (years)", ylab = "Proportion recurrence free")
groups <- names(table(prolifHilo))
legend('bottomleft', groups, fill=1:2)

## ---- eval=FALSE---------------------------------------------------------
## load(url("https://github.com/mikblack/PATH302/raw/master/uppsalaObjects.RData"))

