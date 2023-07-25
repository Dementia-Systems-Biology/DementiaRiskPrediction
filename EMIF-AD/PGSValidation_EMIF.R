library(tidyverse)
library(caret)
library(pROC)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_fil.RData")
load("~/EMIF/PGS_EMIF.RData")

rownames(metaData_fil) <- metaData_fil$Sample_Name
rownames(PGS_all) <- PGS_all$ID

samples <- intersect(metaData_fil$Sample_Name, PGS_all$ID)

dat_fil <- metaData_fil[samples,]
PGS_fil <- PGS_all[samples,]


pROC::auc(pROC::roc(dat_fil$Hypertension[!is.na(dat_fil$Hypertension)],
          PGS_fil$SBPauto[!is.na(dat_fil$Hypertension)]))

pROC::auc(pROC::roc(dat_fil$Obesity[!is.na(dat_fil$Obesity)],
                    PGS_fil$BMI[!is.na(dat_fil$Obesity)]))

pROC::auc(pROC::roc(ifelse(dat_fil$Eduy <= 9,1,0)[!is.na(dat_fil$Eduy)],
                    PGS_fil$EA22[!is.na(dat_fil$Eduy)]))

pROC::auc(pROC::roc(dat_fil$Depression[!is.na(dat_fil$Depression)],
                    PGS_fil$MDD[!is.na(dat_fil$Depression)]))

pROC::auc(pROC::roc(dat_fil$CardiovascularDis[!is.na(dat_fil$CardiovascularDis)],
                    PGS_fil$CAD[!is.na(dat_fil$CardiovascularDis)]))

pROC::auc(pROC::roc(ifelse(dat_fil$Alcohol >=1,1,0)[!is.na(dat_fil$Alcohol)],
                    PGS_fil$CAD[!is.na(dat_fil$Alcohol)]))
          