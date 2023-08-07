# ============================================================================ #
# File: GetStatistics_EMIF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Get statstics of model performance in independent test set of 
#              the EMIF-AD cohort.
# ============================================================================ #

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)
library(mltools)
library(reportROC)


################################################################################

# Epi-MCI Model

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# MCI and NL only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

# Set factors
Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

table(Y_test$Y)
table(Y_train$Y)

# Load models
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData")
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_EN.RData")
load("Models/EMIF_Models/MRS/Data/Fit_EMIF_MCI_sPLS.RData")

# Get case/control threshold (maximum geometric mean of sensitivity and specificity)
predictions <- fit$pred
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

# Make confusion matrix for cross-validation performance
confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value in independent test set
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y

# Confusion matrix
confusionMatrix(pred,obs, positive = "MCI")

# MCC
mcc(pred, obs)

# AUROC
pROC::auc(obs,RF$MCI)


################################################################################

# Epi-AD Model

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# MCI and NL only
X_train <- X_train[(Y_train$Diagnosis == "AD") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "AD")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "AD") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "AD") | (Y_test$Diagnosis == "NL"),]

# Set factors
Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","AD"),
                    levels = c("Control", "AD"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","AD"),
                   levels = c("Control", "AD"))

table(Y_test$Y)
table(Y_train$Y)


# Load models
load("Models/EMIF_Models/MRS/Fit_EMIF_AD_RF.RData")
load("Models/EMIF_Models/MRS/Fit_EMIF_AD_EN.RData")
load("Models/EMIF_Models/MRS/Fit_EMIF_AD_sPLS.RData")

# Get case/control threshold (maximum geometric mean of sensitivity and specificity)
predictions <- fit$pred
test <- pROC::roc(predictions$obs,predictions$AD)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$AD > threshold, "AD", "Control"), levels = c("Control", "AD"))

# Make confusion matrix for cross-validation performance
confusionMatrix(pred,predictions$obs, positive = "AD")

# Get observed and predicted value in independent test set
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$AD > threshold, "AD", "Control"), levels = c("Control", "AD"))
obs <- Y_test$Y

# Confusion matrix
confusionMatrix(pred,obs, positive = "AD")

# MCC
mcc(pred, obs)

# AUROC
pROC::auc(obs,RF$AD)


################################################################################

# MRS/CSF Model

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# Add CSF biomarkers
load("EMIF-AD/Data/metaData_fil.RData")
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]
# MCI and NL only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

# Get samples with MRSs and CSF biomarkers available
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]
X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])


# Load models:

# CSF only
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbioonly.RData")
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbioonly.RData")
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbioonly.RData")

# CSF + MRSs
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbio.RData")
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbio.RData")
load("Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbio.RData")

# Get case/control threshold (maximum geometric mean of sensitivity and specificity)
predictions <- fit$pred
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

# Make confusion matrix for cross-validation performance
confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value in independent test set
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y

# Confusion matrix
confusionMatrix(pred,obs, positive = "MCI")

# MCC
mcc(pred, obs)

# AUROC
pROC::auc(obs,RF$MCI)

################################################################################

# MRS/PGS Model

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# Add PGSs:

# Load data
load("EMIF-AD/Data/metaData_fil.RData")
load("EMIF-AD/Data/PGS_EMIF_AD.RData")

# Get PGSs
output <- inner_join(PGS_all,metaData_fil[,c("X","Sample_Name")],by = c("ID" = "Sample_Name"))
rownames(output) <- output$X
samples <- output$X
output <- output[,2:13]
colnames(output) <- paste0(colnames(output), "_PGS")

# MCI and control (NL) samples only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

# Get samples with MRSs and PGSs available
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]
X_train <- cbind.data.frame(X_train[rownames(Y_train),], output[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], output[rownames(Y_test),])

# Load models

# PGS/MRS model
load("EMIF/PGS_AD/Fit_EMIF_MCI_RF_PGS2.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_EN_PGS2.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_sPLS_PGS2.RData")

# PGS only
load("EMIF/PGS_AD/Fit_EMIF_MCI_RF_PGSonly.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_EN_PGSonly.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_sPLS_PGSonly.RData")

# Get case/control threshold (maximum geometric mean of sensitivity and specificity)
predictions <- fit$pred
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

# Make confusion matrix for cross-validation performance
confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value in independent test set
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y

# Confusion matrix
confusionMatrix(pred,obs, positive = "MCI")

# MCC
mcc(pred, obs)

# AUROC
pROC::auc(obs,RF$MCI)


################################################################################

# MRS/PGS/CSF Model

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# Add PGSs:

# Load data
load("EMIF-AD/Data/metaData_fil.RData")
load("EMIF-AD/Data/PGS_EMIF_AD.RData")

# Get PGSs
output <- inner_join(PGS_all,metaData_fil[,c("X","Sample_Name")],by = c("ID" = "Sample_Name"))
rownames(output) <- output$X
samples_PGS <- output$X
output <- output[,2:13]
colnames(output) <- paste0(colnames(output), "_PGS")

# Add CSF biomarkers
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore")
samples_csf <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                                  (!is.na(CSFbio$AB_Zscore)) &
                                  (!is.na(CSFbio$Ttau_ASSAY_Zscore))]
samples <- intersect(samples_PGS, samples_csf)

# Select MCI and Control (NL) samples only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

# Get samples with MRSs, PGSs, and CSF biomarkers available
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]
X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])
X_train <- cbind.data.frame(X_train[rownames(Y_train),],output[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], output[rownames(Y_test),])

# Load models:

# MRS/PGS/CSF models
load("EMIF/PGS_AD/Fit_EMIF_MCI_RF_CSF_PGS.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_EN_CSF_PGS2.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_sPLS_CSF_PGS.RData")

# PGS/CSF only
load("EMIF/PGS_AD/Fit_EMIF_MCI_RF_CSF_PGSonly.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_EN_CSF_PGSonly.RData")
load("EMIF/PGS_AD/Fit_EMIF_MCI_sPLS_CSF_PGSonly.RData")

# Get case/control threshold (maximum geometric mean of sensitivity and specificity)
predictions <- fit$pred
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

# Make confusion matrix for cross-validation performance
confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value in independent test set
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y

# Confusion matrix
confusionMatrix(pred,obs, positive = "MCI")

# MCC
mcc(pred, obs)

# AUROC
pROC::auc(obs,RF$MCI)

