# ============================================================================ #
# File: ML_EMIF_CSF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Machine learning (ML) using MRSs and CSF biomarkers as 
#              variables in the EMIF-AD cohort.
# ============================================================================ #

###############################################################################

# Make predictions

###############################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)
library(pROC)
library(e1071)
library(ranger)
library(dplyr)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load training and test data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")

# Add CSF biomarkers
load("EMIF/metaData_fil.RData")
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]
# MCI and control (NL) only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))


# Combine MRSs and CSF biomarkers
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]
X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])

table(Y_test$Y)
table(Y_train$Y)

# Create index for repeated cross-validation
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           index = CVindex)


#*****************************************************************************#
# ElasticNet
#*****************************************************************************#

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "glmnet"

#=============================================================================#
# MRS + CSF
#=============================================================================#

# Actual training
which(colnames(X_train) == "ChrAge")
set.seed(123)
fit <- train(x = X_train[,-18],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbio.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)

#=============================================================================#
# CSF + age
#=============================================================================#

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbioage.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)

#=============================================================================#
# CSF only
#=============================================================================#

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbioonly.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)

#=============================================================================#
# CSF + MRS w/o epi-age
#=============================================================================#

# Actual training
set.seed(123)
which(colnames(X_train) == "EpiAge")
fit_woa <- train(x = X_train[,-c(10,18)],
                 y = Y_train$Y,
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)

# Save model
save(fit_woa, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_EN_CSFbio_noage.RData")

# Prediction in test set
testPred <- predict(fit_woa, X_test, type = "prob")
roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)


#*****************************************************************************#
# sPLS
#*****************************************************************************#

# Number of component (K)
K_CV <- 1:10

# Thresholding parameter (eta)
eta_CV <- seq(0.1,0.9,length.out = 20)

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = 0.5

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "spls"

#=============================================================================#
# MRS + CSF
#=============================================================================#

# Actual training
which(colnames(X_train) == "ChrAge")
set.seed(123)
fit <- train(x = X_train[,-18],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbio.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#=============================================================================#
# CSF + age
#=============================================================================#

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbioage.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)


#=============================================================================#
# CSF only
#=============================================================================#

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbioonly.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)

#=============================================================================#
# CSF + MRS w/o epi-age
#=============================================================================#

# Actual training
set.seed(123)
fit_woa <- train(x = X_train[,-c(10,18)],
                 y = Y_train$Y,
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)

# Save model
save(fit_woa, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbio_noage.RData")

# Prediction in test set
testPred <- predict(fit_woa, X_test, type = "prob")
roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)


#*****************************************************************************#
# Random Forest with recursive feature elimination
#*****************************************************************************#

#=============================================================================#
# CSF + MRS
#=============================================================================#

which(colnames(X_train) == "ChrAge")
X_train1 <- X_train[,-18]

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train1)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train1)){
  X_temp <- X_train1[,features]
  
  # Number of randomly selected predictors
  mtry_CV <- 1:length(features)
  
  # split rule
  splitrule_CV <- "gini"
  
  # minimal node size
  min.node.size_CV = 1:length(features)
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # Use MSE as performance metric
  performance_metric = "ROC"
  MLmethod = "ranger"
  
  # Actual training
  set.seed(456)
  fitList[[f]] <- train(x = X_temp,
                        y = Y_train$Y,
                        metric= performance_metric,
                        method = MLmethod,
                        tuneGrid = parameterGrid,
                        trControl = fitControl,
                        maximize = TRUE,
                        importance = "impurity")
  
  # ROC in CV
  fit <- fitList[[f]]
  performance[f] <- fit$results[(fit$results$mtry == fit$bestTune$mtry) &
                                  (fit$results$splitrule == fit$bestTune$splitrule) &
                                  (fit$results$min.node.size == fit$bestTune$min.node.size),"ROC"]
  
  # Variable importance: Gini index
  importance <- varImp(fit)$importance
  
  # Remove feature with lowest importance
  removedFeatures <- c(removedFeatures,rownames(importance)[which.min(importance$Overall)])
  
  # Selected features
  features <- setdiff(features,rownames(importance)[which.min(importance$Overall)])
}

# Get model with maximum performance
fit <- fitList[[which.max(performance)]]

# Save model
save(fit, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbio.RData")


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#=============================================================================#
# CSF + age
#=============================================================================#

X_train1 <- X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")]

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train1)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train1)){
  X_temp <- X_train1[,features]
  
  # Number of randomly selected predictors
  mtry_CV <- 1:length(features)
  
  # split rule
  splitrule_CV <- "gini"
  
  # minimal node size
  min.node.size_CV = 1:length(features)
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # Use MSE as performance metric
  performance_metric = "ROC"
  MLmethod = "ranger"
  
  # Actual training
  set.seed(123)
  fitList[[f]] <- train(x = X_temp,
                        y = Y_train$Y,
                        metric= performance_metric,
                        method = MLmethod,
                        tuneGrid = parameterGrid,
                        trControl = fitControl,
                        maximize = TRUE,
                        importance = "impurity")
  
  # ROC in CV
  fit <- fitList[[f]]
  performance[f] <- fit$results[(fit$results$mtry == fit$bestTune$mtry) &
                                  (fit$results$splitrule == fit$bestTune$splitrule) &
                                  (fit$results$min.node.size == fit$bestTune$min.node.size),"ROC"]
  
  # Variable importance: Gini index
  importance <- varImp(fit)$importance
  
  # Remove feature with lowest importance
  removedFeatures <- c(removedFeatures,rownames(importance)[which.min(importance$Overall)])
  
  # Selected features
  features <- setdiff(features,rownames(importance)[which.min(importance$Overall)])
}

# Get model with maximum performance
fit_wo <- fitList[[which.max(performance)]]

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbioage.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)


#=============================================================================#
# CSF only
#=============================================================================#

X_train1 <- X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")]

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train1)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train1)){
  X_temp <- X_train1[,features]
  
  # Number of randomly selected predictors
  mtry_CV <- 1:length(features)
  
  # split rule
  splitrule_CV <- "gini"
  
  # minimal node size
  min.node.size_CV = 1:length(features)
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # Use MSE as performance metric
  performance_metric = "ROC"
  MLmethod = "ranger"
  
  # Actual training
  set.seed(123)
  fitList[[f]] <- train(x = X_temp,
                        y = Y_train$Y,
                        metric= performance_metric,
                        method = MLmethod,
                        tuneGrid = parameterGrid,
                        trControl = fitControl,
                        maximize = TRUE,
                        importance = "impurity")
  
  # ROC in CV
  fit <- fitList[[f]]
  performance[f] <- fit$results[(fit$results$mtry == fit$bestTune$mtry) &
                                  (fit$results$splitrule == fit$bestTune$splitrule) &
                                  (fit$results$min.node.size == fit$bestTune$min.node.size),"ROC"]
  
  # Variable importance: Gini index
  importance <- varImp(fit)$importance
  
  # Remove feature with lowest importance
  removedFeatures <- c(removedFeatures,rownames(importance)[which.min(importance$Overall)])
  
  # Selected features
  features <- setdiff(features,rownames(importance)[which.min(importance$Overall)])
}

# Get model with maximum performance
fit_wo <- fitList[[which.max(performance)]]

# Save model
save(fit_wo, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbioonly.RData")

# Prediction in test set
testPred <- predict(fit_wo, X_test, type = "prob")
roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)


#=============================================================================#
# CSF + MRS w/o epi-age
#=============================================================================#

X_train1 <- X_train[,-c(10,18)]

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train1)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train1)){
  X_temp <- X_train1[,features]
  
  # Number of randomly selected predictors
  mtry_CV <- 1:length(features)
  
  # split rule
  splitrule_CV <- "gini"
  
  # minimal node size
  min.node.size_CV = 1:length(features)
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # Use MSE as performance metric
  performance_metric = "ROC"
  MLmethod = "ranger"
  
  # Actual training
  set.seed(123)
  fitList[[f]] <- train(x = X_temp,
                        y = Y_train$Y,
                        metric= performance_metric,
                        method = MLmethod,
                        tuneGrid = parameterGrid,
                        trControl = fitControl,
                        maximize = TRUE,
                        importance = "impurity")
  
  # ROC in CV
  fit <- fitList[[f]]
  performance[f] <- fit$results[(fit$results$mtry == fit$bestTune$mtry) &
                                  (fit$results$splitrule == fit$bestTune$splitrule) &
                                  (fit$results$min.node.size == fit$bestTune$min.node.size),"ROC"]
  
  # Variable importance: Gini index
  importance <- varImp(fit)$importance
  
  # Remove feature with lowest importance
  removedFeatures <- c(removedFeatures,rownames(importance)[which.min(importance$Overall)])
  
  # Selected features
  features <- setdiff(features,rownames(importance)[which.min(importance$Overall)])
}

# Get model with maximum performance
fit_woa <- fitList[[which.max(performance)]]

# Save model
save(fit_woa, file = "Models/EMIF_Models/CSF/Fit_EMIF_MCI_RF_CSFbio_noage.RData")

# Prediction in test set
testPred <- predict(fit_woa, X_test, type = "prob")
roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)



#*****************************************************************************#
# Plotting
#*****************************************************************************#

# Load models and prediction:

# Random forest models
load("EMIF/Fit_EMIF_MCI_RF_CSFbio.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_RF_CSFbioonly.RData")
RF_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_RF_CSFbio_noage.RData")
RF_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_RF_CSFbioage.RData")
RF_wa <- predict(fit_wo, X_test, type = "prob")

# ElasticNet models
load("EMIF/Fit_EMIF_MCI_EN_CSFbio.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_EN_CSFbioonly.RData")
EN_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_EN_CSFbio_noage.RData")
EN_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_EN_CSFbioage.RData")
EN_wa <- predict(fit_wo, X_test, type = "prob")

# sPLS-DA models
load("EMIF/Fit_EMIF_MCI_sPLS_CSFbio.RData")
sPLS <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_sPLS_CSFbioonly.RData")
sPLS_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_sPLS_CSFbio_noage.RData")
sPLS_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_sPLS_CSFbioage.RData")
sPLS_wa <- predict(fit_wo, X_test, type = "prob")

# Combine predictions into dataframe
testDF <- data.frame(EN = EN$MCI,
                     sPLS = sPLS$MCI,
                     RF = RF$MCI,
                     EN_wo = EN_wo$MCI,
                     sPLS_wo = sPLS_wo$MCI,
                     RF_wo = RF_wo$MCI,
                     EN_woa = EN_woa$MCI,
                     sPLS_woa = sPLS_woa$MCI,
                     RF_woa = RF_woa$MCI,
                     EN_wa = EN_wa$MCI,
                     sPLS_wa = sPLS_wa$MCI,
                     RF_wa = RF_wa$MCI,
                     EpiAge = X_test$EpiAge,
                     Age = Y_test$Age,
                     Sex= Y_test$Gender,
                     Diagnosis = Y_test$Diagnosis,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "MCI",1,0))



# Prepare data for plotting
score <- c("EN_wo", "sPLS_wo", "RF_wo", "EN", "sPLS", "RF")
scoreName <- c("CSF (EN)", "CSF (sPLS-DA)", "CSF (RF)",
               "CSF + MRS (EN)", "CSF + MRS (sPLS-DA)", "CSF + MRS (RF)")
scoreName1 <- c("CSF (EN)\t\t", "CSF (sPLS-DA)\t","CSF (RF-RFE)\t",
                "MRSs/CSF (EN):\t", "MRSs/CSF (sPLS-DA):", 
                "MRSs/CSF (RF-RFE):")
plotDF <- as.data.frame(testDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
liValue <- rep(NA, length(score))
uiValue <- rep(NA, length(score))

# Get sensitivities and specificities for ROC curve
for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName1[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(test)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(test)[3]),2),nsmall = 2)
}


plotAUC <- data.frame(AUC = paste0(scoreName1,"\t",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName1,
                      X = 0.75,
                      Y = rev(seq(0.05,0.3,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName1)

# Colors for plotting
colors <- c(rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  #ggtitle("MCI vs Control") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "EMIF-AD/ModelPerformance/ROC_MCI_EMIF_combinedPlot.png", width = 7, height = 5)

