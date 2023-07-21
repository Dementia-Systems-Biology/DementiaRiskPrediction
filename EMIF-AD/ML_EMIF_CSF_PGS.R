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

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load training and test data
load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

# Add PGSs
load("EMIF/metaData_fil.RData")
load("~/EMIF/PGS_EMIF_AD.RData")
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


# MCI and NL only
X_train <- X_train[(Y_train$Diagnosis == "MCI") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "MCI")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "MCI") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "MCI") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))


Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]

X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])

X_train <- cbind.data.frame(X_train[rownames(Y_train),],output[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], output[rownames(Y_test),])




table(Y_test$Y)
table(Y_train$Y)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Settings for repeated cross-valBasename
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
# MRS + CSF + PGS
#=============================================================================#

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_EN_CSF_PGS.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)

# CSF + PGS only

# Actual training
set.seed(123)
fit <- train(x = cbind.data.frame(X_train[,str_detect(colnames(X_train), "PGS")],
                                  X_train[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore")]),
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_EN_CSF_PGSonly.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)




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
# MRS + CSF + PGS
#=============================================================================#

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_sPLS_CSF_PGS.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)


# CSF + PGS only

# Actual training
set.seed(123)
fit <- train(x = cbind.data.frame(X_train[,str_detect(colnames(X_train), "PGS")],
                                  X_train[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore")]),
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_sPLS_CSF_PGSonly.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)



#*****************************************************************************#
# RandomForest
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

#=============================================================================#
# CSF + MRS + PGS
#=============================================================================#


X_train1 <- X_train

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
save(fit, file = "EMIF/Fit_EMIF_MCI_RF_CSF_PGS.RData")


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#=============================================================================#
# CSF + PGS
#=============================================================================#


X_train1 <- cbind.data.frame(X_train[,str_detect(colnames(X_train), "PGS")],
                             X_train[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore")])
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
save(fit, file = "EMIF/Fit_EMIF_MCI_RF_CSF_PGSonly.RData")


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#*****************************************************************************#
# Plotting
#*****************************************************************************#
load("EMIF/Fit_EMIF_MCI_EN_CSFbio.RData")
testPred <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_EN_CSFbioonly.RData")
testPred_wo <- predict(fit_wo, X_test, type = "prob")

roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
roc_test_wo <- pROC::roc(Y_test$Y, testPred_wo$MCI)
pROC::roc.test(roc_test, roc_test_wo)


load("EMIF/Fit_EMIF_MCI_RF_CSF_PGS.RData")
RF <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_EN_CSF_PGS.RData")
EN <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_sPLS_CSF_PGS.RData")
sPLS <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_RF_CSF_PGSonly.RData")
RFonly <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_EN_CSF_PGSonly.RData")
ENonly <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_MCI_sPLS_CSF_PGSonly.RData")
sPLSonly <- predict(fit, X_test, type = "prob")



plotDF <- data.frame(EN = EN$MCI,
                     sPLS = sPLS$MCI,
                     RF = RF$MCI,
                     ENonly = ENonly$MCI,
                     sPLSonly = sPLSonly$MCI,
                     RFonly = RFonly$MCI,
                     EpiAge = X_test$EpiAge,
                     Age = Y_test$Age,
                     Sex= Y_test$Gender,
                     Diagnosis = Y_test$Diagnosis,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "MCI",1,0))


# Calculate sensitivites, specificities and AUC values
score <- c("ENonly", "sPLSonly", "RFonly","EN", "sPLS", "RF")
scoreName1 <- c("PGSs/CSF (EN)\t\t\t",
                "PGSs/CSF (sPLS-DA)\t\t","PGSs/CSF (RF-RFE)\t\t",
                "MRSs/PGSs/CSF (EN):\t\t", "MRSs/PGSs/CSF (sPLS-DA):", 
                "MRSs/PGSs/CSF (RF-RFE):\t")
ROCplot <- NULL                       # Data frame with sensitivities and specificities
aucValue <- rep(NA, length(score))    # AUC
liValue <- rep(NA, length(score))     # lower interval value of AUC
uiValue <- rep(NA, length(score))     # upper interval value of AUC

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

# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0(scoreName1,"\t",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName1,
                      X = 0.7,
                      Y = rev(seq(0.05,0.3,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName1)

# Colors for plotting
colors <- rev(c("#084594","#EF3B2C","#CB181D", "#99000D"))

colors <- c(rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold", size = 4) +
  scale_color_manual(values = colors) +
  #ggtitle("MCI vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_CSF_PGSs.png", width = 7, height = 5)
save(p, file = "EMIF/ROC_MCI_EMIF_MRS_CSF_APOE.RData")