
# Load packages
library(prospectr)
library(tidyverse)
library(caret)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_EMIF.RData")
load("EMIF/predictedScore_factors_EMIF.RData")

###############################################################################

# Filter data

###############################################################################

metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X
metaData_EMIF <- metaData_EMIF[rownames(predictedScore_factors),]

# Keep midlife samples only
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,]
range(metaData_EMIF$Age)

# Remove converters (Control to MCI or AD)

# CTR CONVERT
converters1 <-  unique(rownames(metaData_EMIF)[metaData_EMIF$CTR_Convert == 1])[-1]

# CTR to MCI or AD
converters2 <- unique(metaData_EMIF$X[(metaData_EMIF$LastFU_Diagnosis == "MCI") | (metaData_EMIF$LastFU_Diagnosis == "AD")])[-1]
converters2 <- intersect(converters2, metaData_EMIF$X[metaData_EMIF$Diagnosis == "NL"])

# All converters
converters_all <- unique(c(converters1, converters2))
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF), converters_all),]
table(metaData_EMIF$CTR_Convert)

# Save meta data
samples <- intersect(metaData_EMIF$X, rownames(predictedScore_factors))
metaData_fil <- metaData_EMIF[samples,]
save(metaData_fil, file = "EMIF/metaData_fil.RData")

# save predicted scores
predictedScore_factors_fil <- predictedScore_factors[samples,]
save(predictedScore_factors_fil, file = "EMIF/predictedScore_factors_fil.RData")


###############################################################################

# Data splitting with Kennard-Stone algorithm

###############################################################################

# Load packages
library(prospectr)
library(tidyverse)
library(caret)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_fil.RData")
load("EMIF/predictedScore_factors_fil.RData")

# Male
nTrain0 <- round(0.7*nrow(predictedScore_factors_fil[metaData_fil$Gender == 0,]))
selectedSamples0 <- prospectr::kenStone(
  X = predictedScore_factors_fil[metaData_fil$Gender == 0,], 
  k = nTrain0,
  .center = TRUE,
  .scale = TRUE
)
temp0 <- rownames(predictedScore_factors_fil)[metaData_fil$Gender == 0]
test0 <- temp0[selectedSamples0$test]
train0 <- temp0[selectedSamples0$model]

# Female
nTrain1 <- round(0.7*nrow(predictedScore_factors_fil[metaData_fil$Gender == 1,]))
selectedSamples1 <- prospectr::kenStone(
  X = predictedScore_factors_fil[metaData_fil$Gender == 1,], 
  k = nTrain1,
  .center = TRUE,
  .scale = TRUE
)
temp1 <- rownames(predictedScore_factors_fil)[metaData_fil$Gender == 1]
test1 <- temp1[selectedSamples1$test]
train1 <- temp1[selectedSamples1$model]


# Combine into single train and test set
X_train <- predictedScore_factors_fil[c(train0, train1),]
Y_train <- metaData_fil[c(train0, train1),]

X_test <- predictedScore_factors_fil[c(test0, test1),]
Y_test <- metaData_fil[c(test0, test1),]

# No overlap in test and training?
intersect(rownames(X_train), rownames(X_test))

# No duplicated samples?
sum(duplicated(X_test))
sum(duplicated(X_train))

# Distribution of diagnosis
table(Y_test$Diagnosis)
table(Y_train$Diagnosis)

save(X_train, file = "EMIF/X_train_EMIF.RData")
save(Y_train, file = "EMIF/Y_train_EMIF.RData")
save(X_test, file = "EMIF/X_test_EMIF.RData")
save(Y_test, file = "EMIF/Y_test_EMIF.RData")

###############################################################################

# Make predictions: MCI vs Control

###############################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 

load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

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

# Create cross-validation indices
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
save(fit, file = "EMIF/Fit_EMIF_MCI_EN.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)
ci(roc_test)
varImp(fit)

roc_test <- roc(Y_test$Y, X_test$EpiAge)
auc(roc_test)

roc_test <- roc(Y_test$Y, Y_test$AB_Zscore)
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
save(fit, file = "EMIF/Fit_EMIF_MCI_sPLS.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)



#*****************************************************************************#
# RandomForest + recursive feature elimination
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train)){
  X_temp <- X_train[,features]
  
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
save(fitList, performance, removedFeatures, file = "EMIF/FitList_EMIF_MCI_RF.RData")

# Get model with maximum performance
fit <- fitList[[which.max(performance)]]

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_RF.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#*****************************************************************************#
# Make ROC plot
#*****************************************************************************#

# Load models
load("EMIF/Fit_EMIF_MCI_RF.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_EN.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_sPLS.RData")
sPLS <- predict(fit, X_test, type = "prob")

# Combine values in data frame
plotDF <- data.frame(EN = log(EN$MCI/(1-EN$MCI)),
                     sPLS = log(sPLS$MCI/(1-sPLS$MCI)),
                     RF = log(RF$MCI/(1-RF$MCI)),
                     EpiAge = X_test$EpiAge,
                     Age = Y_test$Age,
                     Sex = Y_test$Gender,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "MCI",1,0))


# Check whether the predicted score is significantly associated with MCI status:

# ElasticNet
model <- lm(Ynum ~ EN + Age + Sex, data = plotDF)
summary(model)

# sPLS-DA
model <- lm(Ynum ~ EN + Age + Sex, data = plotDF)
summary(model)

# Random Forest
model <- lm(Ynum ~ EN + Age + Sex, data = plotDF)
summary(model)

# Epigenetic age
model <- lm(Ynum ~ EpiAge + Age + Sex, data = plotDF)
summary(model)

# Calculate sensitivites, specificities and AUC values
score <- c("EN", "sPLS", "RF","EpiAge")
scoreName <- c("ElasticNet", "sPLS-DA", 
               "Random Forest", "Epi-Age")
ROCplot <- NULL                       # Data frame with sensitivities and specificities
aucValue <- rep(NA, length(score))    # AUC
liValue <- rep(NA, length(score))     # lower interval value of AUC
uiValue <- rep(NA, length(score))     # upper interval value of AUC

for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(test)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(test)[3]),2),nsmall = 2)
}

# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0("AUROC = ",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName,
                      X = 0.8,
                      Y = rev(seq(0.05,0.25,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

# Colors for plotting
colors <- rev(c("#084594","#EF3B2C","#CB181D", "#99000D"))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold", size = 4) +
  scale_color_manual(values = colors) +
  ggtitle("MCI vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF.png", width = 7.5, height = 5)

###############################################################################

# Make predictions: AD vs Control

###############################################################################

library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 

load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

# AD and NL only
X_train <- X_train[(Y_train$Diagnosis == "AD") | (Y_train$Diagnosis == "NL"),]
Y_train <- Y_train[(Y_train$Diagnosis == "AD")| (Y_train$Diagnosis == "NL"),]
X_test <- X_test[(Y_test$Diagnosis == "AD") |  (Y_test$Diagnosis == "NL"),]
Y_test <- Y_test[(Y_test$Diagnosis == "AD") | (Y_test$Diagnosis == "NL"),]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","AD"),
                    levels = c("Control", "AD"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","AD"),
                   levels = c("Control", "AD"))

table(Y_test$Y)
table(Y_train$Y)

# Create cross-validation indices
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
save(fit, file = "EMIF/Fit_EMIF_AD_EN.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)
ci(roc_test)
varImp(fit)

roc_test <- roc(Y_test$Y, X_test$EpiAge)
auc(roc_test)

roc_test <- roc(Y_test$Y, Y_test$AB_Zscore)
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
save(fit, file = "EMIF/Fit_EMIF_AD_sPLS.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)



#*****************************************************************************#
# RandomForest - recursive feature elimination
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

removedFeatures <- NULL
fitList <- list()
features <- colnames(X_train)
performance <- rep(NA,length(features))
for (f in 1:ncol(X_train)){
  X_temp <- X_train[,features]
  
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
save(fitList, performance, removedFeatures, file = "EMIF/FitList_EMIF_AD_RF.RData")

# Get model with maximum performance
fit <- fitList[[which.max(performance)]]

# Save model
save(fit, file = "EMIF/Fit_EMIF_AD_RF.RData")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")
roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)


#*****************************************************************************#
# Make ROC plot
#*****************************************************************************#

# Load models
load("EMIF/Fit_EMIF_AD_RF.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_AD_EN.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_AD_sPLS.RData")
sPLS <- predict(fit, X_test, type = "prob")

# Combine values in data frame
plotDF <- data.frame(EN = log(EN$AD/(1-EN$AD)),
                     sPLS = log(sPLS$AD/(1-sPLS$AD)),
                     RF = log(RF$AD/(1-RF$AD)),
                     EpiAge = X_test$EpiAge,
                     Age = Y_test$Age,
                     Sex = Y_test$Gender,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "AD",1,0))


# Check whether the predicted score is significantly associated with MCI status:

# ElasticNet
model <- lm(Ynum ~ EN + Age + Sex, data = plotDF)
summary(model)

# sPLS-DA
model <- lm(Ynum ~ sPLS + Age + Sex, data = plotDF)
summary(model)

# Random Forest
model <- lm(Ynum ~ RF + Age + Sex, data = plotDF)
summary(model)

# Epigenetic age
model <- lm(Ynum ~ EpiAge + Age + Sex, data = plotDF)
summary(model)

# Calculate sensitivites, specificities and AUC values
score <- c("EN", "sPLS", "RF","EpiAge")
scoreName <- c("ElasticNet", "sPLS-DA", 
               "Random Forest", "Epi-Age")
ROCplot <- NULL                       # Data frame with sensitivities and specificities
aucValue <- rep(NA, length(score))    # AUC
liValue <- rep(NA, length(score))     # lower interval value of AUC
uiValue <- rep(NA, length(score))     # upper interval value of AUC

for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(test)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(test)[3]),2),nsmall = 2)
}

# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0("AUROC = ",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName,
                      X = 0.8,
                      Y = rev(seq(0.05,0.25,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

# Colors for plotting
colors <- rev(c("#084594","#EF3B2C","#CB181D", "#99000D"))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold", size = 4) +
  scale_color_manual(values = colors) +
  ggtitle("AD vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_AD_EMIF.png", width = 7.5, height = 5)
