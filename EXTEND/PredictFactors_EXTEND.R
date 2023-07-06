# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("Data/metaData_ageFil.RData")
load("Data/methSet_allNorm_fil.RData")


#################################################################################

# Prepare phenotype data (Y)

################################################################################

#******************************************************************************#
# Discrete variables
#******************************************************************************#

# Low Education
Education <- ifelse(dat$None.of.the.above == 1,"Yes","No")

# Physical Inactivity
Physical <- ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 1,"No","Yes")

# Unhealthy Diet
Diet <- ifelse(as.numeric(dat$Fruit) + as.numeric(dat$Vegtables) > 3, "No", "Yes")


# Depression
Depression <- ifelse(dat$Depression==1,"Yes","No")

# Diabetes
Diabetes <-ifelse(dat$T2.Diabetes==1,"Yes","No")

# Heart Disease
HeartDisease <- ifelse(dat$Heart.Disease==1,"Yes","No")

# Sex
SexMale <- ifelse(dat$Sex == 1, "Yes", "No")



#******************************************************************************#
# Continuous variables
#******************************************************************************#

# High Systolic Blood Pressure
SysBP <- as.numeric(dat$MeanSysBP)

# High total Cholesterol (log transformed)
TotalChol <- as.numeric(dat$Chol)


# Combine
Y_all <- data.frame(SysBP,
                    TotalChol,
                    Education,
                    Physical,
                    Diet,
                    Depression,
                    Diabetes,
                    HeartDisease,
                    SexMale)

rownames(Y_all) <- dat$Basename


Y_train <- Y_all[colnames(methSet_allNorm_fil),] 

save(Y_train, file = "EXTEND/Y_train_factors.RData")

#################################################################################

# Prepare correlation data

################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("Data/methSet_allNorm_fil.RData")

X_train <- methSet_allNorm_fil

all(colnames(X_train) == rownames(Y_train))


#*****************************************************************************#
# correlation-based selection
#*****************************************************************************#

correlations <- matrix(NA, nrow = nrow(X_train), ncol = ncol(Y_train))
for (i in 1:ncol(Y_train)) {
  factor <- Y_train[!is.na(Y_train[,i]),i]
  if(!is.numeric(Y_train[,i])){
    factor <- ifelse(factor == "Yes",1,0)
  }
  dataMatrix <- X_train[,!is.na(Y_train[,i])]
  
  correlations[,i] <- apply(dataMatrix, 1, 
                            function(x){cor(x, 
                                            factor, 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_train)
colnames(correlations) <- colnames(Y_train)
save(correlations, file = "EXTEND/correlations_factors.RData")


# Selected probes
selectedProbes <- list()
for (i in 1:ncol(correlations)){
  selectedProbes[[i]] <- names(tail(sort(abs(correlations[,i])),70000))
}
names(selectedProbes) <- colnames(correlations)
save(selectedProbes, file = "EXTEND/selectedProbes.RData")



################################################################################

# Model Training

################################################################################

#*****************************************************************************#
# Correlation + ElasticNet
#*****************************************************************************#
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("EXTEND/selectedProbes.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil

# Convert to M-values
X_train <- log2(X_train/(1-X_train))

# Test if samples are in correct order
all(colnames(X_train) == rownames(Y_train))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Machine learning method
MLmethod = "glmnet"

# Model training
trainList <- list()
predList <- list()
for (f in 1:ncol(Y_train)){
  
  
  #=============================================================================#
  # Continuous variables
  #=============================================================================#
  
  if (is.numeric(Y_train[,f])){
    
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "MAE"
    
    # Risk factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    # Model training
    trainResults <- list()
    predResults <- list()
    for(i in 1:length(CVindex)) {
      
      # Sample indices for repeated CV
      index <- list(CVindex[[i]])
      
      # Select samples from specific fold
      X_CV <- dataMatrix[,index[[1]]]
      factor_CV_cor <- riskfactor[index[[1]]]
      
      
      # Calculate correlations with factors
      correlations_CV <- apply(X_CV, 1, 
                               function(x){cor(x, factor_CV_cor,method = "spearman")})
      
      names(correlations_CV) <- rownames(X_CV)
      
      # Select top correlated features for each factor
      finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = TRUE,
                                 summaryFunction = regressionSummary,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(dataMatrix[finalProbes,]),
                   y = riskfactor,
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = FALSE)
      
      
      trainResults[[i]] <- fit$results
      predResults[[i]] <- fit$pred
      
    }
    
  } 
  
  #=============================================================================#
  # Discrete variables
  #=============================================================================#
  if (!is.numeric(Y_train[,f])){
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    # Model training
    trainResults <- list()
    for(i in 1:length(CVindex)) {
      
      # Sample indices for repeated CV
      index <- list(CVindex[[i]])
      
      # Select samples from specific fold
      X_CV <- dataMatrix[,index[[1]]]
      factor_CV_cor <- ifelse(riskfactor[index[[1]]] == "Yes",1,0)
      
      
      # Calculate correlations with factors
      correlations_CV <- apply(X_CV, 1, 
                               function(x){cor(x, factor_CV_cor,method = "spearman")})
      
      names(correlations_CV) <- rownames(X_CV)
      
      # Select top correlated features for each factor
      finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = TRUE,
                                 summaryFunction = twoClassSummary,
                                 classProbs = TRUE,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(dataMatrix[finalProbes,]),
                   y = factor(riskfactor, levels = c("No", "Yes")),
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = TRUE)
      
      
      trainResults[[i]] <- fit$results
      predResults[[i]] <- fit$pred
      
    }
    
  }
  
  
  # Save output
  trainList[[f]] <- trainResults
  predList[[f]] <- predResults
}

names(trainList) <- colnames(Y_train)
names(predList) <- colnames(Y_train)

# Save results
save(trainList, predList, file = "EXTEND/OutputList_Cor_EN.RData")


#load("EXTEND/OutputList_Cor_EN.RData")

# Get optimal hyperparameters
PerformanceList <- list()
for (i in 1:length(trainList)){
  if (!is.numeric(Y_train[,i])){
    trainResults <- trainList[[i]]
    perf <- matrix(NA, nrow = 1000, ncol = 25)
    for (j in 1:length(trainResults)){
      perf[,j] <- trainResults[[j]]$ROC
    }
    optPar <- which.max(rowMeans(perf))
    optPerf <- perf[optPar,]
    optAlpha <- trainResults[[1]]$alpha[optPar]
    optLambda <- trainResults[[1]]$lambda[optPar]
    
    PerformanceList[[i]] <- list(optPerf, optAlpha, optLambda)
  }
  if (is.numeric(Y_train[,i])){
    trainResults <- trainList[[i]]
    perf <- matrix(NA, nrow = 1000, ncol = 25)
    for (j in 1:length(trainResults)){
      perf[,j] <- trainResults[[j]]$MAE
    }
    optPar <- which.min(rowMeans(perf))
    optPerf <- perf[optPar,]
    optAlpha <- trainResults[[1]]$alpha[optPar]
    optLambda <- trainResults[[1]]$lambda[optPar]
    
    PerformanceList[[i]] <- list(optPerf, optAlpha, optLambda)
  }
  
}
names(PerformanceList) <- names(trainList)
save(PerformanceList, file = "EXTEND/PerformanceList_Cor_EN.RData")



# Train model with optimal hyperparameters
bestModel <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for lambda
  lambdaCV <- PerformanceList[[f]][[3]]
  
  # Set grid for alpha
  alphaCV <- PerformanceList[[f]][[2]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  colnames(parameterGrid) <- c(".alpha", ".lambda")
  
  # ML method
  MLmethod = "glmnet"
  
  if (is.numeric(Y_train[,f])){
    # Performance metric
    performance_metric = "MAE"
    
    # factor
    factor_cor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Calculate correlations with factors
    correlations_all <- apply(dataMatrix, 1, 
                              function(x){cor(x, factor_cor,method = "spearman")})
    
    names(correlations_all) <- rownames(dataMatrix)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_all)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none")
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix[finalProbes,]),
                            y = factor_cor,
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = FALSE)
    
    
  }
  if (!is.numeric(Y_train[,f])){
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    factor_cor <- ifelse(riskfactor == "Yes",1,0)
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Calculate correlations with factors
    correlations_all <- apply(dataMatrix, 1, 
                              function(x){cor(x, factor_cor,method = "spearman")})
    
    names(correlations_all) <- rownames(dataMatrix)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_all)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none", 
                               classProbs = TRUE)
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix[finalProbes,]),
                            y = factor(riskfactor, levels = c("No", "Yes")),
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = TRUE)
    
    
  }
  
}
names(bestModel) <- colnames(Y_train)

# Save best models
save(bestModel, file = "EXTEND/bestModel_Cor_EN.RData")


# Pred vs Obs in CV
perf_discrete <- rep(NA,7)
for (i in 3:9){
  
  predResults <- NULL
  for (j in 1:length(predList[[i]])){
    test <- predList[[i]][[j]]
    test <- test[(test$alpha == PerformanceList[[i]][[2]]) &
                   (test$lambda == PerformanceList[[i]][[3]]),c("obs", "Yes")]
    predResults <- rbind.data.frame(predResults,test)
  }
  perf_discrete[i-2] <- as.numeric(pROC::auc(pROC::roc(predResults$obs, predResults$Yes)))
}

AUCdf <- data.frame(RiskFactor = colnames(Y_train)[3:9],
                    AUC = perf_discrete)
save(AUCdf, file = "EXTEND/AUCdf_Cor_EN.RData")

# Pred vs Obs in CV
perf_continuous <- rep(NA,2)
for (i in 1:2){
  
  predResults <- NULL
  for (j in 1:length(predList[[i]])){
    test <- predList[[i]][[j]]
    test <- test[(test$alpha == PerformanceList[[i]][[2]]) &
                   (test$lambda == PerformanceList[[i]][[3]]),c("obs", "pred")]
    predResults <- rbind.data.frame(predResults,test)
  }
  perf_continuous[i] <- as.numeric(R2(pred = predResults$pred, obs = predResults$obs))
}

R2df <- data.frame(RiskFactor = colnames(Y_train)[1:2],
                   R2 = perf_continuous)
save(R2df, file = "EXTEND/R2df_Cor_EN.RData")



#*****************************************************************************#
# Correlation + Random forest
#*****************************************************************************#
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(ranger)
library(caret)
library(tidyverse)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("EXTEND/selectedProbes.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil

# Convert to M-values
X_train <- log2(X_train/(1-X_train))
rm(methSet_allNorm_fil)

# Test if samples are in correct order
all(colnames(X_train) == rownames(Y_train))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Machine learning method
MLmethod = "ranger"

# Model training
trainList <- list()
predList <- list()
for (f in 7:ncol(Y_train)){
  
  
  #=============================================================================#
  # Continuous variables
  #=============================================================================#
  
  if (is.numeric(Y_train[,f])){
    
    # Number of randomly selected predictors
    mtry_CV <- c(1000, 2000,3000,4000,5000,6000,7000,8000)
    
    # split rule
    splitrule_CV <- "variance"
    
    # minimal node size
    min.node.size_CV = c(10,20,30,40,50,60)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
    colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
    
    
    # Performance metric
    performance_metric = "MAE"
    
    # Risk factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    # Model training
    trainResults <- list()
    predResults <- list()
    for(i in 1:length(CVindex)) {
      
      # Sample indices for repeated CV
      index <- list(CVindex[[i]])
      
      # Select samples from specific fold
      X_CV <- dataMatrix[,index[[1]]]
      factor_CV_cor <- riskfactor[index[[1]]]
      
      
      # Calculate correlations with factors
      correlations_CV <- apply(X_CV, 1, 
                               function(x){cor(x, factor_CV_cor,method = "spearman")})
      
      names(correlations_CV) <- rownames(X_CV)
      
      # Select top correlated features for each factor
      finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = TRUE,
                                 summaryFunction = regressionSummary,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(dataMatrix[finalProbes,]),
                   y = riskfactor,
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = FALSE)
      
      
      trainResults[[i]] <- fit$results
      predResults[[i]] <- fit$pred
      
    }
    
    
  } 
  
  #=============================================================================#
  # Discrete variables
  #=============================================================================#
  if (!is.numeric(Y_train[,f])){
    # Number of randomly selected predictors
    mtry_CV <- c(100,500,1000,1500,2000,3000,4000)
    
    # split rule
    splitrule_CV <- "gini"
    
    # minimal node size
    min.node.size_CV = c(3,5,10,15,20)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
    colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    # Model training
    trainResults <- list()
    predResults <- list()
    for(i in 1:length(CVindex)) {
      
      # Sample indices for repeated CV
      index <- list(CVindex[[i]])
      
      # Select samples from specific fold
      X_CV <- dataMatrix[,index[[1]]]
      factor_CV_cor <- ifelse(riskfactor[index[[1]]] == "Yes",1,0)
      
      
      # Calculate correlations with factors
      correlations_CV <- apply(X_CV, 1, 
                               function(x){cor(x, factor_CV_cor,method = "spearman")})
      
      names(correlations_CV) <- rownames(X_CV)
      
      # Select top correlated features for each factor
      finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = TRUE,
                                 summaryFunction = twoClassSummary,
                                 classProbs = TRUE,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(dataMatrix[finalProbes,]),
                   y = factor(riskfactor, levels = c("No", "Yes")),
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = TRUE)
      
      
      trainResults[[i]] <- fit$results
      predResults[[i]] <- fit$pred
      
    }
    
  }
  
  # Save output
  trainList[[f]] <- trainResults
  predList[[f]] <- predResults
  
  save(trainList, predList, file = "EXTEND/OutputList_Cor_RF.RData")
  
}

names(trainList) <- colnames(Y_train)
names(predList) <- colnames(Y_train)

# Save results
save(trainList, predList, file = "EXTEND/OutputList_Cor_RF.RData")


#load("EXTEND/OutputList_Cor_EN.RData")

# Get optimal hyperparameters
PerformanceList <- list()
for (i in 1:length(trainList)){
  if (!is.numeric(Y_train[,i])){
    trainResults <- trainList[[i]]
    perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
    for (j in 1:length(trainResults)){
      perf[,j] <- trainResults[[j]]$ROC
    }
    optPar <- which.max(rowMeans(perf))
    optPerf <- perf[optPar,]
    opt_mtry <- trainResults[[1]]$mtry[optPar]
    opt_splitrule <- trainResults[[1]]$splitrule[optPar]
    opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
    
    PerformanceList[[i]] <- list(optPerf, opt_mtry, opt_splitrule, opt_min.node.size)
  }
  if (is.numeric(Y_train[,i])){
    trainResults <- trainList[[i]]
    perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
    for (j in 1:length(trainResults)){
      perf[,j] <- trainResults[[j]]$MAE
    }
    optPar <- which.min(rowMeans(perf))
    optPerf <- perf[optPar,]
    opt_mtry <- trainResults[[1]]$mtry[optPar]
    opt_splitrule <- trainResults[[1]]$splitrule[optPar]
    opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
    
    PerformanceList[[i]] <- list(optPerf, opt_mtry, opt_splitrule, opt_min.node.size)
  }
  
}
names(PerformanceList) <- names(trainList)
save(PerformanceList, file = "EXTEND/PerformanceList_Cor_RF.RData")



# Train model with optimal hyperparameters
bestModel <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for mtry
  mtry_CV <- PerformanceList[[f]][[2]]
  
  # Set grid for splitrule
  splitrule_CV <- PerformanceList[[f]][[3]]
  
  # Set grid for splitrule
  min.node.size_CV <- PerformanceList[[f]][[4]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # ML method
  MLmethod = "ranger"
  
  if (is.numeric(Y_train[,f])){
    # Performance metric
    performance_metric = "MAE"
    
    # factor
    factor_cor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Calculate correlations with factors
    correlations_all <- apply(dataMatrix, 1, 
                              function(x){cor(x, factor_cor,method = "spearman")})
    
    names(correlations_all) <- rownames(dataMatrix)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_all)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none")
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix[finalProbes,]),
                            y = factor_cor,
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = FALSE)
    
    
  }
  if (!is.numeric(Y_train[,f])){
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    factor_cor <- ifelse(riskfactor == "Yes",1,0)
    
    # dataMatrix
    dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Calculate correlations with factors
    correlations_all <- apply(dataMatrix, 1, 
                              function(x){cor(x, factor_cor,method = "spearman")})
    
    names(correlations_all) <- rownames(dataMatrix)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_all)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none", 
                               classProbs = TRUE)
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix[finalProbes,]),
                            y = factor(riskfactor, levels = c("No", "Yes")),
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = TRUE)
    
    
  }
  
}
names(bestModel) <- colnames(Y_train)

# Save best models
save(bestModel, file = "EXTEND/bestModel_Cor_RF.RData")

# Pred vs Obs in CV
perf_discrete <- rep(NA,7)
for (i in 3:9){
  
  predResults <- NULL
  for (j in 1:length(predList[[i]])){
    test <- predList[[i]][[j]]
    test <- test[(test$mtry == PerformanceList[[i]][[2]]) &
                   (test$splitrule == PerformanceList[[i]][[3]]) &
                   (test$min.node.size == PerformanceList[[i]][[4]]),c("obs", "Yes")]
    predResults <- rbind.data.frame(predResults,test)
  }
  perf_discrete[i-2] <- as.numeric(pROC::auc(pROC::roc(predResults$obs, predResults$Yes)))
}

AUCdf <- data.frame(RiskFactor = colnames(Y_train)[3:9],
                    AUC = perf_discrete)
save(AUCdf, file = "EXTEND/AUCdf_Cor_RF.RData")

# Pred vs Obs in CV
perf_continuous <- rep(NA,2)
for (i in 1:2){
  
  predResults <- NULL
  for (j in 1:length(predList[[i]])){
    test <- predList[[i]][[j]]
    test <- test[(test$mtry == PerformanceList[[i]][[2]]) &
                   (test$splitrule == PerformanceList[[i]][[3]]) &
                   (test$min.node.size == PerformanceList[[i]][[4]]),c("obs", "pred")]
    predResults <- rbind.data.frame(predResults,test)
  }
  perf_continuous[i] <- as.numeric(R2(pred = predResults$pred, obs = predResults$obs))
}

R2df <- data.frame(RiskFactor = colnames(Y_train)[1:2],
                   R2 = perf_continuous)
save(R2df, file = "EXTEND/R2df_Cor_RF.RData")

#*****************************************************************************#
# ElasticNet w/o feature selection
#*****************************************************************************#
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil

# Convert to M-values
X_train <- log2(X_train/(1-X_train))

# Test if samples are in correct order
all(colnames(X_train) == rownames(Y_train))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Machine learning method
MLmethod = "glmnet"

# Model training
trainList <- list()
predList <- list()
PerformanceList <- list()
for (f in 7:ncol(Y_train)){
  
  
  #=============================================================================#
  # Continuous variables
  #=============================================================================#
  
  if (is.numeric(Y_train[,f])){
    
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "MAE"
    
    # Risk factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[,!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = regressionSummary,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = riskfactor,
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = FALSE)
    
    optPar <- which.min(fit$results$MAE)
    optAlpha <- fit$results$alpha[optPar]
    optLambda <- fit$results$lambda[optPar]
    PerformanceList[[f]] <- list(optAlpha, optLambda)
    
  } 
  
  #=============================================================================#
  # Discrete variables
  #=============================================================================#
  if (!is.numeric(Y_train[,f])){
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[,!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = factor(riskfactor, levels = c("No", "Yes")),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    optPar <- which.max(fit$results$ROC)
    optAlpha <- fit$results$alpha[optPar]
    optLambda <- fit$results$lambda[optPar]
    PerformanceList[[f]] <- list(optAlpha, optLambda)
    
  }
  
  trainList[[f]] <- fit$results
  predList[[f]] <- fit$pred
  save(trainList,predList, file = "EXTEND/OutputList_None_EN.RData")
  save(PerformanceList, file = "EXTEND/PerformanceList_None_EN.RData")
  
}

names(trainList) <- colnames(Y_train)
names(predList) <- colnames(Y_train)

# Save best models
save(trainList,predList, file = "EXTEND/OutputList_None_EN.RData")

names(PerformanceList) <- names(trainList)
save(PerformanceList, file = "EXTEND/PerformanceList_None_EN.RData")


# Train model with optimal hyperparameters
bestModel <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for lambda
  lambdaCV <- PerformanceList[[f]][[2]]
  
  # Set grid for alpha
  alphaCV <- PerformanceList[[f]][[1]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  colnames(parameterGrid) <- c(".alpha", ".lambda")
  
  # ML method
  MLmethod = "glmnet"
  
  if (is.numeric(Y_train[,f])){
    # Performance metric
    performance_metric = "MAE"
    
    # factor
    factor_cor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[,!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none")
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor_cor,
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = FALSE)
    
    
  }
  if (!is.numeric(Y_train[,f])){
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    factor_cor <- ifelse(riskfactor == "Yes",1,0)
    
    # dataMatrix
    dataMatrix <- X_train[,!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none", 
                               classProbs = TRUE)
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor(riskfactor, levels = c("No", "Yes")),
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = TRUE)
    
    
  }
  
}
names(bestModel) <- colnames(Y_train)

# Save best models
bestModel_34 <- list(bestModel[[3]], bestModel[[4]])
names(bestModel_34) <- colnames(Y_train)[3:4]
save(bestModel_34, file = "EXTEND/bestModel_None_EN_34.RData")



# Pred vs Obs in CV
perf_discrete <- rep(NA,7)
for (i in 3:9){
  
  predResults <- NULL
  
  test <- predList[[i]]
  test <- test[(test$alpha == PerformanceList[[i]][[1]]) &
                 (test$lambda == PerformanceList[[i]][[2]]),c("obs", "Yes")]
  predResults <- rbind.data.frame(predResults,test)
  
  perf_discrete[i-2] <- as.numeric(pROC::auc(pROC::roc(predResults$obs, predResults$Yes)))
}

AUCdf <- data.frame(RiskFactor = colnames(Y_train)[3:9],
                    AUC = perf_discrete)
save(AUCdf, file = "EXTEND/AUCdf_None_EN.RData")

# Pred vs Obs in CV
perf_continuous <- rep(NA,2)
for (i in 1:2){
  
  predResults <- NULL
  
  test <- predList[[i]]
  test <- test[(test$alpha == PerformanceList[[i]][[1]]) &
                 (test$lambda == PerformanceList[[i]][[2]]),c("obs", "pred")]
  predResults <- rbind.data.frame(predResults,test)
  
  perf_continuous[i] <- as.numeric(R2(pred = predResults$pred, obs = predResults$obs))
}

R2df <- data.frame(RiskFactor = colnames(Y_train)[1:2],
                   R2 = perf_continuous)
save(R2df, file = "EXTEND/R2df_None_EN.RData")


#*****************************************************************************#
# Literature ElasticNet
#*****************************************************************************#
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("EXTEND/selectedProbes_lit.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

selectedProbes <- list(SysBP = selectedProbes_lit$BP,
                       TotalChol = selectedProbes_lit$TC,
                       Education = selectedProbes_lit$EA,
                       Physical = selectedProbes_lit$Physical,
                       Diet = selectedProbes_lit$Diet,
                       Depression = selectedProbes_lit$MDD,
                       Diabetes = selectedProbes_lit$T2D,
                       HeartDisease = selectedProbes_lit$CHD)

X_train <- methSet_allNorm_fil

# Convert to M-values
X_train <- log2(X_train/(1-X_train))

# Test if samples are in correct order
all(colnames(X_train) == rownames(Y_train))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Machine learning method
MLmethod = "glmnet"

# Model training
trainList <- list()
predList <- list()
PerformanceList <- list()
for (f in 1:(ncol(Y_train)-1)){
  
  
  #=============================================================================#
  # Continuous variables
  #=============================================================================#
  
  if (is.numeric(Y_train[,f])){
    
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "MAE"
    
    # Risk factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = regressionSummary,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = riskfactor,
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = FALSE)
    
    optPar <- which.min(fit$results$MAE)
    optAlpha <- fit$results$alpha[optPar]
    optLambda <- fit$results$lambda[optPar]
    PerformanceList[[f]] <- list(optAlpha, optLambda)
    
  } 
  
  #=============================================================================#
  # Discrete variables
  #=============================================================================#
  if (!is.numeric(Y_train[,f])){
    # Set grid for lambda
    lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
    
    # Set grid for alpha
    alphaCV <- seq(0.1,1,length.out = 10)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = factor(riskfactor, levels = c("No", "Yes")),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    optPar <- which.max(fit$results$ROC)
    optAlpha <- fit$results$alpha[optPar]
    optLambda <- fit$results$lambda[optPar]
    PerformanceList[[f]] <- list(optAlpha, optLambda)
    
  }
  
  trainList[[f]] <- fit$results
  predList[[f]] <- fit$pred
  save(trainList,predList, file = "EXTEND/OutputList_Lit_EN.RData")
  save(PerformanceList, file = "EXTEND/PerformanceList_Lit_EN.RData")
  
}

names(trainList) <- colnames(Y_train)[-9]
names(predList) <- colnames(Y_train)[-9]

# Save best models
save(trainList,predList, file = "EXTEND/OutputList_Lit_EN.RData")

names(PerformanceList) <- names(trainList)[-9]
save(PerformanceList, file = "EXTEND/PerformanceList_Lit_EN.RData")


# Train model with optimal hyperparameters
bestModel <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for lambda
  lambdaCV <- PerformanceList[[f]][[2]]
  
  # Set grid for alpha
  alphaCV <- PerformanceList[[f]][[1]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  colnames(parameterGrid) <- c(".alpha", ".lambda")
  
  # ML method
  MLmethod = "glmnet"
  
  if (is.numeric(Y_train[,f])){
    # Performance metric
    performance_metric = "MAE"
    
    # factor
    factor_cor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none")
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor_cor,
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = FALSE)
    
    
  }
  if (!is.numeric(Y_train[,f])){
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    factor_cor <- ifelse(riskfactor == "Yes",1,0)
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none", 
                               classProbs = TRUE)
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor(riskfactor, levels = c("No", "Yes")),
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = TRUE)
    
    
  }
  
}
names(bestModel) <- colnames(Y_train)[-9]

# Save best models
save(bestModel, file = "EXTEND/bestModel_Lit_EN.RData")



# Pred vs Obs in CV
perf_discrete <- rep(NA,6)
for (i in 3:8){
  
  predResults <- NULL
  
  test <- predList[[i]]
  test <- test[(test$alpha == PerformanceList[[i]][[1]]) &
                 (test$lambda == PerformanceList[[i]][[2]]),c("obs", "Yes")]
  predResults <- rbind.data.frame(predResults,test)
  
  perf_discrete[i-2] <- as.numeric(pROC::auc(pROC::roc(predResults$obs, predResults$Yes)))
}

AUCdf <- data.frame(RiskFactor = colnames(Y_train)[3:8],
                    AUC = perf_discrete)
save(AUCdf, file = "EXTEND/AUCdf_Lit_EN.RData")

# Pred vs Obs in CV
perf_continuous <- rep(NA,2)
for (i in 1:2){
  
  predResults <- NULL
  
  test <- predList[[i]]
  test <- test[(test$alpha == PerformanceList[[i]][[1]]) &
                 (test$lambda == PerformanceList[[i]][[2]]),c("obs", "pred")]
  predResults <- rbind.data.frame(predResults,test)
  
  perf_continuous[i] <- as.numeric(R2(pred = predResults$pred, obs = predResults$obs))
}

R2df <- data.frame(RiskFactor = colnames(Y_train)[1:2],
                   R2 = perf_continuous)
save(R2df, file = "EXTEND/R2df_Lit_EN.RData")


#*****************************************************************************#
# Literature Random Forest
#*****************************************************************************#
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("EXTEND/Y_train_factors.RData")
load("EXTEND/selectedProbes_lit.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

selectedProbes <- list(SysBP = selectedProbes_lit$BP,
                       TotalChol = selectedProbes_lit$TC,
                       Education = selectedProbes_lit$EA,
                       Physical = selectedProbes_lit$Physical,
                       Diet = selectedProbes_lit$Diet,
                       Depression = selectedProbes_lit$MDD,
                       Diabetes = selectedProbes_lit$T2D,
                       HeartDisease = selectedProbes_lit$CHD)

X_train <- methSet_allNorm_fil

# Convert to M-values
X_train <- log2(X_train/(1-X_train))

# Test if samples are in correct order
all(colnames(X_train) == rownames(Y_train))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Machine learning method
MLmethod = "ranger"

# Model training
trainList <- list()
predList <- list()
PerformanceList <- list()
for (f in 1:(ncol(Y_train)-1)){
  
  
  #=============================================================================#
  # Continuous variables
  #=============================================================================#
  
  if (is.numeric(Y_train[,f])){
    
    # Number of randomly selected predictors
    mtry_CV <- c(1000, 2000,3000,4000,5000,6000,7000,8000)
    
    # split rule
    splitrule_CV <- "variance"
    
    # minimal node size
    min.node.size_CV = c(10,20,30,40,50,60)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
    colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
    
    # Performance metric
    performance_metric = "MAE"
    
    # Risk factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = regressionSummary,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = riskfactor,
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = FALSE)
    
    optPar <- which.min(fit$results$MAE)
    opt_mtry <- trainResults[[1]]$mtry[optPar]
    opt_splitrule <- trainResults[[1]]$splitrule[optPar]
    opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
    
    PerformanceList[[i]] <- list(opt_mtry, opt_splitrule, opt_min.node.size)
    
  } 
  
  #=============================================================================#
  # Discrete variables
  #=============================================================================#
  if (!is.numeric(Y_train[,f])){
    # Number of randomly selected predictors
    mtry_CV <- c(100,500,1000,1500,2000,3000,4000)
    
    # split rule
    splitrule_CV <- "gini"
    
    # minimal node size
    min.node.size_CV = c(3,5,10,15,20)
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
    colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
    
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Create indices for cross-validation
    set.seed(123)
    CVindex <- NULL
    for (r in 1:5){
      temp <- createFolds(1:length(riskfactor),5, returnTrain = TRUE)
      names(temp) <- paste0(names(temp),".Rep",r)
      CVindex <- c(CVindex, temp)
    }
    rm(temp)
    
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = CVindex)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = factor(riskfactor, levels = c("No", "Yes")),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    optPar <- which.max(fit$results$ROC)
    opt_mtry <- trainResults[[1]]$mtry[optPar]
    opt_splitrule <- trainResults[[1]]$splitrule[optPar]
    opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
    
    PerformanceList[[i]] <- list(opt_mtry, opt_splitrule, opt_min.node.size)
    
  }
  
  trainList[[f]] <- fit$results
  predList[[f]] <- fit$pred
  save(trainList,predList, file = "EXTEND/OutputList_Lit_RF.RData")
  save(PerformanceList, file = "EXTEND/PerformanceList_Lit_RF.RData")
  
}

names(trainList) <- colnames(Y_train)[-9]
names(predList) <- colnames(Y_train)[-9]

# Save best models
save(trainList,predList, file = "EXTEND/OutputList_Lit_RF.RData")

names(PerformanceList) <- names(trainList)[-9]
save(PerformanceList, file = "EXTEND/PerformanceList_Lit_RF.RData")


# Train model with optimal hyperparameters
bestModel <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for mtry
  mtry_CV <- PerformanceList[[f]][[1]]
  
  # Set grid for splitrule
  splitrule_CV <- PerformanceList[[f]][[2]]
  
  # Set grid for splitrule
  min.node.size_CV <- PerformanceList[[f]][[3]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # ML method
  MLmethod = "ranger"
  
  if (is.numeric(Y_train[,f])){
    # Performance metric
    performance_metric = "MAE"
    
    # factor
    factor_cor <- Y_train[!is.na(Y_train[,f]),f]
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none")
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor_cor,
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = FALSE)
    
    
  }
  if (!is.numeric(Y_train[,f])){
    
    # Performance metric
    performance_metric = "ROC"
    
    # factor
    riskfactor <- Y_train[!is.na(Y_train[,f]),f]
    factor_cor <- ifelse(riskfactor == "Yes",1,0)
    
    # dataMatrix
    dataMatrix <- X_train[rownames(X_train) %in% selectedProbes[[f]],!is.na(Y_train[,f])]
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "none", 
                               classProbs = TRUE)
    
    # Actual training
    set.seed(123)
    bestModel[[f]] <- train(x = t(dataMatrix),
                            y = factor(riskfactor, levels = c("No", "Yes")),
                            metric= performance_metric,
                            method = MLmethod,
                            tuneGrid = parameterGrid,
                            trControl = fitControl,
                            maximize = TRUE)
    
    
  }
  
}
names(bestModel) <- colnames(Y_train)

# Save best models
save(bestModel, file = "EXTEND/bestModel_Lit_RF.RData")

# Pred vs Obs in CV
perf_discrete <- rep(NA,6)
for (i in 3:8){
  
  predResults <- NULL
  test <- predList[[i]]
  test <- test[(test$mtry == PerformanceList[[i]][[1]]) &
                 (test$splitrule == PerformanceList[[i]][[2]]) &
                 (test$min.node.size == PerformanceList[[i]][[3]]),c("obs", "Yes")]
  predResults <- rbind.data.frame(predResults,test)
  perf_discrete[i-2] <- as.numeric(pROC::auc(pROC::roc(predResults$obs, predResults$Yes)))
}

AUCdf <- data.frame(RiskFactor = colnames(Y_train)[3:8],
                    AUC = perf_discrete)
save(AUCdf, file = "EXTEND/AUCdf_Lit_RF.RData")

# Pred vs Obs in CV
perf_continuous <- rep(NA,2)
for (i in 1:2){
  
  predResults <- NULL
  test <- predList[[i]]
  test <- test[(test$mtry == PerformanceList[[i]][[1]]) &
                 (test$splitrule == PerformanceList[[i]][[2]]) &
                 (test$min.node.size == PerformanceList[[i]][[3]]),c("obs", "pred")]
  predResults <- rbind.data.frame(predResults,test)
  
  perf_continuous[i] <- as.numeric(R2(pred = predResults$pred, obs = predResults$obs))
}

R2df <- data.frame(RiskFactor = colnames(Y_train)[1:2],
                   R2 = perf_continuous)
save(R2df, file = "EXTEND/R2df_Lit_RF.RData")



################################################################################

# Get best models

################################################################################
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)
finalModels <- list()

load("EXTEND/Y_train_factors.RData")

# Systolic blood pressure
load("EXTEND/bestModel_None_EN_12.RData") # Change 12
finalModels[[1]] <- bestModel_12[["SysBP"]]

# Total Cholesterol
load("EXTEND/bestModel_Cor_RF.RData")
finalModels[[2]] <- bestModel[["TotalChol"]]

# Other factors
load("EXTEND/bestModel_Cor_EN.RData")
for (i in 3:9){
  finalModels[[i]] <- bestModel[[i]]
}

all(names(bestModel) == colnames(Y_train))
names(finalModels) <- colnames(Y_train)
save(finalModels, file = "EXTEND/finalModels.RData")



