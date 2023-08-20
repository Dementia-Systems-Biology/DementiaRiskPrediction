# ============================================================================ #
# File: PredictCAIDE1_EXTEND.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Predict CAIDE1 score in the EXTEND cohort.
# ============================================================================ #

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
load("EXTEND/Data/metaData_ageFil.RData")
load("EXTEND/Data/methSet_allNorm_fil.RData")

#################################################################################

# Prepare phenotype data (Y)

################################################################################

#*****************************************************************************#
# Age
#*****************************************************************************#
# Discretize age in 3 classes with score 0, 3, and 4.
dat$age_c <- rep(NA,nrow(dat))
dat$age_c[dat$Age < 47] <- 0
dat$age_c[dat$Age >= 47 & dat$Age <= 53] <- 3
dat$age_c[dat$Age > 53] <- 4
table(dat$age_c)

#*****************************************************************************#
# Sex
#*****************************************************************************#
# Male = 1 (1 in dataset)
# Female = 0 (2 in dataset)

dat$Sex_c <- ifelse(dat$Sex == 1,1,0) # Check!!!
table(dat$Sex_c)


#*****************************************************************************#
# Education
#*****************************************************************************#
# "College.or.Uni.degree"                          "A.level.AS.level.or.equiv"                                    
# "O.level.GCSEs.or.equiv"                         "CSEs.or.equiv"                                                
# "NVQ.HND.HNC.or.equiv"                            "Other.professional.quals"                                     
# "None.of.the.above
dat$Edu_c <- ifelse(dat$None.of.the.above == 1,2,0)
table(dat$Edu_c)


#*****************************************************************************#
# Systolic blood pressure
#*****************************************************************************#
dat$Syst_c <- ifelse(dat$MeanSysBP <= 140,0,2)
table(dat$Syst_c)

#*****************************************************************************#
# BMI
#*****************************************************************************#
dat$BMI_c<-ifelse(dat$BMI <= 30,0,2)
table(dat$BMI_c)

#*****************************************************************************#
# Serum total cholesterol level
#*****************************************************************************#
dat <- dat[-which(dat$Chol_unloged=="."),]
dat$Chol_unloged <- as.numeric(dat$Chol_unloged)
dat$Chol_c <- ifelse(dat$Chol_unloged <= 6.5,0,2)
table(dat$Chol_c)

#*****************************************************************************#
# Physical activity
#*****************************************************************************#
dat$PHYSICAL_c <- ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 1,0,1)
table(dat$PHYSICAL_c)


#*****************************************************************************#
# Calculate scores
#*****************************************************************************#
Y_train <- dat[,c(1:8,11,12,151:157)]
Y_train$X <- NULL

Y_train$CAIDE1 <- rowSums(Y_train[,10:16])

# Save data
Y_train <- unique(Y_train)
save(Y_train, file = "EXTEND/Data/Y_train_CAIDE1.RData")

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
load("EXTEND/Data/Y_train_CAIDE1.RData")
load("EXTEND/Data/methSet_allNorm_fil.RData")

X_train <- methSet_allNorm_fil
all(colnames(X_train) == rownames(Y_train))


#*****************************************************************************#
# correlation-based selection
#*****************************************************************************#

# Calculate correlations
Y_cor <- Y_train[,10:16]
correlations <- matrix(NA, nrow = nrow(X_train), ncol = ncol(Y_cor))
for (i in 1:ncol(Y_cor)) {
  factor <- Y_cor[,i]
  correlations[,i] <- apply(X_train, 1, 
                            function(x){cor(x, 
                                            factor, 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_train)
colnames(correlations) <- colnames(Y_cor)

# Save correlations
save(correlations, file = "Models/CAIDE_LIBRA/correlations_LIBRA.RData")

# Selected probes
selectedProbes <- list()
for (i in 1:ncol(correlations)){
  selectedProbes[[i]] <- names(tail(sort(abs(correlations[,i])),10000))
}
names(selectedProbes) <- colnames(correlations)

# Save selected probes
save(selectedProbes, file = "Models/CAIDE_LIBRA/selectedProbes_CAIDE1.RData")


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
load("EXTEND/Data/Y_train_CAIDE1.RData")
load("Models/CAIDE_LIBRA/selectedProbes_CAIDE1.RData")
load("EXTEND/Data/methSet_allNorm_fil.RData")
source("EXTEND/FUN_MachineLearning.R")

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

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Performance metric
performance_metric = "MAE"

# dataMatrix
dataMatrix <- X_train[selectedProbes,]

# Create indices for cross-validation
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(Y_train),5, returnTrain = TRUE)
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
  
  # Calculate correlations with CAIDE1 factors
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < ncol(correlations_CV)){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
                             summaryFunction = regressionSummary,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(dataMatrix[finalProbes,]),
               y = Y_train$CAIDE1,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  
  trainResults[[i]] <- fit$results
  predResults[[i]] <- fit$pred
}

# Save results
save(trainResults, predResults, file = "Models/CAIDE_LIBRA/CAIDE1_OutputList_Cor_EN.RData")

# Get optimal hyperparameters
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]


# Train model with optimal hyperparameters

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# ML method
MLmethod = "glmnet"

# Performance metric
performance_metric = "MAE"

# Select top correlated features for each factor
probes <- list()
for (p in 1:ncol(correlations)){
  probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
}

# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < ncol(correlations_CV)){
    n = n + 1
  } else {
    n = 1
  }
}

# Settings for repeated cross-validation
fitControl <- trainControl(method = "none")

# Actual training
set.seed(123)
CAIDE1_Model <- train(x = t(dataMatrix[finalProbes,]),
                      y = Y_train$CAIDE1,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)


# Save model
save(CAIDE1_Model, file = "Models/CAIDE_LIBRA/CAIDE1_Model_Cor_EN.RData")




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
load("EXTEND/Data/Y_train_CAIDE1.RData")
load("Models/CAIDE_LIBRA/selectedProbes_CAIDE1.RData")
load("EXTEND/Data/methSet_allNorm_fil.RData")
source("EXTEND/FUN_MachineLearning.R")

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

# dataMatrix
dataMatrix <- X_train[selectedProbes,]

# Create indices for cross-validation
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(Y_train),5, returnTrain = TRUE)
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
  
  # Calculate correlations with CAIDE1 factors
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < ncol(correlations_CV)){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
                             summaryFunction = regressionSummary,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(dataMatrix[finalProbes,]),
               y = Y_train$CAIDE1,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  
  trainResults[[i]] <- fit$results
  predResults[[i]] <- fit$pred
}

# Save results
save(trainResults, predResults, file = "Models/CAIDE_LIBRA/CAIDE1_OutputList_Cor_RF.RData")

# Get optimal hyperparameters
perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

# Train model with optimal hyperparameters

# Set grid for mtry
mtry_CV <- opt_mtry

# Set grid for splitrule
splitrule_CV <- opt_splitrule

# Set grid for splitrule
min.node.size_CV <- opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# ML method
MLmethod = "ranger"

# Performance metric
performance_metric = "MAE"

# Select top correlated features for each factor
probes <- list()
for (p in 1:ncol(correlations)){
  probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
}

# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < ncol(correlations_CV)){
    n = n + 1
  } else {
    n = 1
  }
}

# Settings for repeated cross-validation
fitControl <- trainControl(method = "none")

# Actual training
set.seed(123)
CAIDE1_Model <- train(x = t(dataMatrix[finalProbes,]),
                      y = Y_train$CAIDE1,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)



# Save model
save(CAIDE1_Model, file = "Models/CAIDE_LIBRA/CAIDE1_Model_Cor_RF.RData")


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
load("EXTEND/Data/Y_train_factors.RData")
load("EXTEND/Data/methSet_allNorm_fil.RData")
source("EXTEND/FUN_MachineLearning.R")

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


# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Performance metric
performance_metric = "MAE"

# dataMatrix
dataMatrix <- X_train

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
CAIDE1_Model <- train(x = t(dataMatrix),
                      y = Y_train$CAIDE1,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)

# Save model
save(CAIDE1_Model, file = "Models/CAIDE_LIBRA/CAIDE1_Model_None_EN.RData")

