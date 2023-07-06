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

#*****************************************************************************#
# Adherence to a Mediterranean diet (-1.7)
#*****************************************************************************#

# Reported amount of fruits and vegetables consumed by the participant the previous day. 
# A healthy diet was defined as consuming five or more portions of fruits and vegetables on a daily basis [8].
# [8]	National Health Service (2009) 5 A Day. NHS, London. (CHECK FOR UPDATE)
# Variables in EXTEND: 
# dat$Fruit     ==  Fruit(portions per day)	        0    1-2     3-4     5-6     >6
# dat$Vegetables ==  Vegetables(portions per day)	  0    1-2     3-4     5-6     >6

dat <- dat[-which(dat$Fruit=="."|dat$Vegtables=="."),]
dat$Fruit <- as.numeric(dat$Fruit)
dat$Vegtables <- as.numeric(dat$Vegtables)

dat$MEDITERANIAN <- ifelse(dat$Fruit > 3 | dat$Vegtables > 3, 1,0)
table(dat$MEDITERANIAN)

# Add weight to column
dat$MEDITERANIAN<-ifelse(dat$Fruit > 3 | dat$Vegtables > 3, -1.7,0)


#*****************************************************************************#
# Physical inactivity (+1.1)
#*****************************************************************************#

# Variables in EXTEND: Exercise causing  increased pulse more than 2.5 hours per week
dat <- dat[-which(dat$Exercise.increased.pulse.more.than.2halfhrsawk=="."),]
dat$PHYSICAL_INACTIVITY<-ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 0, 1, 0)
table(dat$PHYSICAL_INACTIVITY)

# Add weight to column
dat$PHYSICAL_INACTIVITY<-ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 0, 1.1,0)


#*****************************************************************************#
# Smoking (+1.5)
#*****************************************************************************#

# Variables in EXTEND: Do you currently smoke
dat <- dat[!is.na(dat$Do.you.currently.smoke),]
dat$SMOKING <- ifelse(dat$Do.you.currently.smoke==1,1,0)
table(dat$SMOKING)

# Add weight to column
dat$SMOKING<-ifelse(dat$Do.you.currently.smoke==1,1.5,0)


#*****************************************************************************#
# Low-to-moderate alcohol intake (-1)
#*****************************************************************************#

# UK: Low-to-moderate alcohol use was defined as 1-14 glasses per week according to recent UK alcohol guidelines [7].
#[7]	Department of Health (2016) Alcohol guidelines review - report from the guidelines development group to the UK Chief Medical Officers. Department of Health.
# Variables in EXTEND:
table(dat$alcoholic.drinks.per.day)
table(dat$How.often.do.you.drink.alcohol)
dat <- dat[-which(dat$How.often.do.you.drink.alcohol=="."),]
x <- dat[which(dat$alcoholic.drinks.per.day=="."),c("alcoholic.drinks.per.day","How.often.do.you.drink.alcohol")]

#alcoholic.drinks.per.day       ==1-2(1)     3-4(2)     5-6(3)    7-9(4)     >10(5)
#How.often.do.you.drink.alcohol ==Never(0)    1 a month(1)     2-4 mth(2)    2-3 wk(3)    =>4 wk(4)

# Low-moderate alcohol intake
dat$LtoMAlcohol <- NULL
dat[dat$How.often.do.you.drink.alcohol==0 ,"LtoMAlcohol"] <- 1
dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==3 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==4 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==5 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1

# High alcohol intake
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3  &  dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==2 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day=="." &    dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0

# Remove NA rows
dat<-dat[!is.na(dat$LtoMAlcohol),]

# Check table
table(dat$LtoMAlcohol)

# Add weight to column
dat$LtoMAlcohol<-ifelse(dat$LtoMAlcohol == 1,-1,0)


#*****************************************************************************#
# Obesity (+1.6)
#*****************************************************************************#

dat$OBESITY <- ifelse(dat$BMI >=30, 1, 0)
table(dat$OBESITY)

# Add weight to column
dat$OBESITY<-ifelse(dat$BMI >=30,1.6,0)


#*****************************************************************************#
# Depression (+2.1)
#*****************************************************************************#

table(dat$Depression)
dat$DEPRESSION<-dat$Depression

# Add weight to column
dat$DEPRESSION <- ifelse(dat$DEPRESSION==1,2.1,0)


#*****************************************************************************#
# Type-2-Diabetes (+1.3)
#*****************************************************************************#

table(dat$T2.Diabetes)
dat$DIABETEII <- dat$T2.Diabetes

# Add weight to column
dat$DIABETEII<-ifelse(dat$DIABETEII==1,1.3,0)

#*****************************************************************************#
# Hypertension (+1.6)
#*****************************************************************************#
# Mean systolic blood pressure >= 140 mm Hg or mean diastolic blood pressure >= 90 mm Hg [3].

# No missing values
table(is.na(dat$MeanSysBP))
table(is.na(dat$MeanDiaBP))

dat$HYPERTENSTION <- ifelse (dat$MeanSysBP >=140 | dat$MeanDiaBP >= 90,1,0)
table(dat$HYPERTENSTION)

# Add weight to column
dat$HYPERTENSTION<-ifelse (dat$MeanSysBP >= 140 | dat$MeanDiaBP >= 90,1.6,0)

#*****************************************************************************#
# High cholesterol(+1.4)
#*****************************************************************************#

# Total cholesterol level of >= 5.0 mmol/L and low-density lipoprotein of >= 3.0 mmol/L, 
# following the guidelines of the National Health Service UK [2].

# Remove missing values
dat <- dat[!is.na(dat$HDL_unloged),]
dat$HDL_unloged <- as.numeric(dat$HDL_unloged)


#https://academic.oup.com/view-large/figure/342663278/cvab164f2.tif
dat$Highcholesterol<-ifelse(dat$HDL_unloged > 2.2, 1,0)
table(dat$Highcholesterol)

# Add weight to column
dat$Highcholesterol<-ifelse(dat$HDL_unloged > 2.2, 1.4,0)


#*****************************************************************************#
# Heart disease (+1.0)
#*****************************************************************************#

dat$Heartdisease <- dat$Heart.Disease

# Add weight to column
dat$Heartdisease<-ifelse(dat$Heartdisease==1,1,0)
table(dat$Heartdisease)

#*****************************************************************************#
# Chronic kidney disease (+1.1)
#*****************************************************************************#

# No missing values
table(is.na(dat$Kidney..Disease))
dat$kidneydisease<-dat$Kidney..Disease

# Add weight to column
dat$kidneydisease<-ifelse(dat$kidneydisease==1,1.1,0)
table(dat$kidneydisease)


#*****************************************************************************#
# Calculate scores
#*****************************************************************************#
Y_train <- dat[,c(1:8,11,12,151:161)]
Y_train$X <- NULL
Y_train$LIBRA <- rowSums(Y_train[,10:20])

# Save data
Y_train <- unique(Y_train)
save(Y_train, file = "EXTEND/Y_train_LIBRA.RData")

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
load("EXTEND/Y_train_LIBRA.RData")
load("Data/methSet_allNorm_fil.RData")

X_train <- methSet_allNorm_fil[,Y_train$Basename]
all(colnames(X_train) == Y_train$Basename)


#*****************************************************************************#
# correlation-based selection
#*****************************************************************************#

Y_cor <- Y_train[,10:20]
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

save(correlations, file = "EXTEND/correlations_LIBRA.RData")


# Selected probes
selectedProbes <- list()
for (i in 1:ncol(correlations)){
  selectedProbes[[i]] <- names(tail(sort(abs(correlations[,i])),10000))
}
names(selectedProbes) <- colnames(correlations)
save(selectedProbes, file = "EXTEND/selectedProbes_LIBRA.RData")


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
load("EXTEND/Y_train_LIBRA.RData")
load("EXTEND/selectedProbes_LIBRA.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil[,Y_train$Basename]
all(colnames(X_train) == Y_train$Basename)

# Convert to M-values
X_train <- log2(X_train/(1-X_train))

# Test if samples are in correct order
rownames(Y_train) <- Y_train$Basename
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
features <- unique(unlist(selectedProbes))
dataMatrix <- X_train[features,]

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
  factors <- Y_train[index[[1]],10:20]
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
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1100))
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
               y = Y_train$LIBRA,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  
  trainResults[[i]] <- fit$results
  predResults[[i]] <- fit$pred
}

# Save results
save(trainResults, predResults, file = "EXTEND/LIBRA_OutputList_Cor_EN.RData")


#load("EXTEND/CAIDE1_OutputList_Cor_EN.RData")

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
load("EXTEND/correlations_LIBRA.RData")
probes <- list()
for (p in 1:ncol(correlations)){
  probes[[p]] <- names(tail(sort(abs(correlations[,p])),1100))
}


# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < ncol(correlations)){
    n = n + 1
  } else {
    n = 1
  }
}

# Settings for repeated cross-validation
fitControl <- trainControl(method = "none")

# Actual training
set.seed(123)
LIBRA_Model <- train(x = t(X_train[finalProbes,]),
                      y = Y_train$LIBRA,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)


# Save model
save(LIBRA_Model, file = "EXTEND/LIBRA_Model_Cor_EN.RData")


ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$alpha == optAlpha) & 
                             (predResults[[i]]$lambda == optLambda),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

R2(ObsvsPred$pred, ObsvsPred$obs)

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
load("EXTEND/Y_train_LIBRA.RData")
load("EXTEND/selectedProbes_CAIDE1.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil[,Y_train$Basename]
all(colnames(X_train) == Y_train$Basename)

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
features <- unique(unlist(selectedProbes))
dataMatrix <- X_train[features,]

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
  factors <- Y_train[index[[1]],10:20]
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
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1100))
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
               y = Y_train$LIBRA,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  
  trainResults[[i]] <- fit$results
  predResults[[i]] <- fit$pred
}

# Save results
save(trainResults, predResults, file = "EXTEND/LIBRA_OutputList_Cor_RF.RData")


#load("EXTEND/CAIDE1_OutputList_Cor_RF.RData")
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
load("EXTEND/correlations_LIBRA.RData")
probes <- list()
for (p in 1:ncol(correlations)){
  probes[[p]] <- names(tail(sort(abs(correlations[,p])),1100))
}


# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < ncol(correlations)){
    n = n + 1
  } else {
    n = 1
  }
}


# Settings for repeated cross-validation
fitControl <- trainControl(method = "none")

# Actual training
set.seed(123)
LIBRA_Model <- train(x = t(X_train[finalProbes,]),
                      y = Y_train$LIBRA,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)



# Save model
save(LIBRA_Model, file = "EXTEND/LIBRA_Model_Cor_RF.RData")


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
load("EXTEND/Y_train_LIBRA.RData")
load("Data/methSet_allNorm_fil.RData")
source("FUN_MachineLearning.R")

X_train <- methSet_allNorm_fil[,Y_train$Basename]

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
  temp <- createFolds(1:nrow(Y_train),5, returnTrain = TRUE)
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
LIBRA_Model <- train(x = t(dataMatrix),
                      y = Y_train$LIBRA,
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = FALSE)

# Save model
save(LIBRA_Model, file = "EXTEND/LIBRA_Model_Cor_None.RData")

