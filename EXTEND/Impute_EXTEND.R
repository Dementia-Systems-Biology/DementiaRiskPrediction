# File: Impute_EXTEND.R
# Author: Jarno Koetsier
# Date: August 6, 2023

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EXTEND/Data/lumi_dpval_EXTEND.RData")
load("EXTEND/Data/methSet_allNorm_EXTEND.RData")
load("EXTEND/Data/metaData_ageFil.RData")

# Convert to M-values
X_EXTEND_M <- log2(methSet_allNorm/(1-methSet_allNorm))

# Set low quality probes to NA
X_EXTEND_mis <- X_EXTEND_M
X_EXTEND_mis[lumi_dpval > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_EXTEND_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_EXTEND_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_EXTEND_copy <- X_EXTEND_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_EXTEND_mis[removeProbes[j,1], removeProbes[j,2]]
  X_EXTEND_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}

# Optimize for number of PCs
pred <- matrix(NA, nrow = 100,ncol = 20)
nPCs_CV <- 10
for (p in 3:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_EXTEND_imp_CV <- imputePCA(t(X_EXTEND_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_EXTEND_imp_CV[removeProbes[j,2], removeProbes[j,1]]
  }
  save(pred, file = "EXTEND/pred_imp.RData")
}

load("EXTEND/pred_imp.RData")

# Get number of PCs with minimal MAE
test <- apply(pred,2,function(x) MAE(obs = values,pred = x))
plot(test)
optPC <- which.min(test)

# Perform imputation
X_EXTEND_imp <- imputePCA(t(X_EXTEND_mis),ncp = optPC)$completeObs

# Save data
save(X_EXTEND_imp, file = "EXTEND/Data/X_EXTEND_imp.RData")
