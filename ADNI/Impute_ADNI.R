# ============================================================================ #
# File: Impute_ADNI.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Impute low quality methylation values in the ADNI cohort.
# ============================================================================ #

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
load("ADNI/Data/methSet_allNorm_ADNI.RData")
load("ADNI/Data/lumi_dpval_ADNI.RData")
load("ADNI/Data/MetaData_ADNI.RData")

# Prepare data
X_ADNI <- methSet_allNorm
all(colnames(X_ADNI) == colnames(lumi_dpval))

# optimize imputation
X_ADNI_M <- log2(X_ADNI/(1-X_ADNI))
test <- lumi_dpval[rownames(X_ADNI_M), colnames(X_ADNI_M)]
X_ADNI_mis <- X_ADNI_M
X_ADNI_mis[test > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_ADNI_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_ADNI_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_ADNI_copy <- X_ADNI_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_ADNI_mis[removeProbes[j,1], removeProbes[j,2]]
  X_ADNI_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}


# Perform imputation for different number of principal components (PCs)
nPCs_CV <- 15
pred <- matrix(NA, nrow = 100,ncol = nPCs_CV)
for (p in 3:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_ADNI_imp_CV <- imputePCA(t(X_ADNI_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_ADNI_imp_CV[removeProbes[j,2], removeProbes[j,1]]
    save(pred, file = "ADNI/pred_imp.RData")
  }
}
save(pred, file = "ADNI/Data/pred_imp.RData")

# Find optimal number of PCs
test <- apply(pred,2,function(x) MAE(obs = values,pred = x))
optPC <- which.min(test)

# Impute the data
set.seed(123)
X_ADNI_imp <- imputePCA(t(X_ADNI_mis),ncp = optPC)$completeObs

# Save the imputed data
save(X_ADNI_imp, file = "ADNI/Data/X_ADNI_imp.RData")






















