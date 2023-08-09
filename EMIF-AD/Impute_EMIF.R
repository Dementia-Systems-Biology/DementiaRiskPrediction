# ============================================================================ #
# File: Impute_EMIF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Impute low quality methylation values in the EMIF-AD cohort.
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
load("EMIF-AD/Data/lumi_dpval_EMIF.RData")
load("EMIF-AD/Data/X_EMIF.RData")
load("EMIF-AD/Data/metaData_EMIF.RData")

# optimize imputation
X_EMIF_M <- log2(X_EMIF/(1-X_EMIF))
test <- lumi_dpval[rownames(X_EMIF_M), colnames(X_EMIF_M)]
X_EMIF_mis <- X_EMIF_M
X_EMIF_mis[test > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_EMIF_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_EMIF_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_EMIF_copy <- X_EMIF_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_EMIF_mis[removeProbes[j,1], removeProbes[j,2]]
  X_EMIF_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}

# Perform imputation for different number of principal components (PCs)
nPCs_CV <- 15
pred <- matrix(NA, nrow = 100,ncol = nPCs_CV)
for (p in 1:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_EMIF_imp_CV <- imputePCA(t(X_EMIF_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_EMIF_imp_CV[removeProbes[j,2], removeProbes[j,1]]
  }
  save(pred, file = "EMIF-AD/Data/predImp.RData")
}

# Find optimal number of PCs
load("EMIF-AD/Data/predImp.RData")
test <- apply(pred,2,function(x) MAE(obs = values,pred = x))
optPC <- which.min(test)

# Impute the data
X_EMIF_imp <- imputePCA(t(X_EMIF_mis),ncp = optPC)$completeObs

# Save the imputed data
save(X_EMIF_imp, file = "EMIF-AD/Data/X_EMIF_imp.RData")









