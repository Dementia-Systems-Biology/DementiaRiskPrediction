
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
load("~/PPMI/lumi_dpval_PPMI.RData")
load("~/PPMI/methSet_allNorm_PPMI.RData")
load("~/PPMI/metaData_ppmi.RData")

X_PPMI <- methSet_allNorm
all(colnames(X_PPMI) == colnames(lumi_dpval))

# optimize imputation
X_PPMI_M <- log2(X_PPMI/(1-X_PPMI))
test <- lumi_dpval[rownames(X_PPMI_M), colnames(X_PPMI_M)]
X_PPMI_mis <- X_PPMI_M
X_PPMI_mis[test > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_PPMI_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_PPMI_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_PPMI_copy <- X_PPMI_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_PPMI_mis[removeProbes[j,1], removeProbes[j,2]]
  X_PPMI_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}



nPCs_CV <- 10
pred <- matrix(NA, nrow = 100,ncol = nPCs_CV)
for (p in 3:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_PPMI_imp_CV <- imputePCA(t(X_PPMI_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_PPMI_imp_CV[removeProbes[j,2], removeProbes[j,1]]
  }
  save(pred, file = "PPMI/pred_imp.RData")
}

# Evaluate which number of PCs is the best
load("PPMI/pred_imp.RData")
test <- apply(pred,2,function(x) MAE(obs = values,pred = x))
optPC <- which.min(test) # 5 PCs

# Perform imputation using optimal number of PCs
set.seed(123)
X_PPMI_imp <- imputePCA(t(X_PPMI_mis),ncp = optPC)$completeObs
save(X_PPMI_imp, file = "PPMI/X_PPMI_imp.RData")









