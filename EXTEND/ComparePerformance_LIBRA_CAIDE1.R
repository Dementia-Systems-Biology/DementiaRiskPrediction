# ============================================================================ #
# File: ComparePerformance_LIBRA_CAIDE1.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Get cross-validation performance of LIBRA and CAIDE1 prediction.
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

#==============================================================================#
# LIBRA model: correlation + EN:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/LIBRA_OutputList_Cor_EN.RData")

# Get CV performance
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}

# Optimal hyperparameters
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

# Observed vs predicted value in CV
ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$alpha == optAlpha) & 
                             (predResults[[i]]$lambda == optLambda),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.03440694

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.284027


#==============================================================================#
# LIBRA model: correlation + RF:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/LIBRA_OutputList_Cor_RF.RData")

# Get CV performance
perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}

# Optimal hyperparameters
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

# Observed vs predicted value in CV
ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$mtry == opt_mtry) & 
                             (predResults[[i]]$splitrule == opt_splitrule) &
                             (predResults[[i]]$min.node.size == opt_min.node.size),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
#0.03729306

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.288517

#==============================================================================#
# LIBRA model: no feature selection + EN:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/LIBRA_Model_None_EN.RData")

# CV performance
temp <- LIBRA_Model$pred

# Observed vs predicted value in CV
ObsvsPred <- temp[(temp$alpha == LIBRA_Model$bestTune$alpha) &
                    (temp$lambda == LIBRA_Model$bestTune$lambda),]

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.01475733

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.304778

#==============================================================================#
# CAIDE1 model: correlation + EN:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/CAIDE1_OutputList_Cor_EN.RData")

# Get CV performance
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}

# Optimal hyperparameters
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

# Observed vs predicted value in CV
ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$alpha == optAlpha) & 
                             (predResults[[i]]$lambda == optLambda),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.4470842

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.555834

#==============================================================================#
# CAIDE1 model: correlation + RF:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/CAIDE1_OutputList_Cor_RF.RData")

# Get CV performance
perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}

# Optimal hyperparameters
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

# Observed vs predicted value in CV
ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$mtry == opt_mtry) & 
                             (predResults[[i]]$splitrule == opt_splitrule) &
                             (predResults[[i]]$min.node.size == opt_min.node.size),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
#0.4650472

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.541346

#==============================================================================#
# CAIDE1 model: no feature selection + EN:
#==============================================================================#

# Load model
load("Models/LIBRA_CAIDE/CAIDE1_Model_Cor_None.RData")

# CV performance
temp <- CAIDE1_Model$pred

# Observed vs predicted value in CV
ObsvsPred <- temp[(temp$alpha == CAIDE1_Model$bestTune$alpha) &
                    (temp$lambda == CAIDE1_Model$bestTune$lambda),]

# Get R-squared
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.4473198

# Get MAE
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.566595
