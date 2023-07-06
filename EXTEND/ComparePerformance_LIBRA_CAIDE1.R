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


load("~/EXTEND/LIBRA_OutputList_Cor_EN.RData")
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$alpha == optAlpha) & 
                             (predResults[[i]]$lambda == optLambda),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.03440694
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.284027

load("~/EXTEND/LIBRA_OutputList_Cor_RF.RData")
perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$mtry == opt_mtry) & 
                             (predResults[[i]]$splitrule == opt_splitrule) &
                             (predResults[[i]]$min.node.size == opt_min.node.size),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

R2(ObsvsPred$pred, ObsvsPred$obs)
#0.03729306
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.288517

load("~/EXTEND/LIBRA_Model_Cor_None.RData")
temp <- LIBRA_Model$pred

ObsvsPred <- temp[(temp$alpha == LIBRA_Model$bestTune$alpha) &
                    (temp$lambda == LIBRA_Model$bestTune$lambda),]
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.01475733

MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.304778





load("~/EXTEND/CAIDE1_OutputList_Cor_EN.RData")
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$alpha == optAlpha) & 
                             (predResults[[i]]$lambda == optLambda),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.4470842
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.555834

load("~/EXTEND/CAIDE1_OutputList_Cor_RF.RData")
perf <- matrix(NA, nrow = nrow(trainResults[[1]]), ncol = 25)
for (j in 1:length(trainResults)){
  perf[,j] <- trainResults[[j]]$MAE
}
optPar <- which.min(rowMeans(perf))
optPerf <- perf[optPar,]
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

ObsvsPred <- NULL
for (i in 1:length(predResults)){
  temp <- predResults[[i]][(predResults[[i]]$mtry == opt_mtry) & 
                             (predResults[[i]]$splitrule == opt_splitrule) &
                             (predResults[[i]]$min.node.size == opt_min.node.size),c("pred", "obs")]
  ObsvsPred <- rbind.data.frame(ObsvsPred, temp)
}

R2(ObsvsPred$pred, ObsvsPred$obs)
#0.4650472
MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.541346

load("~/EXTEND/CAIDE1_Model_Cor_None.RData")
temp <- CAIDE1_Model$pred

ObsvsPred <- temp[(temp$alpha == CAIDE1_Model$bestTune$alpha) &
                    (temp$lambda == CAIDE1_Model$bestTune$lambda),]
R2(ObsvsPred$pred, ObsvsPred$obs)
# 0.4473198

MAE(ObsvsPred$pred, ObsvsPred$obs)
#1.566595
