# ============================================================================ #
# File: TimeAnalysis_ADNI.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform survival (time) analysis in the ADNI cohort.
# ============================================================================ #

# Load packages
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(tidyverse)
library(caret)
library(patchwork)
library(ranger)
library(pROC)
library(mgcv)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("ADNI/Data/predictedScore_factors_ADNI.RData")
load("ADNI/Data/MetaData_ADNI.RData")

# Filter for midlife samples
metaData_fil <- as.data.frame(MetaData_baseline)
rownames(metaData_fil) <- metaData_fil$Basename
midlife_samples <- intersect(metaData_fil$Basename[metaData_fil$Age <= 75], 
                             rownames(predictedScore_factors))
metaData_fil <- metaData_fil[midlife_samples,]


#metaData_fil <- metaData_fil[metaData_fil$MMSE>= 26,]

# prepare data
predictedScore_factors_fil <- predictedScore_factors[metaData_fil$Basename,]
table(metaData_fil$DX)

# Make predictions using epi-MCI (RF-RFE) model
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")

vars <- c("MMSE", "RAVLT.learning", "RAVLT.immediate", "RAVLT.forgetting",
          "RAVLT.perc.forgetting", "ADASQ4","ADAS11","ADAS13","TRABSCOR",
          "TAU","PTAU", "ABETA")
coeffs <- rep(NA, length(vars))
pvalues <- rep(NA, length(vars))

for (v in 1:length(vars)){
  predictDF <- data.frame(RID = metaData_fil$RID, 
                          pred = log(pred_RF$MCI/(1-pred_RF$MCI)),
                          Age = metaData_fil$Age,
                          Sex = ifelse(metaData_fil$Sex== "M",1,0),
                          Score = metaData_fil[,vars[v]])
  
  predictDF <- predictDF[!is.na(predictDF$Score),]
  if (vars[v] == "ABETA"){
    predictDF$Score[predictDF$Score == ">1700"] <- "1701"
  }
  if (vars[v] == "TAU" | vars[v] == "PTAU"){
    predictDF$Score <- predictDF$Score
  }
  
  
  # Get mean prediction for same individual 
  # (some individuals have more than one baseline sample available)
  predictDF <- predictDF %>%
    group_by(RID) %>%
    reframe(RID = RID,
            pred= mean(pred),
            Age = Age,
            Sex = Sex,
            Score = Score)
  
  predictDF <- unique(predictDF)
  
  coeffs[v] <- summary(lm(Score ~ pred + Age + Sex, data = predictDF))$coefficients["pred", "Estimate"]
  pvalues[v] <- summary(lm(Score ~ pred + Age + Sex, data = predictDF))$coefficients["pred", "Pr(>|t|)"]
}

names(pvalues) <- vars
names(coeffs) <- vars
p.adjust(pvalues, method = "fdr")
sign(coeffs)
