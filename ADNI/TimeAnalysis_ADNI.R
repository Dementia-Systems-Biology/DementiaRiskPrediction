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
library(tidyverse)
library(caret)
library(tidyverse)
library(caret)
library(patchwork)
library(ranger)
library(pROC)

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

# Filter for 
keep_samples <- unique(MetaData_allTime$RID[(MetaData_allTime$VISCODE == "m48") | (MetaData_allTime$VISCODE == "m42")])
metaData_fil <- metaData_fil[metaData_fil$RID %in% keep_samples,]

# Filter for MMSE >= 26
metaData_fil <- metaData_fil[metaData_fil$MMSE.bl >= 26,]

# prepare data
predictedScore_factors_fil <- predictedScore_factors[metaData_fil$Basename,]
table(metaData_fil$DX)

# Meta data for all time points
metaData_fil_time <- MetaData_allTime[MetaData_allTime$RID %in% metaData_fil$RID,]

# Make predictions using epi-MCI (RF-RFE) model
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")
predictDF <- data.frame(RID = metaData_fil$RID, 
                        pred = pred_RF$MCI,
                        Age = metaData_fil$Age)

# Get mean prediction for same individual 
# (some individuals have more than one baseline sample available)
predictDF <- predictDF %>%
  group_by(RID) %>%
  reframe(RID = RID,
          pred= mean(pred),
          Age = Age)

predictDF <- unique(predictDF)


# Epi-MCI (RF-RFE) score as predictor
predictDF$predClass <- "Intermediate risk"
predictDF$predClass[predictDF$pred < quantile(predictDF$pred,0.33)] <- "Low risk"
predictDF$predClass[predictDF$pred > quantile(predictDF$pred,0.67)] <- "High risk"
table(predictDF$predClass)


# Histogram of predictions
predictDF$predClass <- factor(predictDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- ggplot(predictDF) + 
  geom_histogram(aes(x = log(pred/(1-pred)), fill = predClass),
                 position = "identity",color = "black", bins= 35) +
  theme_classic() +
  xlab("Epi-MCI Score") +
  ylab("Count") +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(p, file = "ADNI/TimeAnalysis/RiskScoreDistribution_ADNI.png", width = 7, height = 5)


###############################################################################

# MMSE

###############################################################################

# Prepare data for survival analysis
var <- "MMSE"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF$MMSE),]

# Only include cognitive impaired
testDF <- testDF[testDF$MMSE < 24,]
testDF$Status <- ifelse((testDF$MMSE < 24), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]


# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     MMSE = NA,
                     Time = 49,
                     Status = "Normal")

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))


# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("MMSE") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "ADNI/TimeAnalysis/KaplanMeier_ADNI_MMSE.jpg", width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT learning

###############################################################################

# Prepare data for survival analysis
var <- "RAVLT.learning"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")
  
# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (learning)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT forgetting

###############################################################################

# Prepare data for survival analysis
var <- "RAVLT.forgetting"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")
  
# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (forgetting)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 


# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT percent forgetting

###############################################################################

# Prepare data for survival analysis
var <- "RAVLT.perc.forgetting"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (percent forgetting)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)

###############################################################################

# RAVLT immediate

###############################################################################

# Prepare data for survival analysis
var <- "RAVLT.immediate"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (immediate)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADASQ4

###############################################################################

# Prepare data for survival analysis
var <- "ADASQ4"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("ADAS-Q4") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADAS13

###############################################################################

# Prepare data for survival analysis
var <- "ADAS13"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle(var) +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADAS11

###############################################################################

# Prepare data for survival analysis
var <- "ADAS11"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle(var) +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# Log rank test
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# LDELTOTAL

###############################################################################

# Prepare data for survival analysis
var <- "LDELTOTAL"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("Logical Memory - Delayed Recall") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)

###############################################################################

# TRABSCOR: Trails B

###############################################################################

# Prepare data for survival analysis
var <- "TRABSCOR"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]

# remove NA values
testDF <- testDF[!is.na(testDF[,var]),]

# Only include cognitive impaired
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Only include up to 48 months in analysis
testDF <- testDF[testDF$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
colnames(testDF) <- c("RID", "VISCODE", "Var", "Status","Time")

# Set the samples not converted to a cognitively healthy status at month 49
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m49",
                     Var = NA,
                     Status = "Normal",
                     Time = 49)

testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID)) == length(unique(metaData_fil_time$RID))

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("RID", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("Trail Making Test Part B Time") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)