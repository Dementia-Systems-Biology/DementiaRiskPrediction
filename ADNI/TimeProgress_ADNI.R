# ============================================================================ #
# File: TimeAnalysis_ADNI.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Calculate the MRSs in the ADNI cohort.
# ============================================================================ #

# Load packages
library(tidyverse)
library(caret)
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
                 position = "identity",color = "black", bins= 40) +
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

var <- "MMSE"
testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF$MMSE),]
testDF <- testDF[testDF$MMSE < 24,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF$MMSE < 24), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     MMSE = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = "ADNI/TimeAnalysis/KaplanMeier_ADNI_MMSE.jpg", width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT learning

###############################################################################

var <- "RAVLT.learning"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]

normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     RAVLT.learning = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)


kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT forgetting

###############################################################################

var <- "RAVLT.forgetting"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     RAVLT.forgetting = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT percent forgetting

###############################################################################

var <- "RAVLT.perc.forgetting"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     RAVLT.perc.forgetting = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)

###############################################################################

# RAVLT immediate

###############################################################################

var <- "RAVLT.immediate"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     RAVLT.immediate = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)


kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADASQ4

###############################################################################

var <- "ADASQ4"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     ADASQ4 = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)


kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)
###############################################################################

# ADAS13

###############################################################################

var <- "ADAS13"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     ADAS13 = NA,
                     Time = 85,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 


kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)
###############################################################################

# ADAS11

###############################################################################

var <- "ADAS11"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     ADAS11 = NA,
                     Time = max(testDF$Time) + 1,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)
###############################################################################

# LDELTOTAL

###############################################################################

var <- "LDELTOTAL"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var])
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var])

testDF <- testDF[testDF[,var] < mean_var - 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     LDELTOTAL = NA,
                     Time = max(testDF$Time) + 1,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)

###############################################################################

# TRABSCOR: Trails B

###############################################################################

var <- "TRABSCOR"

testDF <- metaData_fil_time[,c("RID", "VISCODE", var)]
testDF <- testDF[!is.na(testDF[,var]),]

sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)

testDF <- testDF[testDF[,var] > mean_var + 2*sd_var,]
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

for (i in unique(testDF$RID)){
  testDF[testDF$RID == i, "Time"] <- min(testDF[testDF$RID == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,c(1,4,5)]),]
normal <- data.frame(RID = setdiff(metaData_fil$RID[!is.na(metaData_fil[,var])], unique(testDF$RID)),
                     VISCODE = "m85",
                     TRABSCOR = NA,
                     Time = max(testDF$Time) + 1,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(testDF$RID))
length(unique(metaData_fil_time$RID))


kaplanDF <- testDF[,c("RID", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("RID" = "RID"))
kaplanDF$Time <- kaplanDF$Time/12


survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)


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

ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,".jpg"), width = 7, height = 5)

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

coxph(Surv(Time, Test) ~ predClass, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 


kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)