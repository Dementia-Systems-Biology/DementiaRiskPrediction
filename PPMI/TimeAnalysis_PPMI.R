# ============================================================================ #
# File: TimeAnalysis_PPMI.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform survival (time) analysis in the PPMI cohort.
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

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
cogcat <- read.csv("PPMI/Data/CogCatWide_Filtered.csv") # Diagnosis over time
cogcat$PATNO <- as.character(cogcat$PATNO)         
load("PPMI/Data/predictedScore_factors_PPMI.RData")     # MRSs
load("PPMI/Data/metaData_ppmi.RData")                   # Meta data

# Remove reverters
reverters <- cogcat
reverters[is.na(reverters)] <- ""

rev_cat <- ""
for (i in 3:11){
  rev_cat <- paste0(rev_cat, reverters[,i]) 
}

reverters <- c(cogcat$PATNO[str_detect(rev_cat, "DementiaNormal")],
               cogcat$PATNO[str_detect(rev_cat, "MCINormal")],
               cogcat$PATNO[str_detect(rev_cat, "DementiaCognitive Complaint")],
               cogcat$PATNO[str_detect(rev_cat, "MCICognitive Complaint")])

cogcat <- cogcat[!(cogcat$PATNO %in% reverters),]

# midlife samples
metaData_all <- metaData_all[(metaData_all$age >= 40) & (metaData_all$age <= 75),]
samples <- intersect(metaData_all$Basename, rownames(predictedScore_factors))
rownames(metaData_all) <- metaData_all$Basename
metaData_fil <- metaData_all[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]

# combine meta data with cogcat
test <- inner_join(metaData_fil, cogcat, by = c("PATNO" = "PATNO"))
table(test$Class, test$CogDecon)
predictedScore_factors_fil <- predictedScore_factors_fil[test$Basename,]

# make predictions:

# Load epi-MCI model (RF-RFE)
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData")

# Make predictions
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")
predictDF <- data.frame(PATNO = test$PATNO, 
                        pred = pred_RF$MCI)

# Split into three classes
predictDF$predClass <- "Intermediate risk"
predictDF$predClass[predictDF$pred < quantile(predictDF$pred,0.33)] <- "Low risk"
predictDF$predClass[predictDF$pred > quantile(predictDF$pred,0.67)] <- "High risk"
table(predictDF$predClass)

# Make histogram of scores
predictDF$predClass <- factor(predictDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- ggplot(predictDF) + 
  geom_histogram(aes(x = log(pred/(1-pred)), fill = predClass),
                 position = "identity",color = "black", bins= 30) +
  theme_classic() +
  xlab("Epi-MCI Score") +
  ylab("Count") +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(p, file = "PPMI/TimeAnalysis/RiskScoreDistribution.png", width = 7, height = 5)

# Prepare data into correct format for survival analysis
testDF <- gather(as.data.frame(cogcat[,3:11]))
testDF$PATNO <- rep(cogcat$PATNO, 9)
testDF <- testDF[!is.na(testDF$value),]
testDF <- testDF[(testDF$value == "Dementia") | (testDF$value == "MCI"),]
testDF$Time <- as.numeric(str_remove(testDF$key, "X"))
testDF$Status <- ifelse((testDF$value == "Dementia") | (testDF$value == "MCI"), "Cognitive Impaired", "Normal")

for (i in unique(testDF$PATNO)){
  testDF[testDF$PATNO == i, "Time"] <- min(testDF[testDF$PATNO == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,3:5]),]
normal <- data.frame(key = "X9",
                     value = "Normal",
                     PATNO = setdiff(cogcat$PATNO, unique(testDF$PATNO)),
                     Time = 9,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(cogcat$PATNO))


# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF[,c("PATNO", "Time", "Status")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("PATNO" = "PATNO"))

# Make kaplan-meier curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FDD0A2","#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("Dementia/MCI") +
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
ggsave(p, file = "KaplanMeier_PPMI_MRSonly.jpg", width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age, 
         data = kaplanDF) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF$Status1 <- factor(kaplanDF$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF$Status1, kaplanDF$pred)
ci(test_roc)
