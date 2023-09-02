# ============================================================================ #
# File: ValidateMPSs_EMIF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Validate the performance of the MRSs in the EMIF-AD cohort.
# ============================================================================ #

# Load packages
library(tidyverse)
library(caret)
library(pROC)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/metaData_fil.RData")                # Meta data
load("EMIF-AD/Data/predictedScore_factors_fil.RData")  # MRSs

# Prepare data
PhenoNames <- c("Hypertension", "Alcohol Consumption", "Cardiovascular Disease",
                "Depression", "Years of Education", "Smoking Status", "Obesity")
testDF <- metaData_fil[,c("Hypertension", "Alcohol", "CardiovascularDis",
                          "Depression", "Eduy", "Smoking", "Obesity")]

MRS_test <- colnames(predictedScore_factors_fil)
MRS_names <- c("Systolic blood pressure",
               "Total cholesterol",
               "Low education",
               "Phsyical inactivity",
               "Unhealthy diet",
               "Depression",
               "Type II diabetes",
               "Heart disease",
               "Sex (male)",
               "Age",
               "Alcohol intake",
               "BMI",
               "HDL cholesterol",
               "Smoking")

# Get pairwise correlations
corDF <- NULL
for (i in 1:ncol(testDF)){
  temp_cor <- apply(predictedScore_factors_fil[,MRS_test],2,function(x){cor(x,as.numeric(testDF[,i]),
                                                     method = "spearman", use = "pairwise.complete.obs")})
  temp_p <- apply(predictedScore_factors_fil[,MRS_test],2,function(x){cor.test(x,as.numeric(testDF[,i]),
                                                     method = "spearman", use = "pairwise.complete.obs")$p.value})
  
  temp <- data.frame(MRS = MRS_names,
                     Pheno = rep(PhenoNames[i],length(temp_cor)),
                     Cor = temp_cor,
                     Pvalue = temp_p)
  corDF <- rbind.data.frame(corDF, temp)
}

corDF$FDR <- p.adjust(corDF$Pvalue, method = "fdr")
corDF$Sig <- ifelse(corDF$FDR< 0.05, "Yes", "No")

# Make correlation plot
p <- ggplot(corDF) +
  geom_tile(aes(x = MRS, y = Pheno, fill = Cor, color = Sig),
            stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  xlab("Methylation Risk Scores") +
  ylab("Observed Phenotype") +
  labs(fill  = "Spearman\nCorrelation")+
  #coord_flip() +
  theme_bw()+
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       oob = scales::squish,
                       limits = c(-0.5,0.5))+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(color = "none")


# Save plot
ggsave(p, file = "EMIF-AD/ModelPerformance/EMIF_Correlations_MRS.png", width = 7, height = 6)


# Prepare data for AUROC calculation
temp <- metaData_fil
temp$Alcohol <- ifelse(temp$Smoking == 0, 0,1)
temp$Smoking <- ifelse(temp$Smoking == 0, 0,1)
temp$Eduy <- ifelse(temp$Eduy <= 9,1,0)
test <- data.frame(Pheno = c("Depression", "Hypertension", "Alcohol", "Smoking", "Eduy", "CardiovascularDis", "Obesity"),
                   MRS = c("Depression", "SysBP", "Alcohol", "Smoking", "Education", "HeartDisease", "BMI"))

ROCplot <- NULL
aucValue <- rep(NA, nrow(test))    # AUC
liValue <- rep(NA, nrow(test))     # lower interval value of AUC
uiValue <- rep(NA, nrow(test))     # upper interval value of AUC


# Get AUROC for each risk factor
factorName <- c("Depression", "Hypertension", "Alcohol Consumption", "Smoking Status", "Education", "Cardiovascular Disease", "Obesity")
for (i in 1:nrow(test)){
  r <- pROC::roc(temp[,test$Pheno[i]], predictedScore_factors_fil[, test$MRS[i]])
  temp_roc <- data.frame(Sensitivity = r$sensitivities,
                         Specificity = r$specificities,
                         Name = rep(factorName[i],length(r$specificities)))
  aucValue[i] <- format(round(as.numeric(auc(r)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(r)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(r)[3]),2),nsmall = 2)
  
  ROCplot <- rbind.data.frame(ROCplot, temp_roc)
}
# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0("AUROC = ",aucValue, " (", liValue, "-", uiValue, ")"),
                      Name = factorName,
                      X = 0.8,
                      Y = rev(seq(0.05,0.3,length.out = length(aucValue))))

ROCplot$Name <- factor(ROCplot$Name, levels = factorName)

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Name), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Name),
            fontface = "bold") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "EMIF-AD/ModelPerformance/ROC_factors_EMIF.png", width = 8, height = 5)
