# Load packages
library(prospectr)
library(tidyverse)
library(caret)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_fil.RData")
load("EMIF/X_EMIF_imp.RData")
load("~/Data/LIBRA_Model_Cor_RF.RData")
load("~/Data/CAIDE1_Model_Cor_RF.RData")



epiCAIDE <- predict(CAIDE1_Model, X_EMIF_imp[metaData_fil$X,])
epiLIBRA <- predict(LIBRA_Model, X_EMIF_imp[metaData_fil$X,])

load("~/EMIF/Y_test_EMIF.RData")
metaData_fil <- Y_test
epiCAIDE <- predict(CAIDE1_Model, X_EMIF_imp[metaData_fil$X,])
epiLIBRA <- predict(LIBRA_Model, X_EMIF_imp[metaData_fil$X,])

plot(epiCAIDE, metaData_fil$Age)

plotDF <- data.frame(epiCAIDE = epiCAIDE,
                     epiLIBRA = epiLIBRA,
                     Y = metaData_fil$Diagnosis)

plotDF <- plotDF[(plotDF$Y == "MCI") | (plotDF$Y == "NL"),]
plotDF$Y <- factor(ifelse(plotDF$Y == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

# Calculate sensitivites, specificities and AUC values
score <- c("epiCAIDE", "epiLIBRA")
scoreName1 <- c("Epi-CAIDE:", "Epi-LIBRA:")
ROCplot <- NULL                       # Data frame with sensitivities and specificities
aucValue <- rep(NA, length(score))    # AUC
liValue <- rep(NA, length(score))     # lower interval value of AUC
uiValue <- rep(NA, length(score))     # upper interval value of AUC

for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])#, direction = "<")
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName1[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(test)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(test)[3]),2),nsmall = 2)
}

# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0(scoreName1,"\t",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName1,
                      X = 0.75,
                      Y = rev(seq(0.05,0.1,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName1)

# Colors for plotting
colors <- rev(c("#084594","#EF3B2C","#CB181D", "#99000D"))

colors <- c(rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold", size = 4) +
  scale_color_manual(values = colors) +
  ggtitle("MCI vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_CAIDE_LIBRA_MCI_rev.png", width = 7.5, height = 5)



plotDF <- data.frame(epiCAIDE = epiCAIDE,
                     epiLIBRA = epiLIBRA,
                     Y = metaData_fil$Diagnosis)

plotDF <- plotDF[(plotDF$Y == "AD") | (plotDF$Y == "NL"),]
plotDF$Y <- factor(ifelse(plotDF$Y == "NL","Control","AD"),
                   levels = c("Control", "AD"))

# Calculate sensitivites, specificities and AUC values
score <- c("epiCAIDE", "epiLIBRA")
scoreName1 <- c("Epi-CAIDE:", "Epi-LIBRA:")
ROCplot <- NULL                       # Data frame with sensitivities and specificities
aucValue <- rep(NA, length(score))    # AUC
liValue <- rep(NA, length(score))     # lower interval value of AUC
uiValue <- rep(NA, length(score))     # upper interval value of AUC

for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])#, direction = "<")
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName1[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
  liValue[i] <- format(round(as.numeric(ci(test)[1]),2),nsmall = 2)
  uiValue[i] <- format(round(as.numeric(ci(test)[3]),2),nsmall = 2)
}

# Combine AUC values into data frame
plotAUC <- data.frame(AUC = paste0(scoreName1,"\t",aucValue, " (", liValue, "-", uiValue, ")"),
                      Score = scoreName1,
                      X = 0.75,
                      Y = rev(seq(0.05,0.1,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName1)

# Colors for plotting
colors <- rev(c("#084594","#EF3B2C","#CB181D", "#99000D"))

colors <- c(rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold", size = 4) +
  scale_color_manual(values = colors) +
  ggtitle("AD vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_CAIDE_LIBRA_AD_rev.png", width = 7.5, height = 5)