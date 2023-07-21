# Clear workspace and console
rm(list = ls())
cat("\014") 

load("~/PGS_EMIF_AD.RData")
load("~/predictedScore_factors_EMIF.RData")
load("~/metaData_fil.RData")

samples <- intersect(PGS_all$ID, metaData_fil$Sample_Name)
rownames(metaData_fil) <- metaData_fil$Sample_Name

rownames(PGS_all) <- PGS_all$ID
PGS_fil <- PGS_all[samples,]
MRS_fil <- predictedScore_factors[metaData_fil[samples, "X"],]

plotDF <- as.data.frame(matrix(NA,nrow = 11, ncol = 4))
colnames(plotDF) <- c("MRS", "PGS", "Cor", "pvalue")

MRS_names <- c("SysBP", "TotalChol", "Education", "Physical", "Diet",
               "Depression", "Diabetes", "HeartDisease", "Alcohol", 
               "BMI", "HDL")
PGS_names <- c("SBPauto", "TC", "EA22", "MVPA", "DC2",
               "MDD", "T2D", "CAD", "Alcohol", 
               "BMI", "HDL")

for (i in 1:length(MRS_names)){
  correl <- cor.test(PGS_fil[,PGS_names[i]], MRS_fil[,MRS_names[i]])
  pvalue <- as.numeric(correl$p.value)
  coeff <- as.numeric(correl$estimate)
  
  plotDF[i,] <- c(MRS_names[i], PGS_names[i], coeff,pvalue)
}

plotDF$pvalue <- as.numeric(plotDF$pvalue)

names(selCpGs)
plotDF <- plotDF[plotDF$MRS %in% names(selCpGs),]

plotDF$Name <- c("Syst. blood pressure",
                 "Low education",
                 "Physical inactivity",
                 "Dietary intake",
                 "Depression",
                 "Type II diabetes",
                 "Alcohol consumption",
                 "HDL cholesterol")

plotDF$Sig <- ifelse(plotDF$pvalue < 0.05, "Yes", "No")

p <- ggplot(plotDF) +
  geom_bar(aes(x = Name, y = abs(as.numeric(Cor)), fill = Sig), stat = "identity",
           position = position_dodge(), color = "black") +
  coord_flip() +
  xlab("") +
  ylab("|Pearson correlation coefficient|") +
  labs(fill = "p-value < 0.05") +
  scale_fill_manual(values = c("grey","#EF3B2C")) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(p, file = "PGS_MRS_Cor.png", width = 7.5, height = 5)
