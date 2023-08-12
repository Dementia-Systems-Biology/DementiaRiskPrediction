# ============================================================================ #
# File: JIVE.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform JIVE on the model's features and their mQTLs.
# ============================================================================ #

# Load packages
library(vcfR)
library(caret)
library(glmnet)
library(tidyverse)
library(r.jive)
library(httr)
library(dplyr)

# Read VCF file
vcfFile <- read.vcfR("EXTEND/Data/vcf_EXTEND.vcf")

# Prepare genotype matrix
genoType <- vcfFile@gt
genoType <- genoType[,-1]
genoType[genoType == "1/1"] <- 2
genoType[genoType == "0/0"] <- 0
genoType[genoType == "1/0"] <- 1
genoType[genoType == "0/1"] <- 1
genoType <- matrix(as.numeric(genoType), ncol = 1036)

# Set row and column names
test <- vcfFile@gt
colnames(genoType) <- colnames(test)[-1]
test <- vcfFile@fix
rownames(genoType) <- test[,3]

# Save genotype data
save(genoType,file = "EXTEND/Data/genoType.RData")

# Load data
load("EXTEND/Data/methSet_allNorm_fil.RData")
load("EXTEND/Data/genoType.RData")
load("Models/ModelInterpretation/selCpGs.RData")

# Make empty data frame for explained variance
varExpl <- as.data.frame(matrix(NA,nrow = length(selCpGs),ncol = 2))
colnames(varExpl) <- c("Methylation", "Genotype")

for (i in 1:length(selCpGs)){
  
  # Get mQTLs of model's features
  selCpGs1 <- selCpGs[[i]]
  query <- list(
    cpgs = selCpGs1,
    pval = 1e-5,
    clumped = 1
  )
  res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
  test <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
  selSNPs <- unique(str_remove(str_remove(str_remove(test$name, "chr"), ":SNP"), ":INDEL"))
  
  # Prepare data for JIVE
  samples <- intersect(colnames(methSet_allNorm_fil), colnames(genoType))
  methData <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% selCpGs1,samples]
  genoType_fil <- genoType[rownames(genoType) %in% selSNPs,samples]
  all(colnames(genoType_fil) == colnames(methData))
  
  # Make JIVE list
  jive_data <- list(Methylation = methData,
                    Genotype = genoType_fil)
  
  # Perform JIVE
  set.seed(123)
  jive_results <- jive(jive_data)
  
  # Get explained variance
  varExpl[i,] <- summary(jive_results)$Variance[1,]
  
  # Save intermediate results
  save(varExpl, file = "Models/ModelInterpretation/JIVE/varExpl.RData")
}

# Set rownames
rownames(varExpl) <- names(selCpGs)

# Prepare data for plotting
plotDF <- data.frame(Name = rep(c("Syst. blood pressure",
                              "Low education",
                              "Physical inactivity",
                              "Dietary intake",
                              "Depression",
                              "Type II diabetes",
                              "Sex",
                              "Age",
                              "Alcohol consumption",
                              "HDL cholesterol"),2),
                     Type = rep(c("Methylation", "Genotype"), each = 10),
                     Value = c(varExpl[,1], varExpl[,2]))

# Make plot
p <- ggplot(plotDF) +
  geom_bar(aes(x = Name, y = Value, fill = Type), stat = "identity",
           position = position_dodge(), color = "black") +
  coord_flip() +
  xlab("") +
  ylab("Joint variation") +
  scale_fill_manual(values = c("#E6AB02","#EF3B2C")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(p, file = "Models/ModelInterpretation/JIVE/jointVariation_JIVE.jpg", width = 7.5, height = 5)

