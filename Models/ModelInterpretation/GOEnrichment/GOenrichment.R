# ============================================================================ #
# File: GOenrichment.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform GO overrepresentation analysis on the model's features.
# ============================================================================ #

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(missMethyl)
library(rrvgo)

# Function to capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load Models
load("Models/MRS_Models/finalModels.RData")          # MRS models
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData") # Epi-MCI

# Get models' features (CpGs):

# EXTEND Models
factors <- rownames(varImp(fit)$importance)[1:7]
selCpGs <- list()
for (i in 1:length(factors)){
  if (finalModels[[factors[i]]]$method == "glmnet"){
    varImportance <- varImp(finalModels[[factors[i]]])$importance
    selCpGs[[i]] <- rownames(varImportance)[varImportance[,1] != 0]
  } else{
    varImportance <- arrange(varImp(finalModels[factors[i]])$importance, by = Overall)
    selCpGs[[i]] <- tail(rownames(varImportance),1000)
  }

}

names(selCpGs) <- factors

# Marioni Models
cpgs <- read.csv("Data/Predictors_Shiny_by_Groups.csv", header = T) 
selCpGs[[8]] <- cpgs$CpG_Site[cpgs$Predictor == "Epigenetic Age (Zhang)"]
selCpGs[[9]] <- cpgs$CpG_Site[cpgs$Predictor == "Alcohol"]
selCpGs[[10]] <- cpgs$CpG_Site[cpgs$Predictor == "HDL Cholesterol"]
names(selCpGs) <- rownames(varImp(fit)$importance)

# Save model's features (CpGs)
save(selCpGs, file = "Models/ModelInterpretation/selCpGs.RData")

# Get all possible CpGs
load("EXTEND/Data/methSet_allNorm_fil.RData")
allCpGs <- rownames(methSet_allNorm_fil)

# Get the union of all model's CpGs
cpgs_models <- unique(unlist(selCpGs))

# Perform GO overrepresentation analysis
set.seed(123)
gst <- gometh(sig.cpg=cpgs_models, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE,
              sig.genes = TRUE)

# Save results
save(gst, file ="Models/ModelInterpretation/GOEnrichment/gst_modelCombined_GO_new.RData")

#******************************************************************************#
# Get genes and CpGs associated with most significantly enriched GO term
#******************************************************************************#

# Load data
load("Models/ModelInterpretation/GOEnrichment/gst_modelCombined_GO_new.RData")
load("Models/ModelInterpretation/selCpGs.RData")
load("Models/ModelInterpretation/Probe2Gene_final_ann.RData")

# Get associated genes
genes <- str_split(gst["GO:0097113", "SigGenesInSet"], ",")[[1]]

# Get associated CpGs
ann_fil <- Probe2Gene_final[Probe2Gene_final$Gene %in% genes,]
ann_fil <- ann_fil[ann_fil$Association != "Intergenic",]
ann_fil <- ann_fil[ann_fil$Association != "ExonBnd",]

# Set risk factor names
selCPG_names <- c("Systolic Blood Pressure", "Low education", "Physical Inactivity",
                  "Unhealthy diet", "Depression", "Type II Diabetes", "Sex", "Age",
                  "L-M Alcohol Intake","HDL Cholesterol")

# For each risk factor model, get the CpGs associated with the GO term
AssociatedGenes <- NULL
for (i in 1:length(selCpGs)){
  AssociatedGenes_i <- ann_fil[ann_fil$CpG %in% selCpGs[[i]],]
  if (nrow(AssociatedGenes_i) != 0){
    AssociatedGenes_i$Model <- selCPG_names[i]
    AssociatedGenes <- rbind.data.frame(AssociatedGenes, AssociatedGenes_i)
  }
}

# Make plot
p <- ggplot(AssociatedGenes)+
  geom_tile(aes(x = CpG, y = Model, fill = Association)) +
  facet_grid(cols = vars(Gene), scales = "free", space = "free") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("") +
  ylab("") +
  ggtitle("AMPA glutamate receptor clustering (GO:0097113)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

# Save plot
ggsave(p, file = "AMPAclustering.jpg", width = 8, height = 5)

# Get statistics
gst_bp <- gst[gst$ONTOLOGY == "BP",]
gst_bp$FDR <- p.adjust(gst_bp$P.DE, method = "fdr")
gst_bp <- arrange(gst_bp, by = P.DE)
gst_bp$SigGenesInSet <- str_replace_all(gst_bp$SigGenesInSet,",", ";")

# Save statistics
write.csv(gst_bp, file = "GOanalysis_all.csv",quote = FALSE)
write.csv(head(gst_bp,20), file = "GOanalysis_fil.csv",quote = FALSE)
