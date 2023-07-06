library(tidyverse)
library(caret)
library(glmnet)
#BiocManager::install("missMethyl")
library(missMethyl)

# Function to capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Clear workspace and console
rm(list = ls())
cat("\014") 


load("~/Data/Fit_EMIF_MCI_RF.RData")
load("Data/methSet_allNorm_fil.RData")

 


load("EXTEND/finalModels.RData")

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


cpgs <- read.csv("Data/Predictors_Shiny_by_Groups.csv", header = T) 

selCpGs[[8]] <- cpgs$CpG_Site[cpgs$Predictor == "Epigenetic Age (Zhang)"]
selCpGs[[9]] <- cpgs$CpG_Site[cpgs$Predictor == "Alcohol"]
selCpGs[[10]] <- cpgs$CpG_Site[cpgs$Predictor == "HDL Cholesterol"]
names(selCpGs) <- rownames(varImp(fit)$importance)
save(selCpGs, file = "selCpGs.RData")


allCpGs <- rownames(methSet_allNorm_fil)

pvalues <- matrix(NA, nrow = 22708, ncol = length(selCpGs))
FDRs <- pvalues
for (i in 1:length(selCpGs)){
  
  gst <- gometh(sig.cpg=selCpGs[[i]], all.cpg=allCpGs, collection="GO",
                array.type = "EPIC",
                plot.bias=TRUE)
  
  pvalues[,i] <- gst$P.DE
  FDRs[,i] <- gst$FDR
  
}

pvalues <- as.data.frame(pvalues)
FDRs <- as.data.frame(FDRs)
colnames(pvalues) <-  names(selCpGs)
rownames(pvalues) <- firstup(paste0(gst$TERM, " (", gst$ONTOLOGY, ")"))
colnames(FDRs) <- names(selCpGs)
rownames(FDRs) <- firstup(paste0(gst$TERM, " (", gst$ONTOLOGY, ")"))
save(FDRs, file = "FDRs_GO.RData")
save(pvalues, file = "pvalues_GO.RData")
save(gst, file = "gst_GO.RData")


load("pvalues_GO.RData")
rownames(pvalues)[rowSums(pvalues < 0.05) > 2]
rownames(pvalues)[which.max(rowSums(-log(pvalues)))]


# Which terms reach significance in more than 3 models?
test <- pvalues[rowSums(pvalues < 0.1) > 3,]
terms <- gst[rownames(test),]
colnames(test) <- c("Systolic Blood Pressure", "Low education", "Physical Inactivity",
                    "Unhealthy diet", "Depression", "Type II Diabetes", "Sex", "Age",
                   "L-M Alcohol Intake","HDL Cholesterol")

# Prepare data for plotting
plotDF <- gather(test)
plotDF$Name <- rep(rownames(test), ncol(test))
plotDF$Sig <- ifelse(plotDF$value < 0.01, "Yes", "No")

# Perform hierarchical clustering on the terms
clusters <- hclust(dist(-log(test)), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

# Make plot (heatmap)
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log(value), color = Sig),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_viridis_c() +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = "-log p-value") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")



cpgs_models <- unique(unlist(selCpGs))
gst <- gometh(sig.cpg=cpgs_models, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE)
save(gst, file =  "gst_modelCombined_GO.RData")

gst <- gometh(sig.cpg=cpgs_models, all.cpg=allCpGs, collection="KEGG",
              array.type = "EPIC",
              plot.bias=TRUE)
save(gst, file =  "gst_modelCombined_KEGG.RData")

pvalues_models <- gst$P.DE
