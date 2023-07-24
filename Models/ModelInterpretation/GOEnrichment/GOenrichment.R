library(tidyverse)
library(caret)
library(glmnet)
#BiocManager::install("missMethyl")
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
rownames(pvalues) <- firstup(paste0(gst$TERM, " (", rownames(gst), ")"))
colnames(FDRs) <- names(selCpGs)
rownames(FDRs) <- firstup(paste0(gst$TERM, " (", rownames(gst), ")"))
save(FDRs, file = "FDRs_GO.RData")
save(pvalues, file = "pvalues_GO.RData")
save(gst, file = "gst_GO.RData")

# Select BP only
load("gst_GO.RData")
gst_bp <- gst[gst$ONTOLOGY == "BP",]



simMatrix_BP <- calculateSimMatrix(rownames(gst_bp),
                                   orgdb = "org.Hs.eg.db",
                                   ont = "BP", 
                                   method = "Rel")


pvalues_id <- pvalues
rownames(pvalues_id) <- rownames(gst)
pvalues_id <- pvalues_id[rownames(gst_bp),]
reduceTerms_BP <- reduceSimMatrix(simMatrix_BP,
                                  rowMeans(-log(pvalues_id)),
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")


load("pvalues_GO.RData")
sel_terms <- firstup(paste0(gst_bp$TERM, " (", rownames(gst_bp), ")"))
sel_terms <- firstup(paste0(reduceTerms_BP$parentTerm, " (", reduceTerms_BP$parent, ")"))
sel_pvalues <- pvalues[rownames(pvalues) %in% sel_terms,]


rownames(sel_pvalues)[rowSums(sel_pvalues < 0.05) > 2]







rownames(pvalues)[which.max(rowSums(-log(pvalues)))]


# Which terms reach significance in more than 3 models?
test <- sel_pvalues[rowSums(sel_pvalues < 0.05) > 2,]
terms <- gst[rownames(test),]
colnames(test) <- c("Systolic Blood Pressure", "Low education", "Physical Inactivity",
                    "Unhealthy diet", "Depression", "Type II Diabetes", "Sex", "Age",
                   "L-M Alcohol Intake","HDL Cholesterol")

# Prepare data for plotting
plotDF <- gather(test)
plotDF$Name <- rep(rownames(test), ncol(test))
plotDF$Sig <- ifelse(plotDF$value < 0.05, "Yes", "No")

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

ggsave(p, file = "GOenrichment_heatmap.png", width = 10, height = 7)

cpgs_models <- unique(unlist(selCpGs))
set.seed(123)
gst <- gometh(sig.cpg=cpgs_models, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE,
              sig.genes = TRUE)
save(gst, file =  "gst_modelCombined_GO_new.RData")


pvalues_models <- gst$P.DE


table(unlist(selCpGs))[which(table(unlist(selCpGs))>1)]



# Exlude X and Y chromosomal proteins
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#biomaRt::listAttributes(ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "hgnc_symbol",
                                  "chromosome_name",
                                  "go_id",
                                  "name_1006"), 
                     filters = 'go',
                     values = "GO:0097113",
                     mart = ensembl)

annotations <- annotations[annotations$go_id == "GO:0097113",]

load("~/gst_modelCombined_GO_new.RData")
load("~/Data/probe_annotation.RData")
load("selCpGs.RData")
load("~/Data/Probe2Gene_final_ann.RData")

genes <- str_split(gst["GO:0097113", "SigGenesInSet"], ",")[[1]]

ann_fil <- Probe2Gene_final[Probe2Gene_final$Gene %in% genes,]
ann_fil <- ann_fil[ann_fil$Association != "Intergenic",]
ann_fil <- ann_fil[ann_fil$Association != "ExonBnd",]

selCPG_names <- c("Systolic Blood Pressure", "Low education", "Physical Inactivity",
                  "Unhealthy diet", "Depression", "Type II Diabetes", "Sex", "Age",
                  "L-M Alcohol Intake","HDL Cholesterol")
AssociatedGenes <- NULL
for (i in 1:length(selCpGs)){
  AssociatedGenes_i <- ann_fil[ann_fil$CpG %in% selCpGs[[i]],]
  if (nrow(AssociatedGenes_i) != 0){
    AssociatedGenes_i$Model <- selCPG_names[i]
    AssociatedGenes <- rbind.data.frame(AssociatedGenes, AssociatedGenes_i)
  }
}


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

ggsave(p, file = "AMPAclustering.jpg", width = 8, height = 5)


gst_bp <- gst[gst$ONTOLOGY == "BP",]
gst_bp$FDR <- p.adjust(gst_bp$P.DE, method = "fdr")
gst_bp <- arrange(gst_bp, by = P.DE)
gst_bp$SigGenesInSet <- str_replace_all(gst_bp$SigGenesInSet,",", ";")

write.csv(gst_bp, file = "GOanalysis_all.csv",quote = FALSE)
write.csv(head(gst_bp,20), file = "GOanalysis_fil.csv",quote = FALSE)
