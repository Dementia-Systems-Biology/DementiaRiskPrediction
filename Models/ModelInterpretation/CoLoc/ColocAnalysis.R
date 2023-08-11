# ============================================================================ #
# File: ColocAnalysis.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform colocalization analysis.
# ============================================================================ #

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(coloc)
library(data.table)
library(httr)
library(ggrepel)

# Load GWAS summary statistics
AD_summary <- fread("Data/ADFamilyHistory_2019_84.txt")
AD_summary <- as.data.frame(AD_summary)
rownames(AD_summary) <- AD_summary$SNP

# Get CpGs of each model
load("Models/ModelInterpretation/selCpGs.RData")

results_all <- NULL
# For each model...
for (i in 6:length(selCpGs)){
  
  # Get the included CpGs
  selCpGs1 <- selCpGs[[i]]
  
  # Get the mQTLs
  query <- list(
    cpgs = selCpGs1,
    pval = 1e-5,
    cistrans = "cis",
    clumped = 0
  )
  res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
  test_all <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
  
  # For each CpG...
  selCpGs1 <- unique(test_all$cpg)
  for (cpg in 1:length(selCpGs1)){
    
    # Get mQTL information for single CpG
    test <- as.data.frame(test_all[test_all$cpg == selCpGs1[cpg],])
    
    # Get common SNPS with GWAS data
    common_snps <- intersect(test$rsid,AD_summary$SNP)
    
    # Select common SNPs only (mQTL data)
    test <- test[test$rsid != "",]
    rownames(test) <- test$rsid
    test <- test[common_snps,]
    
    # Prepare mQTL data
    coloc_data_mQTL <- list(beta = test$beta_a1,
                            varbeta = test$se^2,
                            snp = test$rsid,
                            position = str_remove(str_remove(str_remove(test$name, ":SNP"),"chr.*:"), ":"),
                            type = "quant",
                            N = test$samplesize,
                            MAF = ifelse(test$freq_a1 < 0.5, test$freq_a1, 1-test$freq_a1))
    
    
    
    # Select common SNPs only (GWAS data)
    AD_summary_fil <- AD_summary[common_snps,]
    
    # Prepare GWAS data
    coloc_data_GWAS <- list(beta = AD_summary_fil$BETA,
                            varbeta = AD_summary_fil$SE^2,
                            snp = AD_summary_fil$SNP,
                            position = AD_summary_fil$BP,
                            type = "cc")
    
    # If # common SNPS > 0
    if (length(common_snps) > 0){
      
      # Perform coloc analysis
      my.res <- coloc.abf(dataset1=coloc_data_GWAS,
                          dataset2=coloc_data_mQTL)
      
      # Save summary statistics
      results_cpg <- as.data.frame(t(as.data.frame(my.res$summary)))
      results_cpg$CpG <- selCpGs1[cpg]
      results_cpg$Model <- names(selCpGs)[i]
      rownames(results_cpg) <- NULL
      
      # Combine with previous results
      results_all <- rbind.data.frame(results_all, results_cpg)
    }
  }
}

# Save colocalization analysis results
save(results_all, file = "Models/ModelInterpretation/Coloc/coloc_results1.RData")

# Load results (if needed)
load("Models/ModelInterpretation/Coloc/coloc_results1.RData")

# Change model names
results_all$Model[results_all$Model == "Alcohol"] <- "Alcohol consumption"
results_all$Model[results_all$Model == "Depression"] <- "Depression"
results_all$Model[results_all$Model == "Diabetes"] <- "Type II diabetes"
results_all$Model[results_all$Model == "Diet"] <- "Unhealthy diet"
results_all$Model[results_all$Model == "Education"] <- "Low education"
results_all$Model[results_all$Model == "EpiAge"] <- "Age"
results_all$Model[results_all$Model == "HDL"] <- "HDL cholesterol"
results_all$Model[results_all$Model == "Physical"] <- "Phsyical inactivity"
results_all$Model[results_all$Model == "SexMale"] <- "Sex"
results_all$Model[results_all$Model == "SysBP"] <- "Syst. blood pressure"

# Get p-value
results_all$pvalue = 1-(results_all$PP.H4.abf + results_all$PP.H3.abf)

# Calculate FDR-adjusted p-value
results_all$FDR <- p.adjust(results_all$pvalue, method = "fdr")
results_all$Sig <- ifelse(results_all$FDR < 0.05, "Yes", "No")

# Load probe annotation
load("Models/ModelInterpretation/Probe2Gene_final_ann.RData")
Probe2Gene <- Probe2Gene_final %>%
  group_by(CpG) %>%
  summarise(CpG = CpG,
            Genes = paste0(Gene))

Probe2Gene <- Probe2Gene[!duplicated(Probe2Gene),]

# Combine probe annotation with coloc results
results_all_ann <- left_join(results_all, Probe2Gene,
                              by = c("CpG" = "CpG"))
results_all_ann$pvalue = 1-(results_all_ann$PP.H4.abf + results_all_ann$PP.H3.abf)
results_all_ann$FDR <- p.adjust(results_all_ann$pvalue, method = "fdr")
results_all_ann$Sig <- ifelse(results_all_ann$FDR < 0.05, "Yes", "No")

# Make manhattan-like plot
set.seed(456)
p <- ggplot(results_all) +
  geom_point(aes(x = Model, y = -log10(pvalue+10^-10), color = Model), 
             position = position_jitter(width = 0.3, seed = 123)) +
  geom_point(aes(x = Model, y = -log10(pvalue+10^-10), color = Model, alpha = Sig), 
             position = position_jitter(width = 0.3, seed = 123), 
             shape = 1, color = "#737373", size = 2.5) +
  #geom_text_repel(data = results_all_ann[results_all_ann$FDR < 0.05,],
  #                aes(x = Model, y = 1-pvalue, 
  #                    label = paste0(CpG, " (",Genes, ")")), max.overlaps = 100,
  #                     size = 3, color = "#737373", fontface = "italic") + 
  ylab("-log10 p-value (H3 + H4)") +
  xlab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_manual(values = rep(c("#FB6A4A","#CB181D"),6)) +
  scale_alpha_manual(values = c(0,1))

# Save plot
ggsave(p, file = "Models/ModelInterpretation/Coloc/Coloc_results_v2.jpg", width = 8, height = 5)

# Get annotation of significant CpGs
results_all_ann[results_all_ann$FDR < 0.05,]
Probe2Gene_final[Probe2Gene_final$CpG == "cg19514613",]

