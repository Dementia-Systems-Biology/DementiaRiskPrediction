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

# Load DNA methylation data
load("EXTEND/Data/methSet_allNorm_fil.RData")
load("Models/ModelInterpretation/selCpGs.RData")

# Get mQTLs of model's features
query <- list(
  cpgs = selCpGs,
  pval = 1e-5,
  clumped = 1
)
res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
test <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
selSNPs <- unique(str_remove(str_remove(str_remove(test$name, "chr"), ":SNP"), ":INDEL"))

samples <- intersect(colnames(methSet_allNorm_fil), colnames(genoType))
methData <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% selCpGs,samples]
genoType_fil <- genoType[rownames(genoType) %in% selSNPs,samples]

all(colnames(genoType_fil) == colnames(methData))

jive_data <- list(Methylation = methData,
                  Genotype = genoType_fil)

set.seed(123)
jive_results <- jive(jive_data)
save(jive_results, file = "GenoMeth/jive_results_Depression.RData")

showVarExplained(jive_results)
summary(jive_results)


sum(genoType == "1/2")
sum(genoType == "2/1")
sum(genoType == "0/2")
sum(genoType == "2/0")
sum(genoType == "2/2")

load("~/Data/methSet_allNorm_fil.RData")

load("~/Data/metaData_ageFil.RData")
rownames(dat) <- dat$ID
test <- dat[str_remove(colnames(genoType), "0_"), "Basename"]
colnames(genoType) <- test


