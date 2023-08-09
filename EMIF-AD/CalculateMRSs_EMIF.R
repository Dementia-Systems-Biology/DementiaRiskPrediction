# ============================================================================ #
# File: CalculateMRSs_EMIF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Calculate the MRSs in the EMIF-AD cohort
# ============================================================================ #

# Load packages
library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/metaData_EMIF.RData")  # Meta data
load("EMIF-AD/Data/X_EMIF_imp.RData")     # Imputed methylation data


#################################################################################

# EXTEND Models

################################################################################

# Load EXTEND models
load("Models/MRS_Models/finalModels.RData")

# Predict factors
factors <- names(finalModels)
predictedScore_factors <- matrix(NA, nrow = nrow(X_EMIF_imp), ncol = length(factors))
colnames(predictedScore_factors) <- factors
rownames(predictedScore_factors) <- rownames(X_EMIF_imp)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- finalModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-ncol(model$trainingData)]
  
  # Make predictions
  if (f <= 2){
    predictedScore_factors[,f] <- predict(model, X_EMIF_imp[,features])
  }
  if (f > 2){
    prob <- predict(model, X_EMIF_imp[,features], type = "prob")$Yes
    predictedScore_factors[,f] <- log(prob/(1-prob))
  }
  
}
predictedScore_factors <- as.data.frame(predictedScore_factors)


#################################################################################

# Hillary Models

################################################################################

# Read models
# Models are available at: https://zenodo.org/record/4646300#.ZFJg53ZBxPY
cpgs <- read.csv("Models/MRS_Models/Predictors_Shiny_by_Groups.csv", header = T) 

# Load data
data <- X_EMIF_imp


# Check if Data needs to be Transposed
if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data <- t(data) 
}

# Subset CpG sites to those present on list for predictors 
coef <- data[intersect(rownames(data), cpgs$CpG_Site),]

# Identify CpGs missing from input dataframe, include them and provide values 
# as mean methylation value at that site

coef <- if(nrow(coef) == 3629) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat)
} 


# Convert NAs to Mean Value for all individuals across each probe 
na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))

# Conversion to Beta Values if M Values are likely present  
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}


# Calculating the predictors 
loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = colSums(tmp2)
  }
} 
out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
out$ID <- row.names(out) 
out <- out[,c(ncol(out),1:(ncol(out)-1))] 
colnames(out) <- c("ID", "EpiAge", "Alcohol", "BMI", "BodyFat", "HDL", "Smoking", "WHR")


predictedScore_factors <- cbind.data.frame(predictedScore_factors,
                                           out[,c("EpiAge",
                                                  "Alcohol",
                                                  "BMI", 
                                                  "HDL",
                                                  "Smoking")])

# Save MRSs
save(predictedScore_factors, file = "EMIF-AD/Data/predictedScore_factors_EMIF.RData")

