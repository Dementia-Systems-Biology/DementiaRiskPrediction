#******************************************************************************#
# Load packages and functions
#******************************************************************************#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory to the shiny directory
setwd("Shiny")

# These packages are required. Install if necessary.
library(caret)
library(glmnet)
library(ranger)
library(data.table)

# We also need source the MMRSCalculation_FUN.R file for some custom functions
source("MMRSCalculation_FUN.R")


#******************************************************************************#
# Download models and additional files from Zenodo
#******************************************************************************#

# Check if the required files are in the Shiny directory and download the
# missing files
all_files <- list.files()
files <- c("finalModels.RData", 
           "all_cpgs.RData", 
           "example.csv",
           "Fit_EMIF_MCI_RF.RData")

for (f in files){
  if (!(f %in% all_files)){
    download.file(url = paste0("https://zenodo.org/record/8306113/files/",f), 
                  f, mode = "wb")
  }
}

if (!("Predictors_Shiny_by_Groups.csv" %in% all_files)){
  download.file(url = "https://zenodo.org/record/4646300/files/Predictors_Shiny_by_Groups.csv", 
                "Predictors_Shiny_by_Groups.csv", mode = "wb")
}


#******************************************************************************#
# Load models from Shiny directory
#******************************************************************************#

# Load EXTEND models
load("finalModels.RData")

# Load model's CpGs
load("all_cpgs.RData")

# Load Marioni models
cpgs <- read.csv("Predictors_Shiny_by_Groups.csv", header = T) 

# Load epi-MCI model
load("Fit_EMIF_MCI_RF.RData")


#******************************************************************************#
# Read data
#******************************************************************************#

# Read data (change data file name if necssary)
data_all <- fread("example.csv")

# Format data (set first column as rownames)
data_all <- as.data.frame(data_all)
rownames(data_all) <- data_all[,1]
data_all <- data_all[,-1]

#******************************************************************************#
# Calculate the Methylation Profile Scores (MPSs)
#******************************************************************************#

# MPSs from the EXTEND models
extend_out <- extendModel(data = data_all, 
                          finalModels = finalModels,
                          all_cpgs =  all_cpgs)

# MPSs from Marioni models
marioni_out <- marioniModel(data = data_all, cpgs = cpgs)


# Combine MPSs from EXTEND and Marioni et al. in a single dataframe
predictedScore_factors <- cbind.data.frame(extend_out,
                                           marioni_out[,c("EpiAge",
                                                          "Alcohol",
                                                          "BMI", 
                                                          "HDL",
                                                          "Smoking")])
rownames(predictedScore_factors) <- rownames(extend_out)

# Download the MPSs
write.table(predictedScore_factors,file = "predictedScore_MPSs.csv",row.names = TRUE,sep = ",",quote = FALSE)

#******************************************************************************#
# Calculate the Multivariate Methylation Risk Scores (MMRSs)
#******************************************************************************#

# Predict
pred_MMRS <- predict(fit, predictedScore_factors, type = "prob")

# Combine into data frame
predictDF <- data.frame(SampleID = rownames(predictedScore_factors), 
                        MMRS = log(pred_MMRS$MCI/(1-pred_MMRS$MCI)))

# Download the MMRSs
write.table(predictDF, file = "predictedScore_MMRS.csv", row.names = FALSE, sep = ",",quote = FALSE)


# END