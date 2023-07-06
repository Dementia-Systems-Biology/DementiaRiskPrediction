library(tidyverse)
setwd("E:/Thesis/PublicationScripts/EMIF-AD")

dir <- "E:/Thesis/EMIF/EMIF_AD_imputed_GWAS_data/Predict_AD/"
files <- list.files("E:/Thesis/EMIF/EMIF_AD_imputed_GWAS_data/Predict_AD", 
                    pattern = ".profile",
                    full.names = FALSE)

for (f in files){
  profile <- read.delim(paste0(dir,f))
  profile <- profile[,c("ID2", "Profile_1")]
  if (f != files[1]){
    PGS_all <- inner_join(PGS_all, profile, by = c("ID2" = "ID2"))
  } else{
    PGS_all <- profile
  }
}

colnames(PGS_all) <- c("ID", str_remove_all(str_remove_all(files,".profile"), "EMIF_"))
save(PGS_all, file = "PGS_EMIF_AD.RData")