# ============================================================================ #
# File: FUN_MachineLearning.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Machine learning functions.
# ============================================================================ #

#*****************************************************************************#
# regressionSummary
#*****************************************************************************#

regressionSummary <- function (data, lev = NULL, model = NULL){
  
  # Remove missing values from prediction
  isNA <- is.na(data$pred)
  pred <- data$pred[!isNA]
  obs <- data$obs[!isNA]

  # Calculate correlation coefficient (r)
  if (length(unique(pred)) < 2 || length(unique(obs)) < 
      2) {
    resamplCor <- NA
  }
  else {
    resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), 
                      silent = TRUE)
    if (inherits(resamplCor, "try-error")) 
      resamplCor <- NA
  }
  
  # Calculate Mean Squared Error (MSE)
  mse <- mean((pred - obs)^2)
  
  # Calculate Mean Absolute error (MAE)
  mae <- mean(abs(pred - obs))
  
  # Save MSE, RMSE, R2, and MAE as output
  out <- c(mse, sqrt(mse), resamplCor^2, mae)
  
  names(out) <- c("MSE","RMSE", "Rsquared", "MAE")
  return(out)
}


#*****************************************************************************#
# classificationSummary
#*****************************************************************************#

classificationSummary <- function (data, lev = NULL, model = NULL){
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  caret::requireNamespaceQuietStop("pROC")
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">", 
                             quiet = TRUE), silent = TRUE)
  rocAUC <- if (inherits(rocObject, "try-error")) 
    NA
  else rocObject$auc
  
  sens <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
  spec <-  specificity(data[, "pred"], data[, "obs"], lev[2])
  out <- c(rocAUC, 
           sens, 
           spec,
           sqrt(sens*spec))
  names(out) <- c("ROC", "Sens", "Spec", "GMean")
  out
}



