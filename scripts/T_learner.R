#######################################################################################
# Script: T-Learner function
# Date: 20/03/24
# Author: Matt Pryce 
# Notes: T-learner function for time to event data
#######################################################################################

library(survival)
library(dplyr)
library(timereg)
library(boot)
library(grf)
library(SuperLearner)


#######################################
#--- For single time point setting ---#
#######################################

#UPDATE

#' @param data The data frame containing all required information
#' @param id Identification for individuals
#' @param outcome The name of the outcome of interest
#' @param exposure The name of the exposure of interest
#' @param outcome_observed_indicator Indicator identifying when the outcome variable is observed (1=observed, 0=missing)
#' @param out_method Statistical technique used to run the outcome models
#' @param out_covariates  List containing the names of the variables to be input into each outcome model
#' @param out_SL_lib Library to be used in super learner if selected
#' @param out_SL_strat Indicator for if stratification should be used in the super learner CV (stratifies outcomes - Only use if outcome binary)
#' @param g_method Statistical technique used to run the missingness model
#' @param g_covariates List containing the names of the variables to be input into the missingness model, excluding exposure
#' @param g_SL_lib Library to be used in super learner if selected for missingness model
#' @param imp_covariates Covariates to be used in SL imputation model if SL imputation used
#' @param imp_SL_lib SL libaray for imputation model if SL imputation used
#' @param imp_SL_strat Whether SL CV folds are stratified in the imputation model if SL imputation used  
#' @param nuisance_estimates_input Indicator for whether nuisance estimates provided
#' @param o_0_pred Variable name for unexposed outcome predictions (if provided)
#' @param o_1_pred Variable name for exposed outcome predictions (if provided)
#' @param g_pred Variable name for censoring predictions (if provided)
#' @param newdata New data to create predictions for

#' @return A list containing: CATE estimates, a dataset used to train the learner, outcome models (if run) 


T_learner <- function(data,
                      id,
                      time, 
                      outcome,
                      censor,
                      exposure,
                      truncation = NULL,
                      time_cuts = "N/A",
                      out_method = c("Parametric","Super learner"),
                      out_covariates,
                      out_SL_lib,
                      nuisance_estimates_input = 0,
                      o_0_pred = NA,
                      o_1_pred = NA,
                      newdata
){
  
  #-----------------------#
  #--- Data management ---#
  #-----------------------#
  
  clean_data <- data_manage_surv(data = data,
                                 learner = "T-learner",
                                 id = id,
                                 time=time, 
                                 outcome=outcome,
                                 censor=censor,
                                 exposure=exposure,
                                 truncation=truncation,
                                 time_cuts=time_cuts,
                                 splits = 1,
                                 out_covariates = out_covariates,
                                 newdata = newdata)
  
  analysis_data <- clean_data$data

  #------------------------------#
  #--- Running outcome models ---#
  #------------------------------#

  if (nuisance_estimates_input == 0){
    tryCatch(
      {
        o_data <- subset(analysis_data,select=-s)

        outcome_models <- nuis_mod_surv(model = "Outcome",
                                        data = o_data,
                                        method = out_method,
                                        covariates = out_covariates,
                                        pred_data_long_all = clean_data$newdata_long_all,
                                        evt_times_uni = clean_data$evt_times_uni,
                                        SL_lib = out_SL_lib,
                                        learner = "T-learner",
                                        LT = clean_data$LT_data)

        #Saving survival predictions
        pred_data <- clean_data$newdata_long_all
        if (out_method == "Parametric"){
          pred_data$S_k_pred_0 <- outcome_models$S_k_pred_long_all_0
          pred_data$S_k_pred_1 <- outcome_models$S_k_pred_long_all_1
        }
        else if (out_method == "Super learner"){
          pred_data$S_k_pred_0 <- outcome_models$pred_data_long_all_pred$S_k_pred_0
          pred_data$S_k_pred_1 <- outcome_models$pred_data_long_all_pred$S_k_pred_1
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured in outcome model function",sep=""))
        print(e)
      }
    )
  }

  #-------------------------------#
  #--- Creating CATE estimates ---#
  #-------------------------------#

  if (nuisance_estimates_input == 0){
    pred_data$CATE_est <- pred_data$S_k_pred_1 - pred_data$S_k_pred_0
  }
  else if (nuisance_estimates_input == 1){
    pred_data <- clean_data$newdata_long_all
    pred_data$CATE_est <- clean_data$newdata_long_all[["o_1_pred"]] - clean_data$newdata_long_all[["o_0_pred"]]
  }

  #-----------------------------#
  #--- Returning information ---#
  #-----------------------------#
  if (nuisance_estimates_input == 0){
    output <- list(CATE_est=pred_data$CATE_est,
                   analysis_data = analysis_data,
                   outcome_models=outcome_models,
                   pred_data=pred_data)
  }
  else {
    output <- list(CATE_est=pred_data$CATE_est,
                   pred_data=pred_data)
  }
  return(output)
}



##############################################################


#---------------#
#--- Example ---#
#---------------#
load("~/PhD/DR_Missing_Paper/Data_example/Data/ACTG175_data.RData")

ACTG175_data$censor_ind <- 1 - ACTG175_data$cens
ACTG175_data <- ACTG175_data[1:800,]
ACTG175_data$trunc <- round(runif(800,0,250))
ACTG175_data <- subset(ACTG175_data,ACTG175_data$days > ACTG175_data$trunc)


LT <- 0

if (LT==1){
  event.SL.library <- c("SL.mean",
                        "SL.glm")
}
if (LT == 0){
  event.SL.library <- cens.SL.library <- lapply(c("survSL.km","survSL.expreg"), function(alg) {
    c(alg,"All")
  })
}

start_time <- proc.time()

T_check <- T_learner(data = ACTG175_data,
                     id = "pidnum",
                     time = "days", 
                     outcome = "cens",
                     censor = "censor_ind",
                     exposure = "treat",
                     # truncation = "trunc",
                     time_cuts = seq(from=200,to=1100,by=100),
                     out_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                     out_method = "Super learner",
                     out_SL_lib = event.SL.library,
                     newdata = ACTG175_data)


end_time <- proc.time()
end_time - start_time










