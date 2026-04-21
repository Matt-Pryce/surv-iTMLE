#######################################################################################
# Script: Causal survival forest
# Date: 09/10/24
# Author: Matt Pryce 
# Notes: Code to run nuisance models/CSF 
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


CSF <- function(data,
                id,
                time, 
                outcome,
                censor,
                exposure,
                splits = 10,
                prop_est = TRUE,
                covariates,
                e_method = "Parametric",
                time_cuts = NULL,
                # out_covariates,
                # out_method = "Parametric",
                # g_covariates,
                # g_method = "Parametric",
                # weight_method = "Parametric",
                # pse_approach,
                # pse_method = "Parametric",
                # pse_SL_lib,
                # nuisance_estimates_input = 0,
                # o_0_pred = NA,
                # o_1_pred = NA,
                # g_pred = NA,
                newdata
){
  
  #-----------------------#
  #--- Data management ---#
  #-----------------------#
  
  clean_data <- data_manage_surv(data = data,
                                 learner = "CSF",
                                 id = id,
                                 time=time, 
                                 outcome=outcome,
                                 censor=censor,
                                 exposure=exposure,
                                 time_cuts=time_cuts,
                                 splits = splits,
                                 e_covariates = covariates,
                                 # out_covariates = out_covariates,
                                 # g_covariates = g_covariates,
                                 pse_covariates = covariates,
                                 newdata = newdata)
  
  analysis_data <- clean_data$data

  if (prop_est == TRUE){
    #--------------------------------------#
    #--- Running propensity score model ---#
    #--------------------------------------#

    #--- Iterating over each split (cross-fitting) ---#
    for (i in 0:(splits-1)){

      #--- Collecting data for training ---#
      tryCatch(
        {
          #Data for propensity score
          e_data <- analysis_data

          #Sorting by split specs
          if (splits == 1){
            e_data <- subset(e_data,e_data$s == i)
          }
          else if (splits == 10){
            e_data <- subset(e_data,e_data$s != i)
          }
        },
        #if an error occurs, tell me the error
        error=function(e) {
          stop(paste("An error occured when collecting data for training nuisance models in split ",i,sep=""))
          print(e)
        }
      )

      #--- Collecting data to obtain predictions ---#
      tryCatch(
        {
          #Creating datasets for obtaining model predictions
          po_e_data<- subset(analysis_data, select = c(covariates,"s"))
          po_e_data <- subset(po_e_data, po_e_data$s == i)

          po_data <- analysis_data
          po_data <- subset(po_data,po_data$s == i)

        },
        #if an error occurs, tell me the error
        error=function(e) {
          stop(paste("An error occured when collecting data for nuisance models predictions in split ",i,sep=""))
          print(e)
        }
      )

      #--- Running nuisance models & obtaining predictions ---#
      #Propensity score model
      tryCatch(
        {
          PS_model <- nuis_mod_surv(model = "Propensity score",
                                    data = e_data,
                                    method = e_method,
                                    covariates = covariates,
                                    pred_data = po_e_data)

          po_data$e_pred <- predict(PS_model$e_mod,newdata=po_e_data,type="response")
        },
        #if an error occurs, tell me the error
        error=function(e) {
          stop(paste("An error occured in propensity score model function in split ",i,sep=""))
          print(e)
        }
      )


      #--- Collecting datasets ---#
      if (i==0){
        po_data_complete <- po_data
      }
      else {
        po_data_complete <- rbind(po_data_complete,po_data)
      }
    }


    #----------------------------------------#
    #--- Running CSF and collating output ---#
    #----------------------------------------#
    #Defining covariate set
    X <- as.matrix(subset(po_data_complete,select = covariates))
    
    #Identifying time points to get estimates for 
    time_cuts <- clean_data$evt_times_uni
  
    #---Looping over each time ---#
    for (tp in seq_along(time_cuts)){
      max_time <- time_cuts[tp]
      
      #--- Running CSF ---#
      mod <- causal_survival_forest(X = X,   #Covariates
                                    Y = po_data_complete$time,   #Event time
                                    W = po_data_complete$A,   #Exposure
                                    D = po_data_complete$Y,   #Event indicator
                                    W.hat = po_data_complete$e_pred,   #Prop score estimates
                                    target = c("survival.probability"),
                                    horizon = max_time)#,        #Probabilty estimated at this
                                    # failure.times = time_cuts)  #Failure times

      #--- Collating output ---#
      new_data <- clean_data$newdata
      X.test <- as.matrix(subset(clean_data$newdata,select=covariates))
      pred <- predict(mod, X.test,estimate.variance=TRUE)
      new_data$pred <- pred$predictions
      new_data$pred_var <- pred$variance.estimates
      new_data$pred_LCI <- new_data$pred - 1.96*sqrt(new_data$pred_var)
      new_data$pred_UCI <- new_data$pred + 1.96*sqrt(new_data$pred_var)
      
      t_start <- clean_data$evt_times_uni0[tp]
      t_end <- clean_data$evt_times_uni0[tp+1]
      new_data$t_start <- rep(t_start,dim(new_data)[1])
      new_data$t_end <- rep(t_end,dim(new_data)[1])
      
      if (tp == 1){
        new_data_all <- new_data
      }
      else {
        new_data_all <- rbind(new_data_all,new_data)
      }
    }
  }
  else if (prop_est == FALSE){
    #----------------------------------------#
    #--- Running CSF and collating output ---#
    #----------------------------------------#
    #Defining covariate set
    X <- as.matrix(subset(analysis_data,select = covariates))
    
    #Identifying time points to get estimates for 
    time_cuts <- clean_data$evt_times_uni
    
    #---Looping over each time ---#
    for (tp in seq_along(time_cuts)){
      max_time <- time_cuts[tp]
      
      #--- Running CSF ---#
      mod <- causal_survival_forest(X = X,   #Covariates
                                    Y = analysis_data$time,   #Event time
                                    W = analysis_data$A,   #Exposure
                                    D = analysis_data$Y,   #Event indicator
                                    target = c("survival.probability"),
                                    horizon = max_time)        #Probabilty estimated at this
      
      #--- Collating output ---#
      new_data <- clean_data$newdata
      X.test <- as.matrix(subset(clean_data$newdata,select=covariates))
      pred <- predict(mod, X.test,estimate.variance=TRUE)
      new_data$pred <- pred$predictions
      new_data$pred_var <- pred$variance.estimates
      new_data$pred_LCI <- new_data$pred - 1.96*sqrt(new_data$pred_var)
      new_data$pred_UCI <- new_data$pred + 1.96*sqrt(new_data$pred_var)
      
      t_start <- clean_data$evt_times_uni0[tp]
      t_end <- clean_data$evt_times_uni0[tp+1]
      new_data$t_start <- rep(t_start,dim(new_data)[1])
      new_data$t_end <- rep(t_end,dim(new_data)[1])
      
      if (tp == 1){
        new_data_all <- new_data
      }
      else {
        new_data_all <- rbind(new_data_all,new_data)
      }
    }
  }
  
  #-----------------------------#
  #--- Returning information ---#
  #-----------------------------#
  output <- list(CSF_mod=mod,
                 preds=new_data_all)
  return(output)
}



##############################################################

# 
# #---------------#
# #--- Example ---#
# #---------------#
# load("ACTG175_data.RData")  #Update
# 
# #Defining censoring indicator
# ACTG175_data$censor_ind <- 1 - ACTG175_data$cens
# 
# 
# start_time <- proc.time()
# 
# CSF_check <- CSF(data = ACTG175_data,
#                  id = "pidnum",
#                  time = "days", 
#                  outcome = "cens",
#                  censor = "censor_ind",
#                  exposure = "treat",
#                  time_cuts = seq(from=200,to=1100,by=100),
#                  splits = 10,
#                  prop_est = TRUE,
#                  covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
#                  e_method = "Parametric",
#                  newdata = ACTG175_data)
# 
# 
# end_time <- proc.time()
# end_time - start_time


#Add max horizon in 

