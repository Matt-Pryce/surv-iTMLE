
#######################################################################################
# Script: survEP-Learner function
# Date: 25/03/24
# Author: Matt Pryce 
# Notes: survEP-learner, a function for estimating the conditional difference in survival
#        probabilities between exposed and unexposed groups. 
#######################################################################################

# library(caTools)
# library(tidyverse)
# library(ggplot2)
# library(DAAG)
# library(glmnet)
# library(randomForest)
# library(caret)
# library(grf)
# library(xgboost)
# library(reshape2)
# library(data.table)
# library(SuperLearner)
# library(mice)
library(Sieve)
library(survival)
library(dplyr)
library(timereg)
library(boot)

#######################################
#--- For single time point setting ---#
#######################################

#' @param data The data frame containing all required information
#' @param id Identification for individuals
#' @param time Time of either event for censoring
#' @param outcome Name of the indicator for if time represent an outcome (1=Event, 0=Censored) 
#' @param censor Name of the indicator for if time represents censoring (1=Censored, 0=Event)
#' @param exposure The name of the baseline exposure of interest
#' @param time_cuts The time splits used to predict nuisance predictions for and use in targeting step
#' @param splits Number of splits for cross-fitting (Variations allowed: 10)
#' @param out_covariates  List containing the names of the variables to be input into each outcome model
#' @param out_method Statistical technique used to run the outcome models
#' @param e_covariates  List containing the names of the variables to be input into the propensity score model
#' @param e_method Statistical technique used to run the propensity score model  (Currently only parametric)
#' @param g_covariates List containing the names of the variables to be input into the missingness model, excluding exposure
#' @param g_method Statistical technique used to run the missingness model
#' @param pse_covariates List containing the names of the variables to be input into the pseudo outcome model
#' @param newdata New data to create predictions for
#' 

#' @param out_SL_lib Library to be used in super learner if selected
#' @param e_SL_lib Library to be used in super learner if selected
#' @param g_SL_lib Library to be used in super learner if selected for missingness model
#' @param nuisance_estimates_input Indicator for whether nuisance estimates provided
#' @param o_0_pred Variable name for unexposed outcome predictions (if provided)
#' @param o_1_pred Variable name for exposed outcome predictions (if provided)
#' @param e_pred Variable name for propensity score predictions (if provided)
#' @param g_pred Variable name for censoring predictions (if provided)
#' @param pse_method Statistical technique used to run the pseudo outcome model
#' @param pse_SL_lib Library to be used in super learner if selected for pseudo outcome model
#' 
#' @param sieve_num_basis Dimension of sieve basis
#' @param sieve_interaction_order Interaction order for sieves
#' @param imp_covariates Covariates to be used in SL imputation model if SL imputation used
#' @param imp_SL_lib SL libaray for imputation model if SL imputation used


#' @return A list containing: CATE estimates, a dataset used to train the learner, dataset containing all 
#'         pseudo-outcome predictions (if splits 4) 


survEP_learner <- function(data,
                           id,
                           time, 
                           outcome,
                           censor,
                           exposure,
                           time_cuts = "N/A",
                           splits = 10,
                           e_covariates,
                           e_method = "Parametric",
                           # e_SL_lib,
                           out_covariates,
                           out_method = "Parametric",
                           g_covariates,
                           g_method = "Parametric",
                           pse_covariates,
                           newdata,
                           target_option
){
  #-----------------------#
  #--- Data management ---#
  #-----------------------#

  clean_data <- data_manage_surv(data = data,
                                 learner = "survEP-learner",
                                 id = id,
                                 time=time, 
                                 outcome=outcome,
                                 censor=censor,
                                 exposure=exposure,
                                 time_cuts=time_cuts,
                                 splits = splits,
                                 e_covariates = e_covariates,
                                 out_covariates = out_covariates,
                                 g_covariates = g_covariates,
                                 pse_covariates = pse_covariates,
                                 newdata = newdata)
  
  analysis_data <- clean_data$data
  analysis_data_long_all <- clean_data$data_long_all
  
  
  #-------------------------------#
  #--- Running nuisance models ---#
  #-------------------------------#

  #--- Iterating over each split (cross-fitting) ---#
  for (i in 0:(splits-1)){

    #--- Collecting data for training models ---#
    tryCatch(
      {
        #Data for propensity score, outcome models & censoring models
        e_data <- analysis_data
        o_data <- analysis_data
        g_data <- analysis_data

        #Sorting by split specs
        if (splits == 1){
          e_data <- subset(e_data,e_data$s == i)
          o_data <- subset(o_data,o_data$s == i)
          g_data <- subset(g_data,g_data$s == i)
        }
        else if (splits == 10){
          e_data <- subset(e_data,e_data$s != i)
          o_data <- subset(o_data,o_data$s != i)
          g_data <- subset(g_data,g_data$s != i)
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured when collecting data for training nuisance models in split ",i,sep=""))
        print(e)
      }
    )

    #--- Collecting data to obtain nuisance model predictions for ---#
    tryCatch(
      {
        #Creating datasets for obtaining model predictions
        po_e_data_long_all <- subset(analysis_data_long_all, select = c(e_covariates,"s"))
        po_e_data_long_all <- subset(po_e_data_long_all, po_e_data_long_all$s == i)

        po_o_data_long_all <- subset(analysis_data_long_all, select = c("ID","Y","A","C","tstart","time",out_covariates,"s"))
        po_o_data_long_all <- subset(po_o_data_long_all, po_o_data_long_all$s == i)

        po_g_data_long_all <- subset(analysis_data_long_all, select = c("ID","C","tstart","time",g_covariates,"A","s"))
        po_g_data_long_all <- subset(po_g_data_long_all, po_g_data_long_all$s == i)

        po_data_long_all <- subset(analysis_data_long_all,select = c("ID","Y2","C","A","tstart","time",pse_covariates,"s","at_risk"))   #Figure out what is needed here
        po_data_long_all <- subset(po_data_long_all,po_data_long_all$s == i)

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
                                  covariates = e_covariates,
                                  pred_data = po_e_data_long_all)

        po_data_long_all$e_pred <- predict(PS_model$e_mod,newdata=po_e_data_long_all,type="response")
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured in propensity score model function in split ",i,sep=""))
        print(e)
      }
    )

    #Outcome survival models
    tryCatch(
      {
        outcome_models <- nuis_mod_surv(model = "Outcome",
                                        data = o_data,
                                        method = out_method,
                                        covariates = out_covariates,
                                        pred_data_long_all = po_o_data_long_all,
                                        evt_times_uni = clean_data$evt_times_uni)

        #Saving survival predictions
        po_data_long_all$S_k_pred_0 <- outcome_models$S_k_pred_long_all_0
        po_data_long_all$S_k_pred_1 <- outcome_models$S_k_pred_long_all_1
        po_data_long_all$H_k_pred_0 <- outcome_models$H_k_pred_long_all_0
        po_data_long_all$H_k_pred_1 <- outcome_models$H_k_pred_long_all_1
        po_data_long_all$h_k_pred_0 <- outcome_models$h_k_pred_long_all_0
        po_data_long_all$h_k_pred_1 <- outcome_models$h_k_pred_long_all_1

        #Defining S(k|A=a) and lambda(k|A=a) - All time point data
        po_data_long_all <- po_data_long_all %>% mutate(S_k_pred_a = case_when(
          A == 1 ~ S_k_pred_1,
          A == 0 ~ S_k_pred_0
        ))
        po_data_long_all <- po_data_long_all %>% mutate(H_k_pred_a = case_when(
          A == 1 ~ H_k_pred_1,
          A == 0 ~ H_k_pred_0
        ))
        po_data_long_all <- po_data_long_all %>% mutate(h_k_pred_a = case_when(
          A == 1 ~ h_k_pred_1,
          A == 0 ~ h_k_pred_0
        ))
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured in outcome model function in split ",i,sep=""))
        print(e)
      }
    )

    #Censoring function
    tryCatch(
      {
        censoring_model <- nuis_mod_surv(model = "Censoring",
                                         data = g_data,
                                         method = g_method,
                                         covariates = g_covariates,
                                         pred_data_long_all = po_g_data_long_all)

        po_data_long_all$G_k_pred <- censoring_model$G_k_pred_long_all
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured in censoring model function in split ",i,sep=""))
        print(e)
      }
    )

    #--- Collecting nuisance models ---#
    #Add later if needed, may be easier to make all the predictions up front and store in a list


    #--- Collecting long datasets ---#
    if (i==0){
      po_data_long_all_complete <- po_data_long_all

      Surv_t_pred_long_all_list_complete <- vector("list", length(clean_data$evt_times_uni))
      for (j in 1:length(clean_data$evt_times_uni)){
        #Loading estimates
        est_all <- outcome_models$Surv_t_pred_long_all_list[[j]]
        Surv_t_pred_long_all_list_complete[[j]] <- est_all
      }
    }
    else {
      po_data_long_all_complete <- rbind(po_data_long_all_complete,po_data_long_all)

      for (j in 1:length(clean_data$evt_times_uni)){
        #Merging to previous estimates
        est_all <- outcome_models$Surv_t_pred_long_all_list[[j]]
        Surv_t_pred_long_all_list_complete[[j]] <- rbind(Surv_t_pred_long_all_list_complete[[j]],est_all)
      }
    }
  }


  #------------------------#
  #--- Updating hazards ---#
  #------------------------#
  #Implementing targeting step for each t in {1,...,t_max}
  #Loops over them (only one time for now)
  TMLE_output_all <- list()
  for (t in seq_along(clean_data$evt_times_uni)[1]){

    #-----------------------#
    #--- Data management ---#
    #-----------------------#
    
    #Trimming data to correct length
    po_data_long_complete_all_TMLE <- subset(po_data_long_all_complete,po_data_long_all_complete$time <= clean_data$evt_times_uni[t])
    po_data_t_preds_complete_all_TMLE <- Surv_t_pred_long_all_list_complete[[t]]
    po_data_t_preds_complete_all_TMLE <- subset(po_data_t_preds_complete_all_TMLE,po_data_t_preds_complete_all_TMLE$time <= clean_data$evt_times_uni[t])

    #Adding in S_t_pred variables to data
    po_data_long_complete_all_TMLE$S_t_pred_0 <- po_data_t_preds_complete_all_TMLE$S_t_pred_long_all_0
    po_data_long_complete_all_TMLE$S_t_pred_1 <- po_data_t_preds_complete_all_TMLE$S_t_pred_long_all_1
    po_data_long_complete_all_TMLE$S_t_pred_a <- po_data_t_preds_complete_all_TMLE$S_t_pred_long_all_a


    #----------------------------#
    #--- Creating sieve basis ---#
    #----------------------------#

    #--- All time point data ---#
    sieve_data <- po_data_long_complete_all_TMLE

    #Defining a subset of covariates that we want to explore heterogeneity over
    pse_covariates <- pse_covariates

    X <- as.matrix(subset(sieve_data, select=pse_covariates))
    # basisN <- 3   #Arbitrary selection so we have a small number for now (this will need to be based on sample size/no. covs - See below)
    basisN <- ceiling((nrow(analysis_data))^(1/3)*ncol(X))
    basisN <- basisN + 1
    interaction_order <- 2

    #Defining sieve basis
    sieve_basis <- Sieve::sieve_preprocess(as.matrix(X),
                                           basisN = basisN,
                                           interaction_order = interaction_order,
                                           type = "cosine")$Phi
    sieve_basis <- as.data.frame(sieve_basis)

    #Attaching to data
    num_covs_pre_sieve <- ncol(po_data_long_complete_all_TMLE)
    po_data_long_complete_all_TMLE <- cbind(po_data_long_complete_all_TMLE,sieve_basis)
    num_covs_post_sieve <- ncol(po_data_long_complete_all_TMLE)

    #----------------------#
    #--- Targeting step ---#
    #----------------------#

    #EDIT TO BE A WHILE LOOP

    TMLE_output <- list()
    for (t_step in 1:1){
    # t_step <- 1
    # TMLE_iter_complete <- 0
    # while (TMLE_iter_complete == 0){

      if (t_step == 1){
        targeting_data <- po_data_long_complete_all_TMLE
      }
      else {
        #--- Pulling data from end of previous iteration ---#
        targeting_data <- targeting_data_prev

        #--- Calculating S_k_pred_a_star and h_k_pred_a_star ---#
        targeting_data <- targeting_data %>% mutate(S_k_pred_a_star = case_when(
          A == 1 ~ S_k_pred_1_star,
          A == 0 ~ S_k_pred_0_star
        ))

        targeting_data <- targeting_data %>% mutate(h_k_pred_a_star = case_when(
          A == 1 ~ h_k_pred_1_star,
          A == 0 ~ h_k_pred_0_star
        ))

        #--- Calculating S_t_pred_0_star, S_t_pred_1_star and S_t_pred_a_star ---#
        #Restricting data to time point of interest
        targeting_data_t <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[t])

        #Trimming variables
        targeting_data_t <- subset(targeting_data_t,select=c(ID,S_k_pred_0_star,S_k_pred_1_star,S_k_pred_a_star))
        names(targeting_data_t) <- c("ID","S_t_pred_0_star","S_t_pred_1_star","S_t_pred_a_star")

        #Merging back with all data to give t preds
        targeting_data <- merge(targeting_data,targeting_data_t,by="ID",all.x = T)

        #--- Dropping additional variables which are to be redefined ---#
        keep_covs <- names(targeting_data)[1:(num_covs_pre_sieve)]
        keep_covs <- append(keep_covs,c("S_t_pred_0_star","S_t_pred_1_star","S_t_pred_a_star",
                                        "S_k_pred_0_star","S_k_pred_1_star","S_k_pred_a_star",
                                        "h_k_pred_0_star","h_k_pred_1_star","h_k_pred_a_star"))
        targeting_data <- subset(targeting_data,select=keep_covs)

        targeting_data <- subset(targeting_data,select=-c(S_t_pred_0,S_t_pred_1,S_t_pred_a,
                                                          S_k_pred_0,S_k_pred_1,S_k_pred_a,
                                                          h_k_pred_0,h_k_pred_1,h_k_pred_a))

        #--- Redefining updated survival predictions as current survival prediction ---#
        names(targeting_data)[names(targeting_data) == "S_t_pred_0_star"] <- "S_t_pred_0"
        names(targeting_data)[names(targeting_data) == "S_t_pred_1_star"] <- "S_t_pred_1"
        names(targeting_data)[names(targeting_data) == "S_t_pred_a_star"] <- "S_t_pred_a"
        names(targeting_data)[names(targeting_data) == "S_k_pred_0_star"] <- "S_k_pred_0"
        names(targeting_data)[names(targeting_data) == "S_k_pred_1_star"] <- "S_k_pred_1"
        names(targeting_data)[names(targeting_data) == "S_k_pred_a_star"] <- "S_k_pred_a"
        names(targeting_data)[names(targeting_data) == "h_k_pred_0_star"] <- "h_k_pred_0"
        names(targeting_data)[names(targeting_data) == "h_k_pred_1_star"] <- "h_k_pred_1"
        names(targeting_data)[names(targeting_data) == "h_k_pred_a_star"] <- "h_k_pred_a"

        targeting_data <- cbind(targeting_data,sieve_basis)
      }

      #--- Create clever covariate for each time point ---#
      targeting_data$H_0 <- targeting_data$at_risk *
        ((1-targeting_data$A)/(1-targeting_data$e_pred)) *
        (1/targeting_data$G_k_pred) *
        (targeting_data$S_t_pred_0/targeting_data$S_k_pred_0)

      targeting_data$H_1 <- targeting_data$at_risk *
        (targeting_data$A/targeting_data$e_pred) *
        (1/targeting_data$G_k_pred) *
        (targeting_data$S_t_pred_1/targeting_data$S_k_pred_1)

      targeting_data$H_a <- targeting_data$at_risk *
        ((targeting_data$A/targeting_data$e_pred) + ((1-targeting_data$A)/(1-targeting_data$e_pred))) *
        (1/targeting_data$G_k_pred) *
        (targeting_data$S_t_pred_a/targeting_data$S_k_pred_a)

      #Sieve basis covariate names
      num_covs_pre <- ncol(targeting_data)  #For creating list of cov names
      sieve_names <- names(sieve_basis)
      num_covs_post <- ncol(targeting_data) #For creating list of cov names


      #--- Running targeting model/s ---#
      if (target_option == "Additive"){
        mod_covs <- sieve_names
        fmla <- as.formula(paste("Surv(tstart,time,Y2==1) ~ ", paste(mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

        #Restrict training data to those at risk  (if move to all data)
        training_data <- subset(targeting_data,targeting_data$at_risk == 1)

        #Additive hazard model
        add_mod <- aalen(fmla,
                         training_data,
                         n.sim=100,
                         max.time=clean_data$evt_times_uni[t],
                         offsets=training_data$h_k_pred_a,
                         weights = training_data$H_a,
                         id=training_data$ID)

        #Collecting coefficients from model for each basis variable at each time
        cumulative_coefficients <- as.data.frame(add_mod$cum)
        for (var in 1:length(mod_covs)){
          #Creating non cumulative coefficient estimates
          cumulative_coefficients$new_var <- diff(c(0,cumulative_coefficients[,mod_covs[var]]))
          #Renaming
          names(cumulative_coefficients)[names(cumulative_coefficients) == "new_var"] <- paste(mod_covs[var],"_coef_est",sep="")
        }
        names(cumulative_coefficients)[names(cumulative_coefficients) == "time"] <- "tstart"
        keep_vars <- append("tstart",paste(mod_covs,"_coef_est",sep=""))
        cumulative_coefficients <- subset(cumulative_coefficients,select=keep_vars)

        #--- Creating updated hazard estimates for each person at each time ---#
        #Adding basis coefficients to data
        targeting_data <- merge(targeting_data,cumulative_coefficients,by="tstart",all.x = T)

        #Creating updated hazard estimates
        targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1
        targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0
        for (cov in 1:(ncol(cumulative_coefficients)-1)){
          targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1_star + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
          targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0_star + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
        }
      }

      if (target_option == "Linear"){
        for (log_tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[log_tp])

          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)

          #Setting up model
          mod_covs <- sieve_names
          fmla <- as.formula(paste("Y2 ~ ", paste(mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

          #Running linear regression
          fluc_mod <- lm(fmla,
                         data = training_data_at_risk,
                         offset=training_data_at_risk$h_k_pred_a,
                         weights = training_data_at_risk$H_a)

          #Obtaining model coefficients
          eps <- coef(fluc_mod)
          eps_data <- data.frame(col = (NA))
          for(i in 1:length(eps)) {
            new_col <- as.numeric(eps[i])
            eps_data[ 1, i] <- new_col
            colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
          }

          #Merging with data that has all individuals in
          training_data <- merge(training_data,eps_data,all.x = T)

          #Recollecting full data
          if (log_tp == 1){
            training_data_all <- training_data
            eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        targeting_data <- training_data_all

        #Creating updated hazard estimates
        targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1
        targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0
        for (cov in 1:ncol(eps_data)){
          targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1_star + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
          targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0_star + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
        }
        
        #--- Checks to see if model has converged ---#
        if (t_step < 3){
          TMLE_iter_complete <- 0
        }
        else {
          diff_1 <- abs(targeting_data$h_k_pred_1_star-targeting_data$h_k_pred_1)
          diff_0 <- abs(targeting_data$h_k_pred_0_star-targeting_data$h_k_pred_0)
          if (sum(diff_1)< 1.00e-05 & sum(diff_0)< 1.00e-05){
            TMLE_iter_complete <- 1
          }
          else{
            TMLE_iter_complete <- 0
          }
        }
        
        if (t_step > 30){
          return("MSE did no converge within 30 iterations")
        }
        
        t_step <- t_step + 1
        
        #--- Create updated survival probabilities ---#
        #Creating updated cumulative hazards at each time
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))
        
        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)
        
        targeting_data_prev <- targeting_data
        
        #--- Saving information from iteration ---#   
        iter_output <- list(data=targeting_data,
                            coefs=eps_data_all)
        
        TMLE_output <- append(TMLE_output,list(iter_output))
      
        

        
        #Pick final iteration for each time point 
        #Look at data
        #Only want the end time point from targeting dataset
        #Do this for each time point then merge these to gain one full dataset 
        
      }



      if (target_option == "Logistic"){   #Coeficients end up very large even when collapsing data
        for (log_tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[log_tp])

          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)

          #Setting up model
          mod_covs <- sieve_names
          fmla <- as.formula(paste("Y2 ~ ", paste(mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

          #Running logistic regression
          fluc_mod <- glm(fmla,
                          data = training_data_at_risk,
                          offset=training_data_at_risk$h_k_pred_a,
                          weights = training_data_at_risk$H_a,
                          family = quasibinomial())

          # #Obtaining model coefficients
          # eps <- coef(fluc_mod)
          # eps_data <- data.frame(col = (NA))
          # for(i in 1:length(eps)) {
          #   new_col <- as.numeric(eps[i])
          #   eps_data[ 1, i] <- new_col
          #   colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
          # }

          # #Merging with data that has all individuals in
          # training_data <- merge(training_data,eps_data,all.x = T)
          # 
          # #Creating full term to add and taking inverse logit
          # targeting_data$add_term <- 0
          # for (cov in 1:ncol(eps_data)){
          #   targeting_data$add_term <- targeting_data$add_term + targeting_data[[num_covs_pre_sieve+cov]]*
          #     targeting_data[[num_covs_post+cov]]
          # }
          # 
          # #Recollecting full data
          # if (log_tp == 1){
          #   training_data_all <- training_data
          #   eps_data_all <- eps_data
          # }
          # else{
          #   training_data_all <- rbind(training_data_all,training_data)
          #   eps_data_all <- rbind(eps_data_all,eps_data)
          # }
        }
        # targeting_data <- training_data_all
      # UPDATE
      #   # #Creating updated hazard estimates
      #   # targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1
      #   # targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0
      #   # for (cov in 1:ncol(eps_data)){
      #   #   targeting_data$h_k_pred_1_star <- targeting_data$h_k_pred_1_star + targeting_data[[num_covs_pre_sieve+cov]]*
      #   #     targeting_data[[num_covs_post+cov]]
      #   #   targeting_data$h_k_pred_0_star <- targeting_data$h_k_pred_0_star + targeting_data[[num_covs_pre_sieve+cov]]*
      #   #     targeting_data[[num_covs_post+cov]]
      #   # }
      }

      # #--- Create updated survival probabilities ---#
      # #Creating updated cumulative hazards at each time
      # targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
      # targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))
      # 
      # #Creating updated survival probabilities
      # targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
      # targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)
      # 
      # targeting_data_prev <- targeting_data


      # #--- Saving information from iteration ---#
      # if (target_option == "Additive"){
      #   iter_output <- list(data=targeting_data,
      #                       coefs=cumulative_coefficients)
      # }
      # else if (target_option == "Logistic" | target_option == "Linear"){
      #   iter_output <- list(data=targeting_data,
      #                       coefs=eps_data_all)
      # }
      # 
      # TMLE_output <- append(TMLE_output,list(iter_output))
      # 
      # #--- Checks to see if model has converged ---#
      # if (t_step < 3){
      #   TMLE_iter_complete <- 0
      # }
      # else {
      #   diff_1 <- abs(targeting_data$S_k_pred_1_star-targeting_data$S_k_pred_1)
      #   diff_0 <- abs(targeting_data$S_k_pred_0_star-targeting_data$S_k_pred_0)
      #   if (sum(diff_1)< 1.00e-05 & sum(diff_0)< 1.00e-05){
      #     TMLE_iter_complete <- 1
      #   }
      #   else{
      #     TMLE_iter_complete <- 0
      #   }
      # }

      

    }
    
    # #Collecting final survival probabilities at that time 
    # final_tp <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[log_tp])
    # 
    # if (t == 1){
    #   final_results <- final_tp
    # }
    # else {
    #   final_results <- rbind(final_results,final_tp)
    # }
    # 
    # #--- Saving iteration outputs from each TMLE process ---#
    # TMLE_output_all <- append(TMLE_output_all,list(TMLE_output))
  }


  # return(list(list(TMLE_output_all),final_results))
  return(training_data_at_risk)
}




#---------------#
#--- Example ---#
#---------------#
load("~/PhD/DR_Missing_Paper/Data_example/Data/ACTG175_data.RData")

#Defining censoring indicator
ACTG175_data$censor_ind <- 1 - ACTG175_data$cens


start_time <- proc.time()

survEP_check <- survEP_learner(data = ACTG175_data,
                               id = "pidnum",
                               time = "days", 
                               outcome = "cens",
                               censor = "censor_ind",
                               exposure = "treat",
                               time_cuts = seq(from=200,to=1200,by=100),
                               splits = 10,
                               e_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               e_method = "Parametric",
                               out_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               out_method = "Parametric",
                               g_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               g_method = "Parametric",
                               pse_covariates = c("age"),
                               newdata = ACTG175_data,
                               target_option = "Logistic")


end_time <- proc.time()
end_time - start_time



#--- Notes ---#
#Time to event = "days"
#Outcome indicator = "cens"
#Censoring indicator = "censor_ind"

# Code for running targeting step when estimating P(T^1>t|Z) - P(T^0>t|Z) for one t 
# Needs to be run for each t in (1,...,t_max)                

##################################################

