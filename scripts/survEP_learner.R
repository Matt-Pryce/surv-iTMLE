
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
# library(caret)
# library(xgboost)
# library(reshape2)
# library(data.table)
# library(mice)
library(Sieve)
library(survival)
library(dplyr)
library(timereg)
library(boot)
library(grf)
library(SuperLearner)
library(glmnet)
library(tidyverse)
library(riskRegression)
library(pec)
library(timeROC)
library(survSuperLearner)

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
#' @param e_SL_lib Library to be used in super learner if selected
#' @param g_covariates List containing the names of the variables to be input into the missingness model, excluding exposure
#' @param g_method Statistical technique used to run the missingness model
#' @param pse_covariates List containing the names of the variables to be input into the pseudo outcome model
#' @param pse_approach Technique used for pseudo-outcome regression
#' @param pse_method Statistical technique used to run the pseudo outcome model
#' @param pse_SL_lib Library to be used in super learner if selected for pseudo outcome model
#' @param newdata New data to create predictions for
#' 

#' @param out_SL_lib Library to be used in super learner if selected
#' @param g_SL_lib Library to be used in super learner if selected for missingness model

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
                           time_cuts = NULL,
                           splits = 10,
                           e_covariates,
                           e_method = "Parametric",
                           e_SL_lib,
                           out_covariates,
                           out_method = "Parametric",
                           out_SL_lib,
                           g_covariates,
                           g_method = "Parametric",
                           g_SL_lib,
                           iso_reg = FALSE,
                           pse_covariates,
                           pse_approach,
                           pse_method = "Parametric",
                           pse_SL_lib,
                           newdata,
                           target_option,
                           sieve_dim = NULL,
                           sieve_interaction = NULL
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
  
  analysis_data_long_all <- analysis_data_long_all %>%
    arrange(ID, tstart) %>%
    group_by(ID)
  
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

        po_o_data_long_all <- subset(analysis_data_long_all, select = c("ID","Y","A","C","tstart","time",out_covariates,"s","at_risk"))
        po_o_data_long_all <- subset(po_o_data_long_all, po_o_data_long_all$s == i)

        po_g_data_long_all <- subset(analysis_data_long_all, select = c("ID","C","tstart","time",g_covariates,"A","s"))
        po_g_data_long_all <- subset(po_g_data_long_all, po_g_data_long_all$s == i)

        po_data_long_all <- subset(analysis_data_long_all,select = c("ID","Y2","C","A","tstart","time",pse_covariates,"s","at_risk"))  
        po_data_long_all <- subset(po_data_long_all,po_data_long_all$s == i)

      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured when collecting data for nuisance models predictions in split ",i,sep=""))
        print(e)
      }
    )

    # --- Running nuisance models & obtaining predictions ---#
    #Propensity score model
    tryCatch(
      {
        PS_model <- nuis_mod_surv(model = "Propensity score",
                                  data = e_data,
                                  method = e_method,
                                  covariates = e_covariates,
                                  SL_lib = e_SL_lib,
                                  pred_data = po_e_data_long_all)

        po_data_long_all$e_pred <- PS_model$e_pred
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop(paste("An error occured in propensity score model function in split ",i,sep=""))
        print(e)
      }
    )

    # Outcome survival models
    tryCatch(
      {
        outcome_models <- nuis_mod_surv(model = "Outcome",
                                        data = o_data,
                                        method = out_method,
                                        covariates = out_covariates,
                                        pred_data_long_all = po_o_data_long_all,
                                        evt_times_uni = clean_data$evt_times_uni,
                                        SL_lib = out_SL_lib,
                                        learner = "survEP-learner")

        #Saving survival predictions
        if (out_method == "Parametric"){
          po_data_long_all$S_k_pred_0 <- outcome_models$S_k_pred_long_all_0
          po_data_long_all$S_k_pred_1 <- outcome_models$S_k_pred_long_all_1
          po_data_long_all$H_k_pred_0 <- outcome_models$H_k_pred_long_all_0
          po_data_long_all$H_k_pred_1 <- outcome_models$H_k_pred_long_all_1
          po_data_long_all$h_k_pred_0 <- outcome_models$h_k_pred_long_all_0
          po_data_long_all$h_k_pred_1 <- outcome_models$h_k_pred_long_all_1
        }
        else if (out_method == "Super learner"){
          po_data_long_all$S_k_pred_0 <- outcome_models$pred_data_long_all_pred$S_k_pred_0
          po_data_long_all$S_k_pred_1 <- outcome_models$pred_data_long_all_pred$S_k_pred_1
          po_data_long_all$H_k_pred_0 <- outcome_models$pred_data_long_all_pred$H_k_pred_0
          po_data_long_all$H_k_pred_1 <- outcome_models$pred_data_long_all_pred$H_k_pred_1
          po_data_long_all$h_k_pred_0 <- outcome_models$pred_data_long_all_pred$h_k_pred_0
          po_data_long_all$h_k_pred_1 <- outcome_models$pred_data_long_all_pred$h_k_pred_1
        }

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
                                         evt_times_uni = clean_data$evt_times_uni,
                                         SL_lib = g_SL_lib,
                                         pred_data_long_all = po_g_data_long_all)

        if (g_method == "Super learner"){
          po_data_long_all$G_k_pred <- censoring_model$G_k_pred_long_all$G_k_pred
        }
        else if (g_method == "Parametric"){
          po_data_long_all$G_k_pred <- censoring_model$G_k_pred_long_all
        }
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
  for (t in seq_along(clean_data$evt_times_uni)){

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
    if (is.null(sieve_dim) == TRUE){
      basisN <- ceiling((nrow(analysis_data))^(1/3)*ncol(X))  
    }
    else {
      basisN <- sieve_dim
    }
    basisN <- basisN + 1
    
    if (is.null(sieve_interaction) == TRUE){
      interaction_order <- 2  
    }
    else {
      interaction_order <- sieve_interaction
    }

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

    TMLE_output <- list()
    # for (t_step in 1:1){
    t_step <- 1
    TMLE_iter_complete <- 0
    while (TMLE_iter_complete == 0){
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
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])

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
          if (tp == 1){
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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
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

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }

      if (target_option == "Logistic"){   #Coeficients end up very large even when collapsing data
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])

          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          training_data_at_risk <- subset(training_data_at_risk,training_data_at_risk$h_k_pred_a != 0 & training_data_at_risk$h_k_pred_a != 1)

          #Setting up model
          mod_covs <- sieve_names
          fmla <- as.formula(paste("Y2 ~ ", paste(mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

          #Running linear regression
          fluc_mod <- glm(fmla,
                          data = training_data_at_risk,
                          offset=qlogis(training_data_at_risk$h_k_pred_a),
                          weights = training_data_at_risk$H_a,
                          family = quasibinomial())

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
          if (tp == 1){
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
        targeting_data$h_k_pred_1_star_temp <- qlogis(targeting_data$h_k_pred_1)
        targeting_data$h_k_pred_0_star_temp <- qlogis(targeting_data$h_k_pred_0)
        for (cov in 1:ncol(eps_data)){
          targeting_data$h_k_pred_1_star_temp <- targeting_data$h_k_pred_1_star_temp + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
          targeting_data$h_k_pred_0_star_temp <- targeting_data$h_k_pred_0_star_temp + targeting_data[[num_covs_pre_sieve+cov]]*
            targeting_data[[num_covs_post+cov]]
        }
        targeting_data$h_k_pred_1_star <- plogis(targeting_data$h_k_pred_1_star_temp)
        targeting_data$h_k_pred_0_star <- plogis(targeting_data$h_k_pred_0_star_temp)

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

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
      
      if (target_option == "Lasso - Linear - Option 1"){  
        #Steps:
        #   - If iteration 1, then for each time run lasso and save list of covarites with no zero coefs 
        #   - If iteration >1 then just run linear with saved covs
        #   - Run linear reg and update preds 
        
        if (t_step == 1){
          lasso_cov_list <- vector(mode = "list", length = t)
        }
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])

          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)

          #If iteration 1, run lasso model
          if (t_step == 1){
            #Setting up lasso model
            lasso_mod_covs <- sieve_names
            X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
            
            #Running lasso model
            lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                   intercept = FALSE,
                                   alpha = 1,
                                   offset = training_data_at_risk$h_k_pred_a,
                                   weights = training_data_at_risk$H_a)

            # Choose a lambda value
            lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda

            # Extract coefficients at the chosen lambda
            coef_matrix <- coef(lasso_mod, s = lambda_value)

            # Get the indices of non-zero coefficients
            nonzero_indices <- which(coef_matrix != 0)

            # Extract the names of non-zero coefficients (excluding the intercept)
            nonzero_covariates <- rownames(coef_matrix)[nonzero_indices]
            nonzero_covariates <- nonzero_covariates[nonzero_covariates != "(Intercept)"]  # Remove intercept

            #Adding to list
            lasso_cov_list[[tp]] <- nonzero_covariates
          }
          
          linear_mod_covs <- lasso_cov_list[[tp]]
          
          #Adding in check to see if no covariates were specified
          if (length(linear_mod_covs) != 0){
            #Setting up linear regression model
            fmla <- as.formula(paste("Y2 ~ ", paste(linear_mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

            #Running linear regression
            fluc_mod <- lm(fmla,
                           data = training_data_at_risk,
                           offset=training_data_at_risk$h_k_pred_a,
                           weights = training_data_at_risk$H_a)

            #Obtaining model coefficients
            eps <- coef(fluc_mod)

            #Making sure no coefs are NA
            for (i in 1:length(eps)){
              if (is.na(eps[i])==1){
                eps[i] <- 0
              }
            }
            eps_data <- data.frame(col = (NA))
            for(i in 1:length(eps)) {
              new_col <- as.numeric(eps[i])
              eps_data[ 1, i] <- new_col
              colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
            }
            
            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,linear_mod_covs]

            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]

            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)

            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)

            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)

            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)

            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)

            #Creating updated hazard estimates
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star <- training_data$h_k_pred_1_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star <- training_data$h_k_pred_0_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }

            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+1):(num_covs_post_sieve_coef_new+2))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }
  
          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        
        targeting_data <- training_data_all

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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data,
                            cov_list=lasso_cov_list)

        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }

      if (target_option == "Lasso - Linear - Option 2"){  
        #For each time run lasso and save list of covariates with no zero coefs 
        #Run linear reg and update preds 
        
        #List fro storing covs used
        cov_list <- list()
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])
          
          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          
          #Setting up lasso model
          lasso_mod_covs <- sieve_names
          X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
          
          #Running lasso model
          lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                 intercept = FALSE,
                                 alpha = 1,
                                 offset = training_data_at_risk$h_k_pred_a,
                                 weights = training_data_at_risk$H_a)
          
          # Choose a lambda value
          lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda
          
          # Extract coefficients at the chosen lambda
          coef_matrix <- coef(lasso_mod, s = lambda_value)
          
          # Get the indices of non-zero coefficients
          nonzero_indices <- which(coef_matrix != 0)
            
          # Extract the names of non-zero coefficients (excluding the intercept)
          nonzero_covariates <- rownames(coef_matrix)[nonzero_indices]
          nonzero_covariates <- nonzero_covariates[nonzero_covariates != "(Intercept)"]  # Remove intercept
          
          #Setting up linear regression model
          linear_mod_covs <- nonzero_covariates
          cov_list <- append(cov_list,linear_mod_covs)
          
          #Adding in check to see if no covariates were specified
          if (length(linear_mod_covs) != 0){
            fmla <- as.formula(paste("Y2 ~ ", paste(linear_mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term
            
            #Running linear regression
            fluc_mod <- lm(fmla,
                           data = training_data_at_risk,
                           offset=training_data_at_risk$h_k_pred_a,
                           weights = training_data_at_risk$H_a)
            
            #Obtaining model coefficients
            eps <- coef(fluc_mod)
            
            #Making sure no coefs are NA
            for (i in 1:length(eps)){
              if (is.na(eps[i])==1){
                eps[i] <- 0
              }
            }
            eps_data <- data.frame(col = (NA))
            for(i in 1:length(eps)) {
              new_col <- as.numeric(eps[i])
              eps_data[ 1, i] <- new_col
              colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
            }
            
            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,linear_mod_covs]
            
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)
            
            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)
            
            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)
            
            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)
            
            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)
            
            #Creating updated hazard estimates
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star <- training_data$h_k_pred_1_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star <- training_data$h_k_pred_0_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }
            
            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+1):(num_covs_post_sieve_coef_new+2))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }
          
          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        
        targeting_data <- training_data_all

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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data,
                            cov_list=cov_list)

        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
      
      if (target_option == "Lasso - Linear - Option 3"){  
        #For each time run lasso and save list of covariates with no zero coefs 
        #Run linear reg and update preds 
        cov_list <- list()
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])
          
          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          
          #Setting up lasso model
          lasso_mod_covs <- sieve_names
          X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
          
          #Running lasso model
          lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                 intercept = FALSE,
                                 alpha = 1,
                                 offset = training_data_at_risk$h_k_pred_a,
                                 weights = training_data_at_risk$H_a)
          
          # Choose a lambda value
          lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda
          
          # Extract coefficients at the chosen lambda
          coef_matrix <- coef(lasso_mod, s = lambda_value)
          
          # Get the indices of non-zero coefficients
          nonzero_indices <- which(coef_matrix != 0)
          nonzero_covariates_val <- coef_matrix[nonzero_indices]
          nonzero_covariates_names <- rownames(coef_matrix)[nonzero_indices]
          
          cov_list <- append(cov_list,nonzero_covariates_names)
          
          if (length(nonzero_indices) != 0){
            #Creating dataset
            eps_data <- as.data.frame(matrix(nrow = 1, ncol=length(nonzero_indices)))
            for(i in 1:length(nonzero_indices)) {
              new_col <- nonzero_covariates_val[i]
              eps_data[1,i] <- new_col
              colnames(eps_data)[i] <- paste(nonzero_covariates_names[i],"_coef",sep="")
            }
            
            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,nonzero_covariates_names]
            
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)
            
            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)

            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)

            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)
            
            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)

            #Creating updated hazard estimates
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star <- training_data$h_k_pred_1_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star <- training_data$h_k_pred_0_star + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }

            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+1):(num_covs_post_sieve_coef_new+2))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }

          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }

        targeting_data <- training_data_all

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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data,
                            cov_list=cov_list)
        
        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
      
      if (target_option == "Lasso - Logistic - Option 1"){  
        #Steps:
        #   - If iteration 1, then for each time run lasso and save list of covarites with no zero coefs 
        #   - If iteration >1 then just run logistic with saved covs
        #   - Run linear reg and update preds 
        
        if (t_step == 1){
          lasso_cov_list <- vector(mode = "list", length = t)
        }
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])
          
          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          training_data_at_risk <- subset(training_data_at_risk,training_data_at_risk$h_k_pred_a != 0 & training_data_at_risk$h_k_pred_a != 1)
          
          #If iteration 1, run lasso model
          if (t_step == 1){
            #Setting up lasso model
            lasso_mod_covs <- sieve_names
            X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
            
            #Running lasso model
            lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                   intercept = FALSE,
                                   alpha = 1,
                                   offset = qlogis(training_data_at_risk$h_k_pred_a),
                                   weights = training_data_at_risk$H_a,
                                   family = "binomial")
            
            # Choose a lambda value
            lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda

            # Extract coefficients at the chosen lambda
            coef_matrix <- coef(lasso_mod, s = lambda_value)

            # Get the indices of non-zero coefficients
            nonzero_indices <- which(coef_matrix != 0)

            # Extract the names of non-zero coefficients (excluding the intercept)
            nonzero_covariates <- rownames(coef_matrix)[nonzero_indices]
            nonzero_covariates <- nonzero_covariates[nonzero_covariates != "(Intercept)"]  # Remove intercept

            #Adding to list
            lasso_cov_list[[tp]] <- nonzero_covariates
          }
          
          logistic_mod_covs <- lasso_cov_list[[tp]]

          #Adding in check to see if no covariates were specified
          if (length(logistic_mod_covs) != 0){
            #Setting up linear regression model
            fmla <- as.formula(paste("Y2 ~ ", paste(logistic_mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

            #Running linear regression
            fluc_mod <- glm(fmla,
                            data = training_data_at_risk,
                            offset=qlogis(training_data_at_risk$h_k_pred_a),
                            weights = training_data_at_risk$H_a,
                            family = quasibinomial())

            #Obtaining model coefficients
            eps <- coef(fluc_mod)

            #Making sure no coefs are NA
            for (i in 1:length(eps)){
              if (is.na(eps[i])==1){
                eps[i] <- 0
              }
            }
            eps_data <- data.frame(col = (NA))
            for(i in 1:length(eps)) {
              new_col <- as.numeric(eps[i])
              eps_data[ 1, i] <- new_col
              colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
            }

            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,logistic_mod_covs]

            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]

            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)

            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)

            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)

            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)

            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)

            #Creating updated hazard estimates
            training_data$h_k_pred_1_star_temp <- qlogis(training_data$h_k_pred_1)
            training_data$h_k_pred_0_star_temp <- qlogis(training_data$h_k_pred_0)
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star_temp <- training_data$h_k_pred_1_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star_temp <- training_data$h_k_pred_0_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }
            training_data$h_k_pred_1_star <- plogis(training_data$h_k_pred_1_star_temp)
            training_data$h_k_pred_0_star <- plogis(training_data$h_k_pred_0_star_temp)
            
            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+3):(num_covs_post_sieve_coef_new+4))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]

            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }

          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        
        targeting_data <- training_data_all

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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data,
                            cov_list=lasso_cov_list)
        
        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
      
      if (target_option == "Lasso - Logistic - Option 2"){  
        #For each time run lasso and save list of covariates with no zero coefs 
        #Run linear reg and update preds 
        cov_list <- list()
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])
          
          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          training_data_at_risk <- subset(training_data_at_risk,training_data_at_risk$h_k_pred_a != 0 & training_data_at_risk$h_k_pred_a != 1)
          
          #Setting up lasso model
          lasso_mod_covs <- sieve_names
          X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
          
          #Running lasso model
          lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                 intercept = FALSE,
                                 alpha = 1,
                                 offset = qlogis(training_data_at_risk$h_k_pred_a),
                                 weights = training_data_at_risk$H_a,
                                 family = "binomial")
          
          # Choose a lambda value
          lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda
          
          # Extract coefficients at the chosen lambda
          coef_matrix <- coef(lasso_mod, s = lambda_value)
          
          # Get the indices of non-zero coefficients
          nonzero_indices <- which(coef_matrix != 0)
          
          # Extract the names of non-zero coefficients (excluding the intercept)
          nonzero_covariates <- rownames(coef_matrix)[nonzero_indices]
          nonzero_covariates <- nonzero_covariates[nonzero_covariates != "(Intercept)"]  # Remove intercept
          
          #Setting up linear regression model
          logistic_mod_covs <- nonzero_covariates
          cov_list <- append(cov_list,logistic_mod_covs)
          
          #Adding in check to see if no covariates were specified
          if (length(logistic_mod_covs) != 0){
            fmla <- as.formula(paste("Y2 ~ ", paste(logistic_mod_covs, collapse= "+"),"-1"))    #-1 indicators no intercept term

            #Running linear regression
            fluc_mod <- glm(fmla,
                            data = training_data_at_risk,
                            offset=qlogis(training_data_at_risk$h_k_pred_a),
                            weights = training_data_at_risk$H_a,
                            family = quasibinomial())
            
            #Obtaining model coefficients
            eps <- coef(fluc_mod)
            
            #Making sure no coefs are NA
            for (i in 1:length(eps)){
              if (is.na(eps[i])==1){
                eps[i] <- 0
              }
            }
            eps_data <- data.frame(col = (NA))
            for(i in 1:length(eps)) {
              new_col <- as.numeric(eps[i])
              eps_data[ 1, i] <- new_col
              colnames(eps_data)[i] <- paste0("V", i,"_coef_est")
            }

            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,logistic_mod_covs]
            
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)
            
            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)
            
            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)
            
            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)
            
            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)
            
            #Creating updated hazard estimates
            training_data$h_k_pred_1_star_temp <- qlogis(training_data$h_k_pred_1)
            training_data$h_k_pred_0_star_temp <- qlogis(training_data$h_k_pred_0)
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star_temp <- training_data$h_k_pred_1_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star_temp <- training_data$h_k_pred_0_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }
            training_data$h_k_pred_1_star <- plogis(training_data$h_k_pred_1_star_temp)
            training_data$h_k_pred_0_star <- plogis(training_data$h_k_pred_0_star_temp)
            
            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+3):(num_covs_post_sieve_coef_new+4))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]

            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }

          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        
        targeting_data <- training_data_all

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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data)#,
        # coefs=eps_data_all)

        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
      
      if (target_option == "Lasso - Logistic - Option 3"){  
        #For each time run lasso 
        cov_list <- list()
        for (tp in 1:t){
          #Limiting to that time point
          training_data <- subset(targeting_data,targeting_data$time == clean_data$evt_times_uni[tp])
          
          #Restricting to those at risk
          training_data_at_risk <- subset(training_data,training_data$at_risk == 1)
          training_data_at_risk <- subset(training_data_at_risk,training_data_at_risk$h_k_pred_a != 0 & training_data_at_risk$h_k_pred_a != 1)
          
          #Setting up lasso model
          lasso_mod_covs <- sieve_names
          X <- as.matrix(subset(training_data_at_risk,select = lasso_mod_covs))
          
          #Running lasso model
          lasso_mod <- cv.glmnet(x = X,y = training_data_at_risk$Y2,
                                 intercept = FALSE,
                                 alpha = 1,
                                 offset = qlogis(training_data_at_risk$h_k_pred_a),
                                 weights = training_data_at_risk$H_a,
                                 family = "binomial")
          
          # Choose a lambda value
          lambda_value <- lasso_mod$lambda.min  #Using the minimum lambda
          
          # Extract coefficients at the chosen lambda
          coef_matrix <- coef(lasso_mod, s = lambda_value)
          
          # Get the indices of non-zero coefficients
          nonzero_indices <- which(coef_matrix != 0)
          nonzero_covariates_val <- coef_matrix[nonzero_indices]
          nonzero_covariates_names <- rownames(coef_matrix)[nonzero_indices]
          
          cov_list <- append(cov_list,nonzero_covariates_names)
          
          if (length(nonzero_indices) != 0){
            #Creating dataset
            eps_data <- as.data.frame(matrix(nrow = 1, ncol=length(nonzero_indices)))
            for(i in 1:length(nonzero_indices)) {
              new_col <- nonzero_covariates_val[i]
              eps_data[1,i] <- new_col
              colnames(eps_data)[i] <- paste(nonzero_covariates_names[i],"_coef",sep="")
            }
            
            #Collecting sieve basis to keep    (The next few lines remove surplus covs and allows for us to update hazard easier)
            sieve_basis_keep <- training_data[,nonzero_covariates_names]
            
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Count how many covs now
            num_covs_pre_sieve_new <- ncol(training_data)
            
            #Adding back in only ones which were used
            training_data <- cbind(training_data,sieve_basis_keep)
            
            #Counting how many covs at this point
            num_covs_post_sieve_new <- ncol(training_data)
            
            #Adding in coeficient estimates
            training_data <- merge(training_data,eps_data,all.x = T)
            
            #Counting how many covs at this point
            num_covs_post_sieve_coef_new <- ncol(training_data)
            
            #Creating updated hazard estimates
            training_data$h_k_pred_1_star_temp <- qlogis(training_data$h_k_pred_1)
            training_data$h_k_pred_0_star_temp <- qlogis(training_data$h_k_pred_0)
            for (cov in 1:ncol(eps_data)){
              training_data$h_k_pred_1_star_temp <- training_data$h_k_pred_1_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
              training_data$h_k_pred_0_star_temp <- training_data$h_k_pred_0_star_temp + training_data[[num_covs_pre_sieve_new+cov]]*
                training_data[[num_covs_post_sieve_new+cov]]
            }
            training_data$h_k_pred_1_star <- plogis(training_data$h_k_pred_1_star_temp)
            training_data$h_k_pred_0_star <- plogis(training_data$h_k_pred_0_star_temp)
            
            #Remove sieve covs/coefs
            keep_cov_list <- c(1:num_covs_pre_sieve_new,(num_covs_post_sieve_coef_new+3):(num_covs_post_sieve_coef_new+4))
            training_data <- training_data[,keep_cov_list]
          }
          else {
            #Removing sieve basis covs
            covs_not_sieve <- c(1:num_covs_pre_sieve,(num_covs_post_sieve+1):dim(training_data)[2])
            training_data <- training_data[,covs_not_sieve]
            
            #Do not update hazards
            training_data$h_k_pred_1_star <- training_data$h_k_pred_1
            training_data$h_k_pred_0_star <- training_data$h_k_pred_0
          }
          
          #Recollecting full data
          if (tp == 1){
            training_data_all <- training_data
            # eps_data_all <- eps_data
          }
          else{
            training_data_all <- rbind(training_data_all,training_data)
            # eps_data_all <- rbind(eps_data_all,eps_data)
          }
        }
        
        targeting_data <- training_data_all
        
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
        targeting_data <- targeting_data %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_1_star = cumsum(h_k_pred_1_star))
        targeting_data <- targeting_data %>% group_by(ID) %>% mutate(H_k_pred_0_star = cumsum(h_k_pred_0_star))

        #Creating updated survival probabilities
        targeting_data$S_k_pred_1_star <- exp(-targeting_data$H_k_pred_1_star)
        targeting_data$S_k_pred_0_star <- exp(-targeting_data$H_k_pred_0_star)

        targeting_data_prev <- targeting_data

        #--- Saving information from iteration ---#
        iter_output <- list(data=targeting_data,
                            cov_list=cov_list)
        
        TMLE_output <- append(TMLE_output,list(iter_output))

        if (TMLE_iter_complete == 1){
          #Collecting final tp survival estimates
          final_iter_data <- targeting_data
          max_tps <- max(unique(targeting_data$time))
          final_tp_ests <- subset(final_iter_data,final_iter_data$time==max_tps)
          TMLE_output <- append(TMLE_output,list(final_tp_ests))
        }
      }
    }
    #--- Saving iteration outputs from each TMLE process ---#
    TMLE_output_all <- append(TMLE_output_all,list(TMLE_output))
  }
  #--- Pulling together final dataset for predictions ---#
  num_tps <- length(clean_data$evt_times_uni)
  for (i in 1:num_tps){
    iters <- length(TMLE_output_all[[i]])
    if (i == 1){
      pred_dataset_all <- TMLE_output_all[[1]][[iters]]
    }
    else {
      pred_dataset <- TMLE_output_all[[i]][[iters]]
      pred_dataset_all <- rbind(pred_dataset_all,pred_dataset)
    }
  }


  #---------------------------#
  #--- Isotonic regression ---#
  #---------------------------#
  pred_dataset_all <- pred_dataset_all %>%
    arrange(ID, tstart) %>%
    group_by(ID)

  if (iso_reg == TRUE){
    #Number individuals
    pred_dataset_all <- pred_dataset_all %>%
      group_by(ID) %>%
      mutate(person_number = cur_group_id())

    #Iterate through these
    data_split <- split(pred_dataset_all, pred_dataset_all$person_number)

    # Apply isotonic regression to each person's data using lapply
    results <- lapply(data_split, function(pred_dataset_all_temp) {

      # Defining outcome lists and time
      y0 <- -pred_dataset_all_temp$S_k_pred_0_star   #Negative as only does monotonic increasing
      y1 <- -pred_dataset_all_temp$S_k_pred_1_star
      time <- pred_dataset_all_temp$time

      # Running isotonic regression
      iso_mod0 <- isoreg(time, y0)
      iso_mod1 <- isoreg(time, y1)

      # Assign the isotonic regression results
      pred_dataset_all_temp$S_k_pred_0_star_iso <- -iso_mod0$yf
      pred_dataset_all_temp$S_k_pred_1_star_iso <- -iso_mod1$yf

      return(pred_dataset_all_temp)
    })

    pred_dataset_all <- do.call(rbind, results)

    pred_dataset_all$S_k_pred_0_star_final <- pmin(pred_dataset_all$S_k_pred_0_star_iso,1)
    pred_dataset_all$S_k_pred_1_star_final <- pmin(pred_dataset_all$S_k_pred_1_star_iso,1)
  }
  else if (iso_reg == FALSE){
    pred_dataset_all$S_k_pred_0_star_final <- pmin(pred_dataset_all$S_k_pred_0_star,1)
    pred_dataset_all$S_k_pred_1_star_final <- pmin(pred_dataset_all$S_k_pred_1_star,1)
  }

  #---------------------------------#
  #--- Pseudo-outcome regression ---#
  #---------------------------------#
  #--- Pooled logistic regression ---#
  if (pse_approach == "Pooled logistic"){
    pred_dataset_all$pse_Y <- pred_dataset_all$S_k_pred_1_star_final - pred_dataset_all$S_k_pred_0_star_final
    pse_model <- nuis_mod_surv(model = "Pseudo outcome - Pooled logistic",
                                     data = pred_dataset_all,
                                     method = pse_method,
                                     covariates = pse_covariates,
                                     SL_lib = pse_SL_lib,
                                     pred_data_long_all = clean_data$newdata_long_all)
  }
  else if (pse_approach == "None"){
    pred_dataset_all$pse_Y <- pred_dataset_all$S_k_pred_1_star_final - pred_dataset_all$S_k_pred_0_star_final
  }

  #--------------#
  #--- Output ---#
  #--------------#
  if (pse_approach == "Pooled logistic"){
    output <- list(TMLE_output=list(TMLE_output_all),
                   data_with_preds=pred_dataset_all,
                   predictions=pse_model$pred)
  }
  else {
    output <- list(TMLE_output=list(TMLE_output_all),
                   data_with_preds=pred_dataset_all)
  }
  return(output)
}




#---------------#
#--- Example ---#
#---------------#
load("~/PhD/DR_Missing_Paper/Data_example/Data/ACTG175_data.RData")

#Defining censoring indicator
ACTG175_data$censor_ind <- 1 - ACTG175_data$cens

event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph"), function(alg) {
  c(alg,"All")
})

# event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph", "survSL.rfsrc"), function(alg) {
#   c(alg, "survscreen.glmnet","All")
# })

start_time <- proc.time()

survEP_check <- survEP_learner(data = ACTG175_data,
                               id = "pidnum",
                               time = "days",
                               outcome = "cens",
                               censor = "censor_ind",
                               exposure = "treat",
                               time_cuts = seq(from=200,to=1100,by=100),
                               splits = 1,
                               e_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               e_method = "Parametric",
                               e_SL_lib = c("SL.glm"),
                               out_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               out_method = "Parametric",
                               out_SL_lib = event.SL.library,
                               g_covariates = c("age","wtkg","hemo","homo","drugs","karnof"),
                               g_method = "Parametric",
                               g_SL_lib = event.SL.library,
                               iso_reg = FALSE,
                               pse_covariates = c("age"),#,"wtkg","hemo","homo","drugs","karnof"),
                               pse_approach = "None",
                               pse_method = "Parametric",
                               pse_SL_lib = c("SL.lm"),
                               newdata = ACTG175_data,
                               target_option = "Lasso - Linear - Option 3")#,
                               # sieve_dim = 15,
                               # sieve_interaction = 2)


end_time <- proc.time()
end_time - start_time


#Key end options
# - No isotonic reg and no final reg
# - No isotonic reg and pooled logsitic reg
# - Isotonic reg and no final reg
# - Isotonic reg and pooledl logistic reg 






#--- Notes ---#
#Time to event = "days"
#Outcome indicator = "cens"
#Censoring indicator = "censor_ind"

# Code for running targeting step when estimating P(T^1>t|Z) - P(T^0>t|Z) for one t 
# Needs to be run for each t in (1,...,t_max)                

##################################################

#Logistic is failing when run over all time points (works for only first one)
#Even for first one coeficients never converge and instead estimates go to -inf hence claims convergence 

#Censoring SL cannot currently do full time points as the process requires saving a matrix that is too big 
#Need to update process 




