#######################################################################################
# Script: Data management function for survival data
# Date: 19/08/24
# Author: Matt Pryce 
# Notes: - Used when running any models for difference in survival probabilties
#        - Do not allow for nuisance estimates to be externally estimated
#######################################################################################

#library(caTools)
#library(tidyverse)
#library(ggplot2)
#library(DAAG)
#library(glmnet)
#library(randomForest)
#library(caret)
#library(grf)
#library(xgboost)
#library(reshape2)
#library(data.table)
#library(SuperLearner)
#library(mice)

#' @param data The data frame containing all required information
#' @param leanrer Which learner is data management being run on 
#' @param analysis Type of analysis to be run (If not mDR, complete case or outcome imputation) 
#' @param id Identification for individuals
#' @param outcome The name of the outcome of interest
#' @param exposure The name of the exposure of interest
#' @param time_cuts A vector of timepoints to cut at
#' @param splits The number of splits used in nuisance training 
#' @param outcome_observed_indicator Indicator identifying when the outcome variable is observed (1=observed, 0=missing)
#' @param out_covariates  List containing the names of the variables to be input into each outcome model
#' @param e_covariates  List containing the names of the variables to be input into the propensity score model
#' @param g_covariates List containing the names of the variables to be input into the missingness model, excluding exposure
#' @param imp_covariates Covariates to be used in SL imputation model if SL imputation used
#' @param pse_covariates List containing the names of the variables to be input into the pseudo outcome model
#' @param newdata New data to create predictions for

#' @return Cleaned data for training and testing, along with indicators for whether the outcome is binary or continuous 


data_manage_surv <- function(data,
                             learner,
                             analysis_type = "N/A",
                             id,
                             time, 
                             outcome,
                             censor,
                             exposure,
                             time_cuts,
                             splits,
                             out_covariates,
                             e_covariates,
                             g_covariates,
                             pse_covariates,
                             newdata
){

  #---------------------------------------------#
  #--- Checking and keeping appropriate data ---#
  #---------------------------------------------#
  tryCatch( 
    {
      vars <- c(id, time, outcome, censor, exposure)
      
      #For all analysis choices and both learners
      if (learner != "CSF"){
        vars <- append(vars,out_covariates)  
      }
      if (learner == "survEP-learner" | learner == "M-learner"){     #Add others if functions when coding other functions
        vars <- append(vars,e_covariates)
        vars <- append(vars,g_covariates)
        vars <- append(vars,pse_covariates)
      }
      if (learner == "CSF"){
        vars <- append(vars,e_covariates)
        vars <- append(vars,pse_covariates)
      }

      vars <- vars[!duplicated(vars)]

      data <- subset(data,select = vars)
    },
    error=function(e) {
      stop('An error occured when selecting the appropriate variables')
      print(e)
    }
  )
  
  #------------------------#
  #--- Rename variables ---#
  #------------------------#
  tryCatch(
    {
      names(data)[names(data) == id] <- "ID"
      names(data)[names(data) == time] <- "time"
      names(data)[names(data) == exposure] <- "A"
      names(data)[names(data) == outcome] <- "Y"
      names(data)[names(data) == censor] <- "C"
    },
    error=function(e) {
      stop('An error occured when renaming variables')
      print(e)
    }
  )

  tryCatch(
    {
      #Checking number of splits is correct
      if (splits != 1 & splits != 10){
        stop("Number of splits not compatible")
      }

      #Creating splits
      data$s <- rep(1:length(data$Y),1) %% splits
    },
    #if an error occurs, tell me the error
    error=function(e) {
      stop('An error occured when splitting the data')
      print(e)
    }
  )


  #-------------------------------#
  #--- Creating long form data ---#
  #-------------------------------#

  #--- Identifying event times ---#
  evts <- subset(data,data$Y == 1)                 #Keeping people who had an event
  evt_times <- evts$time                           #Collecting event dates
  evt_times_uni <- sort(unique(evt_times))         #Collecting unique dates
  evt_times_uni0 <- c(0,evt_times_uni)             #Event times including 0

  #--- Creating data with one observation per person per time ---# (used in the targeting step)

  if (is.null(time_cuts) == 1){
    time_cuts <- evt_times_uni
    time_cuts0 <- c(0,evt_times_uni)             #Event times including 0
  }
  else {
    time_cuts0 <- c(0,time_cuts)
  }

  #Pre-allocate a list to store dataframes
  data_long_list <- vector("list", length(time_cuts))

  #Loop through evt_times_uni and store results in the list
  for (i in seq_along(time_cuts)) {
    data_long_temp <- data
    data_long_temp$time2 <- data_long_temp$time
    data_long_temp$tstart <- time_cuts0[i]
    data_long_temp$time <- time_cuts0[i + 1]

    # Store the result in the list
    data_long_list[[i]] <- data_long_temp
  }

  # Combine all dataframes into one
  data_long_all <- do.call(rbind, data_long_list)

  if (is.null(time_cuts) == 1){
    #Defining at risk indicators for this dataset
    data_long_all <- data_long_all %>% mutate(at_risk = case_when(
      time2 >= time ~ 1,
      time2 <  time ~ 0
    ))

    data_long_all <- data_long_all %>% mutate(Y2 = case_when(
      time2==time & Y == 1 ~ 1,
      time2==time & Y == 0 ~ 0,
      time2>time          ~ 0,
      time2<time          ~ as.numeric(NA)
    ))
  }
  else {
    #Defining indicator for outcome occurring in that interval
    data_long_all <- data_long_all %>% mutate(Y2 = case_when(
      time2<=time & time2>tstart & Y == 1 ~ 1,
      time2<=time & time2>tstart & Y == 0 ~ 0,
      time2>time          ~ 0,
      time2<time          ~ as.numeric(NA)
    ))
    #Defining at risk indicators for this dataset
    data_long_all <- data_long_all %>% mutate(at_risk = case_when(
      is.na(Y2) == 1 ~ 0,
      is.na(Y2) == 0 ~ 1
    ))
  }


  #-----------------------#
  #--- Format new data ---#
  #-----------------------#
  tryCatch(
    {
      #NEW DATA MUST BE INPUT IN SHORT FORMAT
      #--- Formating newdata variable names ---#
      names(newdata)[names(newdata) == id] <- "ID"
      names(newdata)[names(newdata) == outcome] <- "Y"
      names(newdata)[names(newdata) == time] <- "time"
      names(newdata)[names(newdata) == exposure] <- "A"

      #--- Creating dataset with one observation per person at each time ---#
      if (is.null(time_cuts) == 1){
        time_cuts <- evt_times_uni
        time_cuts0 <- c(0,evt_times_uni)             #Event times including 0
      }
      else {
        time_cuts0 <- c(0,time_cuts)
      }

      #Pre-allocate a list to store dataframes
      newdata_long_list <- vector("list", length(time_cuts))

      #Loop through evt_times_uni and store results in the list
      for (i in seq_along(time_cuts)) {
        newdata_long_temp <- newdata
        newdata_long_temp$time2 <- newdata_long_temp$time
        newdata_long_temp$tstart <- time_cuts0[i]
        newdata_long_temp$time <- time_cuts0[i + 1]

        # Store the result in the list
        newdata_long_list[[i]] <- newdata_long_temp
      }

      # Combine all dataframes into one
      newdata_long_all <- do.call(rbind, newdata_long_list)

      if (is.null(time_cuts) == 1){
        #Defining at risk indicators for this dataset
        data_long_all <- data_long_all %>% mutate(at_risk = case_when(
          time2 >= time ~ 1,
          time2 <  time ~ 0
        ))

        data_long_all <- data_long_all %>% mutate(Y2 = case_when(
          time2==time & Y == 1 ~ 1,
          time2==time & Y == 0 ~ 0,
          time2>time          ~ 0,
          time2<time          ~ as.numeric(NA)
        ))
      }
      else {
        #Defining indicator for outcome occurring in that interval
        data_long_all <- data_long_all %>% mutate(Y2 = case_when(
          time2<=time & time2>tstart & Y == 1 ~ 1,
          time2<=time & time2>tstart & Y == 0 ~ 0,
          time2>time          ~ 0,
          time2<time          ~ as.numeric(NA)
        ))
        #Defining at risk indicators for this dataset
        data_long_all <- data_long_all %>% mutate(at_risk = case_when(
          is.na(Y2) == 1 ~ 0,
          is.na(Y2) == 0 ~ 1
        ))
      }

      if (learner == "survEP-learner" | learner == "M-learner"){
        new_data_vars <- c("ID",pse_covariates,"time")
        newdata_long_all <- subset(newdata_long_all,select=new_data_vars)
      }
      if (learner == "T-learner"){
        new_data_vars <- c("ID",out_covariates,"tstart","time","A","Y")
        newdata_long_all <- subset(newdata_long_all,select=new_data_vars)
      }
    },
    error=function(e) {
      stop('An error occured when formatting the new data')
      print(e)
    }
  )


  #--------------------#
  #--- Logic checks ---#
  #--------------------#
  tryCatch(
    {
      #Checking if exposure is binary
      if (as.numeric(all(data$A %in% 0:1)) == 0){
        stop("Exposure must be binary")
      }
      #Checking outcome is binary
      if (as.numeric(all(data$Y %in% 0:1)) == 0){
        stop("Outcome must be binary")
      }
      #Checking censoring indicator is binary
      if (as.numeric(all(data$C %in% 0:1)) == 0){
        stop("Censoring indicator must be binary")
      }
    },
    error=function(e) {
      stop('The exposure, outcome indicator or censoring indicator are non-binary')
      print(e)
    }
  )

  #-----------------------------#
  #--- Returning information ---#
  #-----------------------------#
  if (learner == "survEP-learner" | learner == "M-learner" | learner == "T-learner" | learner == "CSF"){
    if (time_cuts[1] == "N/A"){
      output <- list(data=data,
                     data_long_all=data_long_all,
                     evt_times_uni=evt_times_uni,
                     evt_times_uni0=evt_times_uni0,
                     newdata=newdata,
                     newdata_long_all=newdata_long_all)
    }
    else {
      output <- list(data=data,
                     data_long_all=data_long_all,
                     evt_times_uni=time_cuts,
                     evt_times_uni0=time_cuts0,
                     newdata=newdata,
                     newdata_long_all=newdata_long_all)
    }
  }
  return(output)
}


