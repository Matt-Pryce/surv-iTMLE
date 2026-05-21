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

#' @description Prepares survival data for use in the surv-iTMLE framework. Selects and renames
#'   variables to a standardised internal naming scheme, assigns cross-fitting splits, converts
#'   wide-format data to the long interval format required by the nuisance models, and applies
#'   logic checks on the input variables.
#'
#' @param data A data frame containing all variables required for the analysis.
#' @param learner The learner calling this function. One of \code{"surviTMLE-learner"},
#'   \code{"M-learner"}, \code{"T-learner"}, \code{"CSF"}.
#' @param analysis_type Type of analysis to be run. Default \code{"N/A"}.
#' @param id Name of the column containing individual identifiers.
#' @param time Name of the column containing the observed event or censoring time.
#' @param outcome Name of the binary event indicator (1 = event, 0 = censored).
#' @param censor Name of the binary censoring indicator (1 = censored, 0 = event).
#' @param exposure Name of the binary baseline exposure variable.
#' @param truncation Name of the column containing left-truncation times, or \code{NULL} if
#'   there is no left truncation.
#' @param time_cuts A numeric vector of time points at which the long-format intervals are cut,
#'   or \code{"N/A"} to use all unique observed event times.
#' @param splits Number of cross-fitting splits. Valid values: \code{1} or \code{10}.
#' @param out_covariates Character vector of covariate names for the outcome model.
#' @param e_covariates Character vector of covariate names for the propensity score model.
#' @param g_covariates Character vector of covariate names for the censoring model, excluding
#'   the exposure variable.
#' @param h_covariates Character vector of covariate names for the left-truncation weight model,
#'   or \code{NULL} if not required.
#' @param pse_covariates Character vector of covariate names for the pseudo-outcome regression model.
#' @param newdata A data frame on which final predictions are made. Typically the same as
#'   \code{data}.
#'
#' @return A list containing the processed training dataset, long-format data, unique event times,
#'   formatted prediction data, and a left-truncation indicator.


data_manage_surv <- function(data,
                             learner,
                             analysis_type = "N/A",
                             id,
                             time, 
                             outcome,
                             censor,
                             exposure,
                             truncation = NULL,
                             time_cuts,
                             splits,
                             out_covariates,
                             e_covariates,
                             g_covariates,
                             h_covariates,
                             pse_covariates,
                             newdata
){

  #------------------------------------#
  #--- Sanitize variable names ---#
  #------------------------------------#
  tryCatch({
    # Replace spaces and invalid characters with . in all column names
    names(data)    <- make.names(names(data))
    names(newdata) <- make.names(names(newdata))

    # Update parameter strings to match sanitized column names
    id       <- make.names(id)
    time     <- make.names(time)
    outcome  <- make.names(outcome)
    censor   <- make.names(censor)
    exposure <- make.names(exposure)
    if (!is.null(truncation)) truncation <- make.names(truncation)

    # Update covariate vectors
    out_covariates <- make.names(out_covariates)
    e_covariates   <- make.names(e_covariates)
    g_covariates   <- make.names(g_covariates)
    pse_covariates <- make.names(pse_covariates)
    if (!is.null(h_covariates)) h_covariates <- make.names(h_covariates)
  },
  error = function(e) {
    stop('An error occurred when sanitizing variable names')
    print(e)
  })

  #---------------------------------------------#
  #--- Checking and keeping appropriate data ---#
  #---------------------------------------------#
  tryCatch(
    {
      #Checking if data is left truncated
      if(is.null(truncation) == 0){
        LT_data <- 1
      }
      else {
        LT_data <- 0
      }
      
      vars <- c(id, time, outcome, censor, exposure)
      if (LT_data == 1){
        vars <- append(vars,truncation)
      }
      
      #For all analysis choices and both learners
      if (learner != "CSF"){
        vars <- append(vars,out_covariates)
      }
      if (learner == "surviTMLE-learner" | learner == "M-learner"){     #Add others if functions when coding other functions
        vars <- append(vars,e_covariates)
        vars <- append(vars,g_covariates)
        vars <- append(vars,pse_covariates)
      }
      if (learner == "CSF"){
        vars <- append(vars,e_covariates)
        vars <- append(vars,pse_covariates)
      }
      if (learner == "surviTMLE-learner" & LT_data == 1){
        vars <- append(vars,h_covariates)
      }

      vars <- vars[!duplicated(vars)]

      data <- subset(data,select = vars)
    },
    error=function(e) {
      stop('An error occurred when selecting the appropriate variables')
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
      if (LT_data == 1){
        names(data)[names(data) == truncation] <- "Q"
      }
    },
    error=function(e) {
      stop('An error occurred when renaming variables')
      print(e)
    }
  )

  #------------------------------------------------------------#
  #--- Check pse_covariates are numeric; create indicators ---#
  #------------------------------------------------------------#
  tryCatch({
    non_numeric_pse <- pse_covariates[
      !sapply(data[pse_covariates], is.numeric)
    ]

    if (length(non_numeric_pse) > 0) {
      message(
        "The following pse_covariates are non-numeric and will be expanded ",
        "to binary indicators: ",
        paste(non_numeric_pse, collapse = ", ")
      )

      # Keep numeric pse_covariates unchanged
      new_pse_covariates <- pse_covariates[
        sapply(data[pse_covariates], is.numeric)
      ]

      for (var in non_numeric_pse) {
        # Union of levels across data and newdata to keep indicators consistent
        lvls <- sort(unique(c(as.character(data[[var]]),
                              as.character(newdata[[var]]))))
        # Drop reference level (first alphabetically)
        for (lvl in lvls[-1]) {
          col_name <- paste0(var, ".", make.names(lvl))
          data[[col_name]]    <- as.numeric(data[[var]]    == lvl)
          newdata[[col_name]] <- as.numeric(newdata[[var]] == lvl)
          new_pse_covariates  <- c(new_pse_covariates, col_name)
        }
      }

      pse_covariates <- new_pse_covariates
    }
  },
  error = function(e) {
    stop('An error occurred when checking/converting non-numeric pse_covariates')
    print(e)
  })

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
      stop('An error occurred when splitting the data')
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
    if (LT_data == 1){
      data_long_all <- data_long_all %>% mutate(at_risk = case_when(
        Q>tstart           ~ as.numeric(NA),
        time2 >= time ~ 1,
        time2 <  time ~ 0
      ))
    }
    else {
      data_long_all <- data_long_all %>% mutate(at_risk = case_when(
        time2 >= time ~ 1,
        time2 <  time ~ 0
      ))
    }
    
    data_long_all <- data_long_all %>% mutate(Y2 = case_when(
      time2==time & Y == 1 ~ 1,
      time2==time & Y == 0 ~ 0,
      time2>time          ~ 0,
      time2<time          ~ as.numeric(NA)
    ))
  }
  else {
    #Defining indicator for outcome occurring in that interval
    if (LT_data == 1){
      data_long_all <- data_long_all %>% mutate(Y2 = case_when(
        Q>tstart           ~ as.numeric(NA),
        Q<tstart & time2<=time & time2>tstart & Y == 1 ~ 1,
        Q<tstart & time2<=time & time2>tstart & Y == 0 ~ 0,
        Q<tstart & time2>time          ~ 0,
        Q<tstart & time2<time          ~ as.numeric(NA)
      ))
    }
    else {
      data_long_all <- data_long_all %>% mutate(Y2 = case_when(
        time2<=time & time2>tstart & Y == 1 ~ 1,
        time2<=time & time2>tstart & Y == 0 ~ 0,
        time2>time          ~ 0,
        time2<time          ~ as.numeric(NA)
      ))
    }
   
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

      # if (is.null(time_cuts) == 1){
      #   #Defining at risk indicators for this dataset
      #   data_long_all <- data_long_all %>% mutate(at_risk = case_when(
      #     time2 >= time ~ 1,
      #     time2 <  time ~ 0
      #   ))
      # 
      #   data_long_all <- data_long_all %>% mutate(Y2 = case_when(
      #     time2==time & Y == 1 ~ 1,
      #     time2==time & Y == 0 ~ 0,
      #     time2>time          ~ 0,
      #     time2<time          ~ as.numeric(NA)
      #   ))
      # }
      # else {
      #   #Defining indicator for outcome occurring in that interval
      #   data_long_all <- data_long_all %>% mutate(Y2 = case_when(
      #     time2<=time & time2>tstart & Y == 1 ~ 1,
      #     time2<=time & time2>tstart & Y == 0 ~ 0,
      #     time2>time          ~ 0,
      #     time2<time          ~ as.numeric(NA)
      #   ))
      #   #Defining at risk indicators for this dataset
      #   data_long_all <- data_long_all %>% mutate(at_risk = case_when(
      #     is.na(Y2) == 1 ~ 0,
      #     is.na(Y2) == 0 ~ 1
      #   ))
      # }

      if (learner == "surviTMLE-learner" | learner == "M-learner"){
        new_data_vars <- c("ID",pse_covariates,"time")
        newdata_long_all <- subset(newdata_long_all,select=new_data_vars)
      }
      if (learner == "T-learner"){
        new_data_vars <- c("ID",out_covariates,"tstart","time","A","Y")
        newdata_long_all <- subset(newdata_long_all,select=new_data_vars)
      }
    },
    error=function(e) {
      stop('An error occurred when formatting the new data')
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
      
      if (LT_data == 1){
        #Checking if observed time is after truncation time
        if (sum((data$Y < data$Q)) == 0){
          stop("Truncation times must be smaller or equal to event/censoring times")
        }
      }
    },
    error=function(e) {
      stop('Either the exposure indicator is non-binary, the outcome indicator is non-binary, the censoring indicator is non-binary or the truncation times are larger than the event/censoring times')
      print(e)
    }
  )

  #-----------------------------#
  #--- Returning information ---#
  #-----------------------------#
  if (learner == "surviTMLE-learner" | learner == "M-learner" | learner == "T-learner" | learner == "CSF"){
    if (time_cuts[1] == "N/A"){
      output <- list(data=data,
                     data_long_all=data_long_all,
                     evt_times_uni=evt_times_uni,
                     evt_times_uni0=evt_times_uni0,
                     newdata=newdata,
                     newdata_long_all=newdata_long_all,
                     LT_data=LT_data,
                     e_covariates=e_covariates,
                     out_covariates=out_covariates,
                     g_covariates=g_covariates,
                     h_covariates=h_covariates,
                     pse_covariates=pse_covariates)
    }
    else {
      output <- list(data=data,
                     data_long_all=data_long_all,
                     evt_times_uni=time_cuts,
                     evt_times_uni0=time_cuts0,
                     newdata=newdata,
                     newdata_long_all=newdata_long_all,
                     LT_data=LT_data,
                     e_covariates=e_covariates,
                     out_covariates=out_covariates,
                     g_covariates=g_covariates,
                     h_covariates=h_covariates,
                     pse_covariates=pse_covariates)
    }
  }
  return(output)
}


