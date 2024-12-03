#######################################################################################
# Script: Nuisance model function
# Date: 21/03/24
# Author: Matt Pryce 
# Notes: - Used to run each nuisance model
#        - Runs parametric models, RF or SL 
#######################################################################################

library(survival)
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


#' @param model Nuisance model to be run
#' @param data The data frame containing all required information
#' @param method Statistical technique used to run the model
#' @param covariates  List containing the names of the variables to be input into the model
#' @param pred_data New data to create predictions for

#UPDATE

#' @param SL_lib Library to be used in super learner if selected
#' @param Y_bin Indicator for when the outcome is binary
#' @param Y_cont Indicator for when the outcome is continuous 

#' @return Outcome models and their predictions for pred_data 


nuis_mod_surv <- function(model,
                          data,
                          method = c("Parametric","Random forest","Super learner"),
                          covariates,
                          pred_data,
                          pred_data_long_all,
                          evt_times_uni,
                          SL_lib,
                          learner = NA,
                          CI_tuned_params,
                          LT = NULL
){
  
  #Allowing models which do not account for left truncation to continue to run 
  if (is.null(LT) == 1){
    LT <- 0
  }
  
  #------------------------------#
  #--- Creating training data ---#
  #------------------------------#
  
  tryCatch(
    {
      if (model == "Outcome" & LT == 0){
        train_data0 <- subset(data,data$A==0)
        train_data1 <- subset(data,data$A==1)
        train_data0 <- subset(train_data0,select = c("ID","time","Y",covariates))
        train_data1 <- subset(train_data1,select = c("ID","time","Y",covariates))
        
        covs_train0 <- train_data0[,covariates]
        covs_train1 <- train_data1[,covariates]
        
        pred_data_long_all <- pred_data_long_all %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Outcome" & LT == 1){
        train_data0 <- subset(data,data$A==0)
        train_data1 <- subset(data,data$A==1)
        train_data0 <- subset(train_data0,select = c("ID","time","Y",covariates))      
        train_data1 <- subset(train_data1,select = c("ID","time","Y",covariates))        
        
        covs_train0 <- train_data0[,covariates]
        covs_train1 <- train_data1[,covariates]
        
        pred_data_long_all <- pred_data_long_all %>%
          arrange(ID, tstart) %>%
          group_by(ID)
        
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Propensity score"){
        train_data <- as.data.frame(subset(data,select = c(covariates,"A","s")))
      }
      else if (model == "Censoring" & LT == 0){
        train_data <- as.data.frame(subset(data,select = c("time","C",covariates,"A")))
        covs_train <- train_data[,covariates]
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Censoring" & LT == 1){
        train_data <- as.data.frame(subset(data,select = c("time","C","Q",covariates,"A")))  
        covs_train <- train_data[,covariates]
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Truncation"){
        train_data <- as.data.frame(subset(data,select = c("time","Q",covariates,"A")))  
        covs_train <- train_data[,covariates]
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Pseudo outcome - Pooled" | model == "Pseudo outcome - Pooled - CI"){
        data <- data %>% mutate(pse_Y_cond = pse_Y - lag(pse_Y, default = 0))
        data <- subset(data,select=-c(pse_Y))
        data <- data %>% group_by(ID) %>% 
          summarize(across(everything(), last),time=0,tstart=0,pse_Y_cond=0, .groups = "drop") %>%
          bind_rows(data, .) %>% arrange(ID)
        data <- data %>% arrange(ID, time) %>% group_by(ID)
        train_data <- as.data.frame(subset(data,select = c("pse_Y_cond","time",covariates)))
      }
    },
    error=function(e) {
      stop('An error occured when creating analysis data')
      print(e)
    }
  )
  
  #----------------------#
  #--- Running models ---#
  #----------------------#

  if (model == "Outcome" & LT == 0){
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          mod_0 <- coxph(Surv(time, Y) ~ ., data = train_data0)
          mod_1 <- coxph(Surv(time, Y) ~ ., data = train_data1)
        }
        if (method == "Super learner"){
          # set.seed(1)  may need to set seed
          mod_0 <- survSuperLearner(time = train_data0$time, event = train_data0$Y, X = covs_train0, new.times = 10,
                                   event.SL.library = SL_lib, cens.SL.library = SL_lib, verbose = TRUE,
                                   cvControl=list(V = 5),
                                   control=list(saveFitLibrary = TRUE))
          mod_1 <- survSuperLearner(time = train_data1$time, event = train_data1$Y, X = covs_train1, new.times = 10,
                                    event.SL.library = SL_lib, cens.SL.library = SL_lib, verbose = TRUE,
                                    cvControl=list(V = 5),
                                    control=list(saveFitLibrary = TRUE))
        }
      },
      error=function(e) {
        stop('An error occured when running outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Outcome" & LT == 1){
    tryCatch(
      {
        if (method == "Parametric"){   
          #Running first outcome models
          mod_0 <- coxph(Surv(time, Y) ~ ., data = train_data0)
          mod_1 <- coxph(Surv(time, Y) ~ ., data = train_data1)
        }
        if (method == "Super learner"){   #GLobal survival stacking
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X0 <- train_data0[,c(covariates)]
          train_X1 <- train_data1[,c(covariates)]
          
          mod_0 <- survML::stackG(time = train_data0$time,
                                  event = train_data0$Y,
                                  entry = train_data0$Q,
                                  X = train_X0,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,
                                  time_basis = "continuous",
                                  time_grid_approx = sort(unique(train_data0$time)),
                                  surv_form = "exp",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))
          
          mod_1 <- survML::stackG(time = train_data1$time,
                                  event = train_data1$Y,
                                  entry = train_data0$Q,
                                  X = train_X1,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,
                                  time_basis = "continuous",
                                  time_grid_approx = sort(unique(train_data1$time)),
                                  surv_form = "exp",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))
        }
      },
      error=function(e) {
        stop('An error occured when running outcome models')
        print(e)
      }
    )
  }

  if (model == "Censoring" & LT == 0){
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          mod <- coxph(Surv(time, C) ~ ., data = train_data)
        }
        else if (method == "Super learner"){
          mod <- survSuperLearner(time = train_data$time, event = train_data$C, X = covs_train, new.times = 10,  #Check new times
                                    event.SL.library = SL_lib, cens.SL.library = SL_lib, verbose = TRUE,
                                    cvControl=list(V = 5),
                                    control=list(saveFitLibrary = TRUE))
        }
      },
      error=function(e) {
        stop('An error occured when running censoring model')
        print(e)
      }
    )
  }
  
  if (model == "Censoring" & LT == 1){   
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          mod <- coxph(Surv(time, C) ~ ., data = train_data)
        }
        else if (method == "Super learner"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X <- train_data[,c(covariates)]
          
          mod <- survML::stackG(time = train_data$time,
                                event = train_data$C,
                                entry = train_data$Q,
                                X = train_X,
                                newX = test_X,
                                newtimes = evt_times_uni,
                                bin_size = 0.025,
                                time_basis = "continuous",
                                time_grid_approx = sort(unique(train_data$time)),
                                surv_form = "exp",
                                SL_control = list(SL.library = SL_lib,
                                                  V = 5))
        }
      },
      error=function(e) {
        stop('An error occured when running censoring model')
        print(e)
      }
    )
  }
  
  if (model == "Truncation"){  #Survival stacking alg doesnt work as no censoring, but apply it manually
    tryCatch(
      {
        if (method == "Parametric"){     #Estimates probabily of Q greater than t, hence we take 1-pred later
          #Running first outcome models
          train_data$Q_ind <- 1     #Set event indicator to 1 as no one is censored
          mod <- coxph(Surv(Q, Q_ind) ~ ., data = train_data)
        }
        else if (method == "Super learner"){   #Estimates probabily of Q or equal to t, but will output conditional probs 
          #--- Manually implementing global survival stacking as it doesn't work when no censoring ---#
          evt_times_uni0 <- c(0,evt_times_uni)
          for (i in seq_along(evt_times_uni)){
            temp_data <- train_data
            temp_data$Q_ind <- 0
            for (j in 1:dim(temp_data)[1]){
              if (temp_data$Q[j] <= evt_times_uni0[i+1] & temp_data$Q[j] > evt_times_uni0[i]){
                temp_data$Q_ind[j] <- 1
              }  
            }
            temp_data$time <- evt_times_uni[i]
            if (i == 1){
              temp_data_long <- temp_data
            }
            else {
              temp_data_long <- rbind(temp_data_long,temp_data)
            }
          }
          
          #Runs pooled model (includes time as factor)
          temp_data_long$time <- as.factor(temp_data_long$time)
          
          mod <- SuperLearner(Y = temp_data_long$Q_ind, X = data.frame(subset(temp_data_long, select = c(covariates,"time"))),
                              method = "method.NNLS",
                              family = binomial(),
                              cvControl = list(V = 10, stratifyCV=FALSE),
                              SL.library = SL_lib)
        }
      },
      error=function(e) {
        stop('An error occured when running censoring model')
        print(e)
      }
    )
  }

  if (model == "Propensity score"){
    if (method == "Random forest"){
      X <- as.matrix(subset(train_data, select = covariates))
        mod <- regression_forest(X, train_data$A, honesty = FALSE,tune.parameters = "all")
    }
    else if (method == "Parametric"){
      fit_data <- subset(train_data,select = -c(s))
      mod <- glm(A ~ . , data = fit_data, family = binomial())
    }
    else if (method == "Super learner"){
      sums <- table(train_data$A)
      cv_folds <- min(10,sums[1],sums[2])
      mod <- SuperLearner(Y = train_data$A, X = data.frame(subset(train_data, select = covariates)),
                          method = "method.NNLS",
                          family = binomial(),
                          cvControl = list(V = cv_folds, stratifyCV=TRUE),
                          SL.library = SL_lib)
    }
  }

  if (model == "Pseudo outcome - Pooled"){
    if (method == "Random forest"){
      train_data$time <- as.numeric(train_data$time)
      X <- as.matrix(subset(train_data, select = c(covariates,"time")))
      mod <- regression_forest(X, train_data$pse_Y_cond)
    }
    else if (method == "Parametric"){
      train_data$time <- as.factor(train_data$time)
      mod <- lm(pse_Y_cond ~ . , data = train_data)
    }
    else if (method == "Super learner"){
      train_data$time <- as.factor(train_data$time)
      mod <- SuperLearner(Y = train_data$pse_Y_cond, X = data.frame(subset(train_data, select = c(covariates,"time"))),
                          method = "method.NNLS",
                          family = gaussian(),
                          cvControl = list(V = 10, stratifyCV=FALSE),
                          SL.library = SL_lib)
    }
  }

  if (model == "Pseudo outcome - Pooled - CI"){
    if (method == "Random forest"){
      #Defining tuning parameter values
      tuned_sample.fraction <- CI_tuned_params$sample.fraction
      tuned_mtry <- CI_tuned_params$mtry
      tuned_min.node.size <- CI_tuned_params$min.node.size
      tuned_honesty.fraction <- CI_tuned_params$honesty.fraction
      tuned_honesty.prune.leaves <- CI_tuned_params$honesty.prune.leaves
      tuned_alpha <- CI_tuned_params$alpha
      tuned_imbalance.penalty <- CI_tuned_params$imbalance.penalty

      train_data$time <- as.numeric(train_data$time)
      X <- as.matrix(subset(train_data, select = c(covariates,"time")))
      mod <- regression_forest(X, train_data$pse_Y_cond,
                               sample.fraction = tuned_sample.fraction,
                               mtry = tuned_mtry,
                               min.node.size = tuned_min.node.size,
                               honesty.fraction = tuned_honesty.fraction,
                               honesty.prune.leaves = tuned_honesty.prune.leaves,
                               alpha = tuned_alpha,
                               imbalance.penalty = tuned_imbalance.penalty)
    }
    else if (method == "Super learner"){
      #Could be added
    }
  }


  #-----------------------------------#
  #--- Obtaining model predictions ---#
  #-----------------------------------#

  if (model == "Outcome" & LT == 0){
    tryCatch(
      {
        #--- Obtaining preds from outcome models ---#
        if (method == "Parametric"){
          #Survival estimates at the time point (Training and all data points)
          S_k_pred_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="survival")
          S_k_pred_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="survival")

          if (learner == "survEP-learner" | learner == "M-learner"){
            #Hazard estimates at the time point (All time point data)
            pred_data_long_all$H_k_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="expected")
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, tstart) %>%
              group_by(ID) %>%
              mutate(h_k_long_all_0 = c(H_k_long_all_0[1], diff(H_k_long_all_0)))

            H_k_pred_long_all_0 <- pred_data_long_all$H_k_long_all_0
            h_k_pred_long_all_0 <- pred_data_long_all$h_k_long_all_0

            pred_data_long_all$H_k_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="expected")
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, tstart) %>%
              group_by(ID) %>%
              mutate(h_k_long_all_1 = c(H_k_long_all_1[1], diff(H_k_long_all_1)))

            H_k_pred_long_all_1 <- pred_data_long_all$H_k_long_all_1
            h_k_pred_long_all_1 <- pred_data_long_all$h_k_long_all_1


            #Survival predictions for time point of interest
            temp_data_list <- vector("list", length(evt_times_uni))
            for (t in seq_along(evt_times_uni)){
              temp_data <- pred_data_long_all
              temp_data$time2 <- temp_data$time      #Need time to be set to a number but don't want to lose original value
              temp_data$time <- evt_times_uni[t]

              temp_data_list[[t]] <- temp_data
            }

            S_t_pred_long_all_0_list <- lapply(temp_data_list, function(x){predict(mod_0,newdata=x,type="survival")})
            S_t_pred_long_all_1_list <- lapply(temp_data_list, function(x){predict(mod_1,newdata=x,type="survival")})

            for (t in seq_along(evt_times_uni)){
              temp_data <- temp_data_list[[t]]
              temp_data$S_t_pred_long_all_0 <- S_t_pred_long_all_0_list[[t]]
              temp_data$S_t_pred_long_all_1 <- S_t_pred_long_all_1_list[[t]]
              temp_data$time <- temp_data$time2       #Reassigning original value

              temp_data <- temp_data %>% mutate(S_t_pred_long_all_a = case_when(
                A == 1 ~ S_t_pred_long_all_1,
                A == 0 ~ S_t_pred_long_all_0
              ))

              temp_data_list[[t]] <- temp_data
            }
          }
        }
        if (method == "Super learner"){
          #Number of blocks
          num_blocks <- ceiling(dim(pred_data_long_all)[1]/10000)

          #Give observations numbers by 10,000
          number_list <- rep(1:num_blocks, each = 10000)
          pred_data_long_all$iter <- number_list[1:dim(pred_data_long_all)[1]]

          #Loop over the amount of 10000s
          temp_data_list_all <- replicate(length(evt_times_uni), replicate(num_blocks, list(), simplify = FALSE), simplify = FALSE)   #Make matrix
          for (i in 1:num_blocks){
            #Creating subset
            pred_data_long_all_temp <- subset(pred_data_long_all,pred_data_long_all$iter == i)

            #Creating pred data
            covs_test <- pred_data_long_all_temp[,covariates]

            #Creating survival predication
            S_k_pred_long_all_0 <- as.data.frame(predict.survSuperLearner(mod_0, newdata = covs_test, new.times=evt_times_uni)$event.SL.predict)
            S_k_pred_long_all_1 <- as.data.frame(predict.survSuperLearner(mod_1, newdata = covs_test, new.times=evt_times_uni)$event.SL.predict)

            #Rename columns of this matrix
            names(S_k_pred_long_all_0)[seq_along(evt_times_uni)] <- paste("S0_", evt_times_uni,sep="")
            names(S_k_pred_long_all_1)[seq_along(evt_times_uni)] <- paste("S1_", evt_times_uni,sep="")

            #Merging with long data
            pred_data_long_all_temp <- cbind(pred_data_long_all_temp,S_k_pred_long_all_0)
            pred_data_long_all_temp <- cbind(pred_data_long_all_temp,S_k_pred_long_all_1)

            if (learner == "T-learner"){
              #Defining S_k predictions for the dataset
              pred_data_long_all_temp$time_seq <- match(pred_data_long_all_temp$time, evt_times_uni)
              col_index_offset <- 6 + length(covariates)
              pred_data_long_all_temp$S_k_pred_0 <- apply(pred_data_long_all_temp, 1, function(row) {
                row[col_index_offset + row['time_seq']]
              })

              col_index_offset <- 6 + length(covariates) + length(evt_times_uni)
              pred_data_long_all_temp$S_k_pred_1 <- apply(pred_data_long_all_temp, 1, function(row) {
                row[col_index_offset + row['time_seq']]
              })
            }
            if (learner == "survEP-learner" | learner == "M-learner"){
              #Defining S_k predictions for the dataset
              pred_data_long_all_temp$time_seq <- match(pred_data_long_all_temp$time, evt_times_uni)
              col_index_offset <- 9 + length(covariates)
              pred_data_long_all_temp$S_k_pred_0 <- apply(pred_data_long_all_temp, 1, function(row) {
                row[col_index_offset + row['time_seq']]
              })

              col_index_offset <- 9 + length(covariates) + length(evt_times_uni)
              pred_data_long_all_temp$S_k_pred_1 <- apply(pred_data_long_all_temp, 1, function(row) {
                row[col_index_offset + row['time_seq']]
              })

              #Calculating cumulative hazard estimates
              pred_data_long_all_temp$H_k_pred_0 <- -log(pred_data_long_all_temp$S_k_pred_0)
              pred_data_long_all_temp$H_k_pred_1 <- -log(pred_data_long_all_temp$S_k_pred_1)

              #Calculating hazards at each time
              pred_data_long_all_temp <- pred_data_long_all_temp %>%
                arrange(ID, tstart) %>%
                group_by(ID) %>%
                mutate(h_k_pred_0 = c(H_k_pred_0[1], diff(H_k_pred_0)))
              pred_data_long_all_temp <- pred_data_long_all_temp %>%
                arrange(ID, tstart) %>%
                group_by(ID) %>%
                mutate(h_k_pred_1 = c(H_k_pred_1[1], diff(H_k_pred_1)))

              #Survival predictions for time point of interest
              for (t in seq_along(evt_times_uni)){
                temp_data <- pred_data_long_all_temp
                temp_data$S_t_pred_long_all_0 <- temp_data[[9 + length(covariates) + t]]
                temp_data$S_t_pred_long_all_1 <- temp_data[[9 + length(covariates) + length(evt_times_uni) + t]]

                temp_data <- temp_data %>% mutate(S_t_pred_long_all_a = case_when(
                  A == 1 ~ S_t_pred_long_all_1,
                  A == 0 ~ S_t_pred_long_all_0
                ))

                keep_covs_temp <- append(c("ID","Y","A","C","tstart","time","s",
                                           "S_k_pred_0","S_k_pred_1","H_k_pred_0","H_k_pred_1","h_k_pred_0","h_k_pred_1","at_risk",
                                           "S_t_pred_long_all_0","S_t_pred_long_all_1","S_t_pred_long_all_a"),covariates)
                temp_data <- subset(temp_data,select = keep_covs_temp)

                #Allocate to correct bit of matrix
                temp_data_list_all[[t]][[i]] <- temp_data
              }
            }
            if (i == 1){
              pred_data_long_all_final <- pred_data_long_all_temp
            }
            else {
              pred_data_long_all_final <- rbind(pred_data_long_all_final,pred_data_long_all_temp)
            }
          }
          #Loop over things to create full datasets
          temp_data_list <- vector("list", length(evt_times_uni))
          for (tp in 1:length(evt_times_uni)){
            temp_data_list[[tp]] <- do.call(rbind, temp_data_list_all[[tp]])
          }

          pred_data_long_all <- pred_data_long_all_final

          if (learner == "survEP-learner" | learner == "M-learner"){
          keep_covs <- append(c("ID","Y","A","C","tstart","time","s",
                                "S_k_pred_0","S_k_pred_1","H_k_pred_0","H_k_pred_1","h_k_pred_0","h_k_pred_1","at_risk"),covariates)
          }
          else if (learner == "T-learner"){
            keep_covs <- append(c("ID","Y","tstart","time",
                                  "S_k_pred_0","S_k_pred_1"),covariates)
          }
          pred_data_long_all <- subset(pred_data_long_all,select = keep_covs)
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }
  
  
  if (model == "Outcome" & LT == 1){    
    tryCatch(
      {
        #--- Obtaining preds from outcome models ---#
        if (method == "Parametric"){
          #Survival estimates at the time point (Training and all data points)
          S_k_pred_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="survival")
          S_k_pred_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="survival")

          if (learner == "survEP-learner"){
            #Hazard estimates at the time point (All time point data)
            pred_data_long_all$H_k_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="expected")
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, tstart) %>%
              group_by(ID) %>%
              mutate(h_k_long_all_0 = c(H_k_long_all_0[1], diff(H_k_long_all_0)))

            H_k_pred_long_all_0 <- pred_data_long_all$H_k_long_all_0
            h_k_pred_long_all_0 <- pred_data_long_all$h_k_long_all_0

            pred_data_long_all$H_k_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="expected")
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, tstart) %>%
              group_by(ID) %>%
              mutate(h_k_long_all_1 = c(H_k_long_all_1[1], diff(H_k_long_all_1)))

            H_k_pred_long_all_1 <- pred_data_long_all$H_k_long_all_1
            h_k_pred_long_all_1 <- pred_data_long_all$h_k_long_all_1


            #Survival predictions for time point of interest
            temp_data_list <- vector("list", length(evt_times_uni))
            for (t in seq_along(evt_times_uni)){
              temp_data <- pred_data_long_all
              temp_data$time2 <- temp_data$time      #Need time to be set to a number but don't want to lose original value
              temp_data$time <- evt_times_uni[t]

              temp_data_list[[t]] <- temp_data
            }

            S_t_pred_long_all_0_list <- lapply(temp_data_list, function(x){predict(mod_0,newdata=x,type="survival")})
            S_t_pred_long_all_1_list <- lapply(temp_data_list, function(x){predict(mod_1,newdata=x,type="survival")})

            for (t in seq_along(evt_times_uni)){
              temp_data <- temp_data_list[[t]]
              temp_data$S_t_pred_long_all_0 <- S_t_pred_long_all_0_list[[t]]
              temp_data$S_t_pred_long_all_1 <- S_t_pred_long_all_1_list[[t]]
              temp_data$time <- temp_data$time2       #Reassigning original value

              temp_data <- temp_data %>% mutate(S_t_pred_long_all_a = case_when(
                A == 1 ~ S_t_pred_long_all_1,
                A == 0 ~ S_t_pred_long_all_0
              ))

              temp_data_list[[t]] <- temp_data
            }
          }
        }
        if (method == "Super learner"){
          uni_time <- rep(evt_times_uni, times = nrow(pred_data))
          
          #Obtaining predictions from outcome models 
          S_k_pred_wide_all_0 <- as.data.frame(mod_0$S_T_preds)
          S_k_pred_wide_all_1 <- as.data.frame(mod_1$S_T_preds)
          
          #Making into long form 
          ncols_temp <- ncol(pred_data)
          pred_data0 <- cbind(pred_data,S_k_pred_wide_all_0)
          pred_data1 <- cbind(pred_data,S_k_pred_wide_all_1)
          
          pred_data0 <- pred_data0 %>%
            pivot_longer(
              cols = starts_with("V"),  
              names_to = "time2",           
              values_to = "S_k_pred_0"         
            )
          pred_data1 <- pred_data1 %>%
            pivot_longer(
              cols = starts_with("V"),  
              names_to = "time2",         
              values_to = "S_k_pred_1"        
            )
          pred_data_long_all <- pred_data0
          pred_data_long_all$S_k_pred_1 <- pred_data1$S_k_pred_1
          pred_data_long_all$time <- uni_time
          S_k_pred_long_all_0 <- pred_data_long_all$S_K_pred_0
          S_k_pred_long_all_1 <- pred_data_long_all$S_K_pred_1
          
          if (learner == "survEP-learner"){
            #Calculating cumulative hazard estimates
            pred_data_long_all$H_k_pred_0 <- -log(pred_data_long_all$S_k_pred_0)
            pred_data_long_all$H_k_pred_1 <- -log(pred_data_long_all$S_k_pred_1)
            
            #Calculating hazards at each time
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, time) %>%
              group_by(ID) %>%
              mutate(h_k_pred_0 = c(H_k_pred_0[1], diff(H_k_pred_0)))
            pred_data_long_all <- pred_data_long_all %>%
              arrange(ID, time) %>%
              group_by(ID) %>%
              mutate(h_k_pred_1 = c(H_k_pred_1[1], diff(H_k_pred_1)))
            
            #Survival predictions for time point of interest
            temp_data_list <- vector("list", length(evt_times_uni))
            preds0 <- cbind(ID=pred_data$ID,S_k_pred_wide_all_0)
            preds1 <- cbind(ID=pred_data$ID,S_k_pred_wide_all_1)
            for (t in seq_along(evt_times_uni)){
              #Collecting appropriate survival preds 
              pred0 <- preds0[,c(1,t+1)]
              names(pred0)[2] <- "S_t_pred_long_all_0"
              pred1 <- preds1[,c(1,t+1)]
              names(pred1)[2] <- "S_t_pred_long_all_1"
              
              #Merging with long pred data
              pred_data_long_all_temp <- merge(pred_data_long_all,pred0,by="ID",all.x=TRUE)
              pred_data_long_all_temp <- merge(pred_data_long_all_temp,pred1,by="ID",all.x=TRUE)
              
              pred_data_long_all_temp <- pred_data_long_all_temp %>% mutate(S_t_pred_long_all_a = case_when(
                A == 1 ~ S_t_pred_long_all_1,
                A == 0 ~ S_t_pred_long_all_0
              ))
              
              # Saving data to list
              temp_data_list[[t]] <- pred_data_long_all_temp
            }
          }
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }

  if (model == "Censoring" & LT == 0){
    tryCatch(
      {
        #--- Obtaining preds from censoring model ---#
        if (method == "Parametric"){
          #Censoring estimates (Training and all data points)
          G_k_pred_long_all <- predict(mod,newdata=pred_data_long_all,type="survival")
        }
        if (method == "Super learner"){
          if (dim(pred_data_long_all)[1] <= 10000){
            #Survival estimates at the time point (Training and all data points)
            G_k_pred_long_all <- as.data.frame(predict.survSuperLearner(mod, newdata = covs_test, new.times=evt_times_uni)$event.SL.predict)

            #Rename columns of this matrix
            names(G_k_pred_long_all)[seq_along(evt_times_uni)] <- paste("G_", evt_times_uni,sep="")

            #Merging with long data
            pred_data_long_all <- cbind(pred_data_long_all,G_k_pred_long_all)

            #Defining G_k predictions for the dataset
            pred_data_long_all$time_seq <- match(pred_data_long_all$time, evt_times_uni)
            col_index_offset <- 6 + length(covariates)
            pred_data_long_all$G_k_pred <- apply(pred_data_long_all, 1, function(row) {
              row[col_index_offset + row['time_seq']]
            })
          }
          else{
            #Number of blocks
            num_blocks <- ceiling(dim(pred_data_long_all)[1]/10000)

            #Give observations numbers by 10,000
            number_list <- rep(1:num_blocks, each = 10000)
            pred_data_long_all$iter <- number_list[1:dim(pred_data_long_all)[1]]

            #Loop over the amount of 10000s
            for (i in 1:num_blocks){
              #Creating subset
              pred_data_long_all_temp <- subset(pred_data_long_all,pred_data_long_all$iter == i)

              #Creating pred data
              covs_test <- pred_data_long_all_temp[,covariates]

              #Creating survival predication
              G_k_pred_long_all <- as.data.frame(predict.survSuperLearner(mod, newdata = covs_test, new.times=evt_times_uni)$event.SL.predict)

              #Rename columns of this matrix
              names(G_k_pred_long_all)[seq_along(evt_times_uni)] <- paste("G_", evt_times_uni,sep="")

              #Merging with long data
              pred_data_long_all_temp <- cbind(pred_data_long_all_temp,G_k_pred_long_all)

              #Defining S_k predictions for the dataset
              pred_data_long_all_temp$time_seq <- match(pred_data_long_all_temp$time, evt_times_uni)
              col_index_offset <- 7 + length(covariates)
              pred_data_long_all_temp$G_k_pred <- apply(pred_data_long_all_temp, 1, function(row) {
                row[col_index_offset + row['time_seq']]
              })

              keep_covs <- c("ID","A","C","tstart","time","s",covariates,"G_k_pred")
              pred_data_long_all_temp <- subset(pred_data_long_all_temp,select = keep_covs)

              if (i == 1){
                pred_data_long_all_final <- pred_data_long_all_temp
              }
              else {
                pred_data_long_all_final <- rbind(pred_data_long_all_final,pred_data_long_all_temp)
              }
            }
            pred_data_long_all <- pred_data_long_all_final
          }
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Censoring" & LT == 1){    
    tryCatch(
      {
        #--- Obtaining preds from censoring model ---#
        if (method == "Parametric"){
          #Censoring estimates (Training and all data points)
          G_k_pred_long_all <- predict(mod,newdata=pred_data_long_all,type="survival")
        }
        else if (method == "Super learner"){
          uni_time <- rep(evt_times_uni, times = nrow(pred_data))
          
          #Obtaining predictions from outcome models 
          G_k_pred_wide_all <- as.data.frame(mod$S_T_preds)
          
          #Making into long form 
          ncols_temp <- ncol(pred_data)
          pred_data <- cbind(pred_data,G_k_pred_wide_all)
          
          pred_data_temp <- pred_data %>%
            pivot_longer(
              cols = starts_with("V"),  
              names_to = "time2",           
              values_to = "G_k_pred"         
            )
          pred_data_long_all <- pred_data_temp
          pred_data_long_all$time <- uni_time
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Truncation"){  
    tryCatch(
      {
        #--- Obtaining preds from censoring model ---#
        if (method == "Parametric"){
          #Truncation estimates (Training and all data points)
          pred_data_long_all$Q_ind <- 1
          pred_data_long_all$Q <- pred_data_long_all$time
          trunc_temp_k_pred_long_all <- predict(mod,newdata=pred_data_long_all,type="survival")
          trunc_k_pred_long_all <- 1 - trunc_temp_k_pred_long_all
        }
        else if (method == "Super learner"){
          #Making predictions for at each time for each person
          pred_data <- subset(pred_data_long_all, select = c(covariates,"time"))
          pred_data$time <- as.factor(pred_data$time)
          pred_data_long_all$pred_cond <- as.vector(predict(mod, pred_data)$pred)
          
          #--- Creating Cumulative predictions ---#
          #Need to number time points
          tps <- unique(pred_data_long_all$time)
          num_tps <- length(tps)

          #Re-numbering time point
          tps2_list <- 0:(num_tps - 1)
          pred_data_long_all$time2 <- 0
          for (i in 0:(num_tps-1)){
            for (j in 1:length(pred_data_long_all$ID)){
              if (pred_data_long_all$time[j] == tps[i+1]){
                pred_data_long_all$time2[j] <- tps2_list[i+1]
              }
            }
          }
          pred_data_long_all$time_ordered <- pred_data_long_all$time2 + 1

          pred_data_long_all$pred <- pred_data_long_all$pred_cond
          for (k in 1:length(pred_data_long_all$ID)){
            if (pred_data_long_all$time_ordered[k] > 1){
              pred_data_long_all$pred[k] <- pred_data_long_all$pred[k-1]+pred_data_long_all$pred_cond[k]
            }
          }
          pred_data_long_all$pred[pred_data_long_all$pred > 1] <- 1
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }

  if (model == "Propensity score"){
    if (method == "Random forest"){
      pred <- predict(mod, pred_data)$predictions
    }
    else if (method == "Parametric"){
      pred <- predict(mod, as.data.frame(pred_data), type = "response")
    }
    else if (method == "Super learner"){
      pred <- as.vector(predict(mod, as.data.frame(pred_data))$pred)
    }
  }

  if (model == "Pseudo outcome - Pooled" | model == "Pseudo outcome - Pooled - CI"){
    pred_data_long_all <- pred_data_long_all[order(pred_data_long_all$ID,pred_data_long_all$time),]
    pred_data <- subset(pred_data_long_all, select = c(covariates,"time"))
    pred_data$time <- as.factor(pred_data$time)
    if (method == "Parametric"){
      #Obtaining conditional predictions
      pred_data_long_all$pred_cond <- predict(mod, pred_data, type = "response")
    }
    else if (method == "Random forest"){
      pred_data$time <- as.numeric(pred_data$time)
      pred_data_matrix <- as.matrix(pred_data)
      pred_data_long_all$pred_cond <- predict(mod, pred_data_matrix)$predictions
    }
    else if (method == "Super learner"){
      pred_data_long_all$pred_cond <- predict(mod, pred_data)$pred
    }

    #--- Creating Cumulative predictions ---#
    #Need to number time points
    tps <- unique(pred_data_long_all$time)
    num_tps <- length(tps)

    #Re-numbering time point
    tps2_list <- 0:(num_tps - 1)
    pred_data_long_all$time2 <- 0
    for (i in 0:(num_tps-1)){
      for (j in 1:length(pred_data_long_all$ID)){
        if (pred_data_long_all$time[j] == tps[i+1]){
          pred_data_long_all$time2[j] <- tps2_list[i+1]
        }
      }
    }
    pred_data_long_all$time_ordered <- pred_data_long_all$time2 + 1

    pred_data_long_all$pred <- pred_data_long_all$pred_cond
    for (k in 1:length(pred_data_long_all$ID)){
      if (pred_data_long_all$time_ordered[k] > 1){
        pred_data_long_all$pred[k] <- pred_data_long_all$pred[k-1]+pred_data_long_all$pred_cond[k]
      }
    }
  }



  #-----------------------------#
  #--- Returning information ---#    #Make sure LT output aligns or gets specified 
  #-----------------------------#

  if (model == "Outcome" & method == "Super learner"  & (learner == "survEP-learner" | learner == "M-learner")){
      output <- list(out_mod_0 = mod_0,
                     out_mod_1 = mod_1,
                     pred_data_long_all_pred = pred_data_long_all,
                     Surv_t_pred_long_all_list=temp_data_list)
  }
  else if (model == "Outcome" & method != "Super learner"  & (learner == "survEP-learner" | learner == "M-learner")){
    output <- list(S_k_pred_long_all_0 = S_k_pred_long_all_0,
                   S_k_pred_long_all_1 = S_k_pred_long_all_1,
                   H_k_pred_long_all_0 = H_k_pred_long_all_0,
                   h_k_pred_long_all_0 = h_k_pred_long_all_0,
                   H_k_pred_long_all_1 = H_k_pred_long_all_1,
                   h_k_pred_long_all_1 = h_k_pred_long_all_1,
                   out_mod_0 = mod_0,
                   out_mod_1 = mod_1,
                   Surv_t_pred_long_all_list=temp_data_list)
  }
  else if (model == "Outcome" & learner == "T-learner" & method == "Parametric"){
    output <- list(S_k_pred_long_all_0 = S_k_pred_long_all_0,
                   S_k_pred_long_all_1 = S_k_pred_long_all_1,
                   out_mod_0 = mod_0,
                   out_mod_1 = mod_1)
  }
  else if (model == "Outcome" & learner == "T-learner" & method == "Super learner"){
    output <- list(pred_data_long_all_pred = pred_data_long_all,
                   out_mod_0 = mod_0,
                   out_mod_1 = mod_1)
  }
  else if (model == "Propensity score"){
    output <- list(e_pred = pred,
                   e_mod = mod)
  }
  else if (model == "Censoring" & method == "Parametric"){
    output <- list(G_k_pred_long_all = G_k_pred_long_all,
                   G_mod = mod)
  }
  else if (model == "Censoring" & method == "Super learner"){
    output <- list(G_k_pred_long_all = pred_data_long_all,
                   G_mod = mod)
  }
  else if (model == "Truncation" & method == "Parametric"){
    output <- list(trunc_k_pred_long_all = trunc_k_pred_long_all,
                   trunc_mod = mod)
  }
  else if (model == "Truncation" & method == "Super learner"){
    output <- list(trunc_k_pred_long_all = pred_data_long_all,
                   trunc_mod = mod)
  }
  else if (model == "Pseudo outcome - Pooled"){
    output <- list(po_mod = mod,
                   pred = pred_data_long_all)
  }
  else if (model == "Pseudo outcome - Pooled - CI"){
    output <- list(pred = pred_data_long_all$pred)
  }

  return(output)
}


