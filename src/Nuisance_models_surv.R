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
                          method = c("Parametric","Random forest","Super learner","Local survival stack","Global survival stack"),
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
        train_data0 <- subset(train_data0,select = c("ID","time","Y","Q",covariates))      
        train_data1 <- subset(train_data1,select = c("ID","time","Y","Q",covariates))        
        
        covs_train0 <- train_data0[,covariates]
        covs_train1 <- train_data1[,covariates]
        
        if (learner != "ltrc"){
          pred_data_long_all <- pred_data_long_all %>%
            arrange(ID, tstart) %>%
            group_by(ID)
        }
        
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Outcome - Diff" & LT == 1){
        train_data <- data
        train_data <- subset(train_data,select = c("ID","time","cen_ind","Q",covariates))
        train_data$D <- train_data$time - train_data$Q

        covs_train <- train_data[,covariates]

        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Truncation weight" & LT == 1){
        train_data <- subset(data)
        train_data <- subset(train_data,select = c("ID","time","Y","Q",covariates))      
        
        covs_train <- train_data[,covariates]

        if (learner != "ltrc"){
          pred_data_long_all <- pred_data_long_all %>%
            arrange(ID, tstart) %>%
            group_by(ID)
        }
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Propensity score"){
        if (LT == 0){
          train_data <- as.data.frame(subset(data,select = c(covariates,"A","s")))
        }
        else if (LT == 1){
          train_data <- as.data.frame(subset(data,select = c(covariates,"A","trunc_weight_pred","s")))
        }
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
        train_data <- as.data.frame(subset(data,select = c("time","Q","trunc_weight_pred",covariates,"A")))  
        covs_train <- train_data[,covariates]
        covs_test <- pred_data_long_all[,covariates]
      }
      else if (model == "Pseudo outcome - Pooled - Factor" | model == "Pseudo outcome - Pooled - Continuous" | model == "Pseudo outcome - Pooled - CI"){
        data <- data %>% mutate(pse_Y_cond = pse_Y - lag(pse_Y, default = 0))
        data <- subset(data,select=-c(pse_Y))
        data_temp <- subset(data,data$tstart==0)
        data_temp$time <-  0
        data_temp$tstart <- 0
        data_temp$pse_Y_cond <- 0
        data <- rbind(data,data_temp)
        data <- data %>% arrange(ID, time) %>% group_by(ID)
        train_data <- as.data.frame(subset(data,select = c("ID","pse_Y_cond","time",covariates)))
      }
      else if (model == "Pseudo outcome - One - DR"){
        train_data <- as.data.frame(subset(data,select = c("pse_Y","weights_DR",covariates)))
      }
      else if (model == "Pseudo outcome - One - R"){
        train_data <- as.data.frame(subset(data,select = c("pse_Y","weights_R",covariates)))
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
        if (method == "Local survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X0 <- train_data0[,c(covariates)]
          train_X1 <- train_data1[,c(covariates)]
          train_data0$Q <- 0  #Setting truncation time to 0
          train_data1$Q <- 0
          
          mod_0 <- survML::stackL(time = train_data0$time,
                                  event = train_data0$Y,
                                  entry = train_data0$Q,
                                  X = train_X0,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,  
                                  time_basis = "continuous",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))

          mod_1 <- survML::stackL(time = train_data1$time,
                                  event = train_data1$Y,
                                  entry = train_data1$Q,
                                  X = train_X1,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,  
                                  time_basis = "continuous",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))
        }
        if (method == "Global survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X0 <- train_data0[,c(covariates)]
          train_X1 <- train_data1[,c(covariates)]
          
          mod_0 <- survML::stackG(time = train_data0$time,
                                  event = train_data0$Y,
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
  
  if (model == "Outcome" & LT == 1){
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          train_data0 <- subset(train_data0,select = c("time","Y",covariates))      
          train_data1 <- subset(train_data1,select = c("time","Y",covariates))  
          
          mod_0 <- coxph(Surv(time, Y) ~ ., data = train_data0)
          mod_1 <- coxph(Surv(time, Y) ~ ., data = train_data1)
        }
        if (method == "Local survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X0 <- train_data0[,c(covariates)]
          train_X1 <- train_data1[,c(covariates)]

          mod_0 <- survML::stackL(time = train_data0$time,
                                  event = train_data0$Y,
                                  entry = train_data0$Q,
                                  X = train_X0,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,
                                  time_basis = "continuous",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))

          mod_1 <- survML::stackL(time = train_data1$time,
                                  event = train_data1$Y,
                                  entry = train_data1$Q,
                                  X = train_X1,
                                  newX = test_X,
                                  newtimes = evt_times_uni,
                                  bin_size = 0.025,
                                  time_basis = "continuous",
                                  SL_control = list(SL.library = SL_lib,
                                                    V = 5))
        }
        if (method == "Global survival stack"){
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
        if (method == "pcox"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.matrix(pred_data[,c(covariates)])
          
          yss_0 = Surv(train_data0[,"time"], train_data0[,"Y"])
          yss_1 = Surv(train_data1[,"time"], train_data1[,"Y"])
          
          cv.fit_0 <- cv.glmnet(as.matrix(covs_train0), yss_0, family = "cox")
          cv.fit_1 <- cv.glmnet(as.matrix(covs_train1), yss_1, family = "cox")

          fit.pred_0 = survival::survfit(cv.fit_0, s = "lambda.min", x = as.matrix(covs_train0), y = yss_0, newx = test_X)
          fit.pred_0 <- summary(fit.pred_0, times = evt_times_uni) 
          fit.pred_1 = survival::survfit(cv.fit_1, s = "lambda.min", x = as.matrix(covs_train1), y = yss_1, newx = test_X)
          fit.pred_1 <- summary(fit.pred_1, times = evt_times_uni)
        }
      },
      error=function(e) {
        stop('An error occured when running outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Outcome - Diff" & LT == 1){
    tryCatch(
      {
        if (method == "Local survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X <- train_data[,c(covariates)]

          mod <- survML::stackL(time = train_data$D,
                                event = train_data$cen_ind,
                                X = train_X,
                                newX = test_X,
                                newtimes = evt_times_uni,
                                bin_size = 0.025,
                                time_basis = "continuous",
                                SL_control = list(SL.library = SL_lib,
                                                  V = 5))
        }
      },
      error=function(e) {
        stop('An error occured when running outcome diff model')
        print(e)
      }
    )
  }
  
  if (model == "Truncation weight" & LT == 1){
    tryCatch(
      {
        if (method == "Local survival stack"){
          train_X <- train_data[,c(covariates)]

          #Listing all truncation times
          trunc_times <- sort(unique(train_data$Q))

          mod <- survML::stackL(time = train_data$time,
                                event = train_data$Y,
                                entry = train_data$Q,
                                X = train_X,
                                newX = train_X,
                                newtimes = trunc_times,
                                bin_size = 0.025,
                                time_basis = "continuous",
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
        else if (method == "Local survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X <- train_data[,c(covariates)]
          train_data$Q <- 0
          
          mod <- survML::stackL(time = train_data$time,
                                event = train_data$C,
                                entry = train_data$Q,
                                X = train_X,
                                newX = test_X,
                                newtimes = evt_times_uni,
                                bin_size = 0.025,
                                time_basis = "continuous",
                                SL_control = list(SL.library = SL_lib,
                                                  V = 5))
        }
        else if (method == "Global survival stack"){
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X <- train_data[,c(covariates)]
          
          mod <- survML::stackG(time = train_data$time,
                                event = train_data$C,
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
  
  if (model == "Censoring" & LT == 1){   
    tryCatch(
      {
        if (method == "Parametric"){   #Does not include Q - Needs updating
          #Running first outcome models
          mod <- coxph(Surv(time, C) ~ ., data = train_data)
        }
        else if (method == "Local survival stack"){
          #Renaming truncation time in train data to align with pred
          train_data$Q2 <- train_data$Q
          
          #Need a dataset for predictions where each person at each time has Q set to that time
          pred_data_long_all$Q2 <- pred_data_long_all$time
          
          #Set up
          train_data$entry <- 0  #Setting entry to be 0 as censoring not truncated
          test_X <- as.data.frame(pred_data_long_all[,c(covariates,"Q2")])
          train_X <- train_data[,c(covariates,"Q2")]
          #pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]  #Used before edit
          
          mod <- survML::stackL(time = train_data$time,
                                event = train_data$C,
                                entry = train_data$entry,
                                X = train_X,
                                newX = test_X,
                                newtimes = evt_times_uni,
                                bin_size = 0.025,
                                time_basis = "continuous",
                                SL_control = list(SL.library = SL_lib,
                                                  V = 5))
        }
        else if (method == "Global survival stack"){
          #Renaming truncation time in train data to align with pred
          train_data$Q2 <- train_data$Q
          
          #Need a dataset for predictions where each person at each time has Q set to that time
          pred_data_long_all$Q2 <- pred_data_long_all$time
          
          #Set up
          train_data$entry <- 0  #Setting entry to be 0 as censoring not truncated
          test_X <- as.data.frame(pred_data_long_all[,c(covariates,"Q2")])
          train_X <- train_data[,c(covariates,"Q2")]
          
          mod <- survML::stackG(time = train_data$time,
                                event = train_data$C,
                                entry = train_data$entry,
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
          fit_data <- subset(train_data,select = -c(trunc_weight_pred))
          mod <- coxph(Surv(Q, Q_ind) ~ ., data = fit_data, weights = train_data$trunc_weight_pred)
        }
        else if (method == "Local survival stack"){   #Estimates probabily of Q or equal to t, but will output conditional probs
          pred_data <- pred_data_long_all[!duplicated(pred_data_long_all$ID), ]
          test_X <- as.data.frame(pred_data[,c(covariates)])
          train_X <- train_data[,c(covariates)]
          train_data$Q_ind <- 1
          train_data$no_trunc <- 0

          mod <- survML::stackL(time = train_data$Q,
                                event = train_data$Q_ind,
                                entry = train_data$no_trunc,
                                X = train_X,
                                newX = test_X,
                                newtimes = evt_times_uni,
                                bin_size = 0.025,
                                time_basis = "continuous",
                                SL_control = list(SL.library = SL_lib,
                                                  V = 5,
                                                  obsWeights = train_data$trunc_weight_pred))
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
      if (LT == 0){
        mod <- regression_forest(X, train_data$A, honesty = FALSE,tune.parameters = "all")
      }
      if (LT == 1){
        mod <- regression_forest(X, train_data$A, honesty = FALSE,tune.parameters = "all",sample.weights = train_data$trunc_weight_pred)
      }
    }
    else if (method == "Parametric"){
      if (LT == 0){
        fit_data <- subset(train_data,select = -c(s))
        mod <- glm(A ~ . , data = fit_data, family = binomial())
      }
      if (LT == 1){
        fit_data <- subset(train_data,select = -c(s,trunc_weight_pred))
        mod <- glm(A ~ . , data = fit_data, family = quasibinomial(), weights = train_data$trunc_weight_pred)
      }
    }
    else if (method == "Super learner"){
      sums <- table(train_data$A)
      cv_folds <- min(10,sums[1],sums[2])
      if (LT == 0){
        mod <- SuperLearner(Y = train_data$A, X = data.frame(subset(train_data, select = covariates)),
                            method = "method.NNLS",
                            family = binomial(),
                            cvControl = list(V = cv_folds, stratifyCV=TRUE),
                            SL.library = SL_lib)
      }
      if (LT == 1){
        mod <- SuperLearner(Y = train_data$A, X = data.frame(subset(train_data, select = covariates)),
                            method = "method.NNLS",
                            family = binomial(),
                            cvControl = list(V = cv_folds, stratifyCV=TRUE),
                            SL.library = SL_lib,
                            obsWeights = train_data$trunc_weight_pred)
      }
    }
  }

  if (model == "Pseudo outcome - Pooled - Factor" | model == "Pseudo outcome - Pooled - Continuous"){  
    if (method == "Parametric"){
      if (method == "Random forest"){   #Only continuous
        train_data$time <- as.numeric(train_data$time)
        X <- as.matrix(subset(train_data, select = c(covariates,"time")))
        mod <- regression_forest(X, train_data$pse_Y_cond)
      }
      if (model == "Pseudo outcome - Pooled - Factor"){
        train_data$time <- as.factor(train_data$time)
      }
      else if (model == "Pseudo outcome - Pooled - Continuous"){
        train_data$time <- as.numeric(train_data$time)
      }
      mod <- lm(pse_Y_cond ~ . , data = train_data)
    }
    else if (method == "Super learner"){
      if (model == "Pseudo outcome - Pooled - Factor"){
        train_data$time <- as.factor(train_data$time)
      }
      else if (model == "Pseudo outcome - Pooled - Continuous"){
        train_data$time <- as.numeric(train_data$time)
      }
      mod <- SuperLearner(Y = train_data$pse_Y_cond, X = data.frame(subset(train_data, select = c(covariates,"time"))),
                          method = "method.NNLS",
                          family = gaussian(),
                          id=train_data$ID,
                          cvControl = list(V = 10, stratifyCV=FALSE),
                          SL.library = SL_lib)
    }
    else if (method == "GAM"){
      if (model == "Pseudo outcome - Pooled - Factor"){
        train_data$time <- as.factor(train_data$time)
      }
      else if (model == "Pseudo outcome - Pooled - Continuous"){
        train_data$time <- as.numeric(train_data$time)
      }

      #Removing time from X
      data_train2 <- subset(train_data, select = covariates)
      
      find_binary_vars <- function(df) {
        binary_vars <- sapply(df, function(col) {
          vals <- unique(na.omit(col))
          length(vals) == 2
        })
        names(df)[binary_vars]
      }
      bin_vars <- find_binary_vars(data_train2)

      cont_vars <- setdiff(
        names(data_train2)[sapply(data_train2, function(v) is.numeric(v) | is.integer(v))],
        bin_vars
      )

      # Start formula string
      form_str <- "Y ~ s(time, k = 4)"

      # Add binary terms
      if (length(bin_vars) > 0) {
        bin_terms <- paste0("s(time, by = ",bin_vars,",k=4)")
        form_str <- paste(form_str, paste(bin_terms, collapse = " + "), sep = " + ")
      }

      # Add continuous terms
      if (length(cont_vars) > 0) {
        cont_terms <- paste0("te(time, ", cont_vars, ", k = c(4, 4))")
        form_str <- paste(form_str, paste(cont_terms, collapse = " + "), sep = " + ")
      }

      final_formula <- as.formula(form_str)
      
      train_data$Y <- train_data$pse_Y_cond
      
      # Fit the gamm model with random intercept by ID and REML smoothing
      mod <- mgcv::gam(
        formula = final_formula,
        data = train_data,
        family = gaussian,
        method = "REML"
      )
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
  
  if (model == "Pseudo outcome - One - DR"){
    if (method == "Super learner"){
      mod <- SuperLearner(Y = train_data$pse_Y, X = data.frame(subset(train_data, select = covariates)),
                          method = "method.NNLS",
                          family = gaussian(),
                          cvControl = list(V = 10, stratifyCV=FALSE),
                          SL.library = SL_lib,
                          obsWeights = train_data$weights_DR)
    }
  }
  if (model == "Pseudo outcome - One - R"){
    if (method == "Super learner"){
      mod <- SuperLearner(Y = train_data$pse_Y, X = data.frame(subset(train_data, select = covariates)),
                          method = "method.NNLS",
                          family = gaussian(),
                          cvControl = list(V = 10, stratifyCV=FALSE),
                          SL.library = SL_lib,
                          obsWeights = train_data$weights_R)
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

          if (learner == "surviTMLE-learner" | learner == "M-learner"){
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
            if (learner == "surviTMLE-learner" | learner == "M-learner"){
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

          if (learner == "surviTMLE-learner" | learner == "M-learner"){
          keep_covs <- append(c("ID","Y","A","C","tstart","time","s",
                                "S_k_pred_0","S_k_pred_1","H_k_pred_0","H_k_pred_1","h_k_pred_0","h_k_pred_1","at_risk"),covariates)
          }
          else if (learner == "T-learner"){
            keep_covs <- append(c("ID","Y","tstart","time",
                                  "S_k_pred_0","S_k_pred_1"),covariates)
          }
          pred_data_long_all <- subset(pred_data_long_all,select = keep_covs)
        }
        if (method == "Local survival stack" | method == "Global survival stack"){
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
          S_k_pred_long_all_0 <- pred_data_long_all$S_k_pred_0
          S_k_pred_long_all_1 <- pred_data_long_all$S_k_pred_1

          if (learner == "surviTMLE-learner"){
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
  
  
  if (model == "Outcome" & LT == 1){
   tryCatch(
     {
       #--- Obtaining preds from outcome models ---#
       if (method == "Parametric"){
         #Survival estimates at the time point (Training and all data points)
         S_k_pred_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="survival")
         S_k_pred_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="survival")

         if (learner == "surviTMLE-learner"){
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
       if (method == "Local survival stack" | method == "Global survival stack"){
         uni_time <- rep(evt_times_uni, times = nrow(pred_data))

         #Obtaining predictions from outcome models
         S_k_pred_wide_all_0 <- as.data.frame(mod_0$S_T_preds)
         S_k_pred_wide_all_1 <- as.data.frame(mod_1$S_T_preds)

         if (learner != "ltrc" | is.na(learner) == TRUE){
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
           S_k_pred_long_all_0 <- pred_data_long_all$S_k_pred_0
           S_k_pred_long_all_1 <- pred_data_long_all$S_k_pred_1

           if (learner == "surviTMLE-learner"){
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
       }
       if (method == "pcox"){
         uni_time <- rep(evt_times_uni, times = nrow(pred_data))
         
         #Obtaining predictions from outcome models
         S_k_pred_wide_all_0 <- as.data.frame(fit.pred_0$surv)
         S_k_pred_wide_all_1 <- as.data.frame(fit.pred_1$surv)
         
         if (learner != "ltrc" | is.na(learner) == TRUE){
           #Making into long form
           ncols_temp <- ncol(pred_data)
           S_k_pred_wide_all_0 <- as.data.frame(t(S_k_pred_wide_all_0))
           S_k_pred_wide_all_1 <- as.data.frame(t(S_k_pred_wide_all_1))
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
           S_k_pred_long_all_0 <- pred_data_long_all$S_k_pred_0
           S_k_pred_long_all_1 <- pred_data_long_all$S_k_pred_1

           # if (learner == "surviTMLE-learner"){
           #   #Calculating cumulative hazard estimates
           #   pred_data_long_all$H_k_pred_0 <- -log(pred_data_long_all$S_k_pred_0)
           #   pred_data_long_all$H_k_pred_1 <- -log(pred_data_long_all$S_k_pred_1)
           # 
           #   #Calculating hazards at each time
           #   pred_data_long_all <- pred_data_long_all %>%
           #     arrange(ID, time) %>%
           #     group_by(ID) %>%
           #     mutate(h_k_pred_0 = c(H_k_pred_0[1], diff(H_k_pred_0)))
           #   pred_data_long_all <- pred_data_long_all %>%
           #     arrange(ID, time) %>%
           #     group_by(ID) %>%
           #     mutate(h_k_pred_1 = c(H_k_pred_1[1], diff(H_k_pred_1)))
           # 
           #   #Survival predictions for time point of interest
           #   temp_data_list <- vector("list", length(evt_times_uni))
           #   preds0 <- cbind(ID=pred_data$ID,S_k_pred_wide_all_0)
           #   preds1 <- cbind(ID=pred_data$ID,S_k_pred_wide_all_1)
           #   for (t in seq_along(evt_times_uni)){
           #     #Collecting appropriate survival preds
           #     pred0 <- preds0[,c(1,t+1)]
           #     names(pred0)[2] <- "S_t_pred_long_all_0"
           #     pred1 <- preds1[,c(1,t+1)]
           #     names(pred1)[2] <- "S_t_pred_long_all_1"
           # 
           #     #Merging with long pred data
           #     pred_data_long_all_temp <- merge(pred_data_long_all,pred0,by="ID",all.x=TRUE)
           #     pred_data_long_all_temp <- merge(pred_data_long_all_temp,pred1,by="ID",all.x=TRUE)
           # 
           #     pred_data_long_all_temp <- pred_data_long_all_temp %>% mutate(S_t_pred_long_all_a = case_when(
           #       A == 1 ~ S_t_pred_long_all_1,
           #       A == 0 ~ S_t_pred_long_all_0
           #     ))
           # 
           #     # Saving data to list
           #     temp_data_list[[t]] <- pred_data_long_all_temp
           #   }
           # }
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
  
  if (model == "Outcome - Diff" & LT == 1){
    tryCatch(
      {
        #--- Obtaining preds from outcome models ---#
        if (method == "Local survival stack" | method == "Global survival stack"){
          #Obtaining predictions from outcome models
          G_diff_k_pred_wide_all <- as.data.frame(mod$S_T_preds)
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Truncation weight" & LT == 1){
    tryCatch(
      {
        #--- Obtaining preds from outcome models ---#
        if (method == "Local survival stack" | method == "Global survival stack"){
          #Need a survival prediction at each persons truncation time

          #Collecting preds from model
          preds_wide_all <- as.data.frame(mod$S_T_preds)

          #Identifying peoples pred at their truncation time
          train_data$index <- match(train_data$Q, trunc_times)
          train_data$trunc_weight_pred <- 9999
          for (i in 1:dim(train_data)[1]){
            train_data$trunc_weight_pred[i] <- preds_wide_all[i,train_data$index[i]]
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
        if (method == "Local survival stack" | method == "Global survival stack"){
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
  
  if (model == "Censoring" & LT == 1){
    tryCatch(
      {
        #--- Obtaining preds from censoring model ---#
        if (method == "Parametric"){   #Not updated to include G or obtain correct estimates 
          #Censoring estimates (Training and all data points)
          G_k_pred_long_all <- predict(mod,newdata=pred_data_long_all,type="survival")
        }
        else if (method == "Local survival stack" | method == "Global survival stack"){
          uni_time <- rep(evt_times_uni, times = nrow(pred_data_long_all))

          #Obtaining predictions from outcome models
          G_k_pred_wide_all <- as.data.frame(mod$S_T_preds)
          
          #Making into long form
          ncols_temp <- ncol(pred_data_long_all)
          
          #Creating dataset with ID and preds (but they need transforming as time along columns)
          temp <- cbind(ID=pred_data_long_all$ID,G_k_pred_wide_all)
          
          #Splitting matrices for each person 
          matrix_list <- split(temp[, -1], temp$ID)  # Remove ID column, split by ID
          matrix_list <- lapply(matrix_list, as.matrix)

          #Transforming matrices
          transformed_list <- lapply(matrix_list, t)
          
          #Re-collating into a long dataset with preds
          n_rows <- sum(sapply(transformed_list, nrow))
          n_cols <- ncol(transformed_list[[1]])
          col_names <- paste0("G_k_Q", seq_len(n_cols))
          
          stacked_matrix <- as.data.frame(matrix(NA, nrow = n_rows, ncol = n_cols))
          colnames(stacked_matrix) <- col_names

          # Fill the pre-allocated data frame
          row_counter <- 1
          for (mat in transformed_list) {
            current_rows <- nrow(mat)
            stacked_matrix[row_counter:(row_counter + current_rows - 1), ] <- as.data.frame(mat)
            row_counter <- row_counter + current_rows
          }
          
          pred_data_long_all <- cbind(pred_data_long_all,stacked_matrix)
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
        else if (method == "Local survival stack"){
          uni_time <- rep(evt_times_uni, times = nrow(pred_data))

          #Obtaining predictions from truncation model
          trunc_k_pred_wide_all <- as.data.frame(mod$S_T_preds)

          if (learner != "ltrc"){     #May cause issues for surv-iTMLE
            #Making into long form
            ncols_temp <- ncol(pred_data)
            pred_data <- cbind(pred_data,trunc_k_pred_wide_all)
            
            pred_data_temp <- pred_data %>%
              pivot_longer(
                cols = starts_with("V"),
                names_to = "time2",
                values_to = "trunc_k_pred"
              )
            pred_data_temp$trunc_k_pred <- 1 - pred_data_temp$trunc_k_pred
            pred_data_long_all <- pred_data_temp
            pred_data_long_all$time <- uni_time
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

  if (model == "Pseudo outcome - Pooled - Factor" | model == "Pseudo outcome - Pooled - Continuous" | model == "Pseudo outcome - Pooled - CI"){
    pred_data_long_all <- pred_data_long_all[order(pred_data_long_all$ID,pred_data_long_all$time),]
    pred_data <- subset(pred_data_long_all, select = c(covariates,"time"))
    if (model == "Pseudo outcome - Pooled - Factor"){
      pred_data$time <- as.factor(pred_data$time)
    }
    else if (model == "Pseudo outcome - Pooled - Continuous"){
      pred_data$time <- as.numeric(pred_data$time)
    }
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
    else if (method == "GAM"){
      pred_data_long_all$pred_cond <- predict(mod, newdata = pred_data)
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
  
  if (model == "Pseudo outcome - One - DR" | model == "Pseudo outcome - One - R"){
    if (method == "Super learner"){
      pred_data <- subset(pred_data, select = c(covariates))
      pred <- predict(mod, pred_data)$pred
    }
  }



  #-----------------------------#
  #--- Returning information ---#    
  #-----------------------------#

  if (model == "Outcome" & (method == "Super learner" | method == "Global survival stack" | method == "Local survival stack") & (learner == "surviTMLE-learner" | learner == "M-learner")){
    output <- list(out_mod_0 = mod_0,
                  out_mod_1 = mod_1,
                  pred_data_long_all_pred = pred_data_long_all,
                  Surv_t_pred_long_all_list=temp_data_list)
  }
  if (model == "Outcome" & (method == "Super learner" | method == "Global survival stack" | method == "Local survival stack") & (learner == "ltrc")){
    output <- list(S_0_pred = S_k_pred_wide_all_0,
                   S_1_pred = S_k_pred_wide_all_1)
  }
  if (model == "Outcome - Diff"){
    output <- G_diff_k_pred_wide_all
  }
  else if (model == "Outcome" & (method != "Super learner" & method != "Local survival stack" & method != "Global survival stack")  & (learner == "surviTMLE-learner" | learner == "M-learner")){
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
  else if (model == "Outcome" & learner == "T-learner" & (method == "Super learner" | method == "Local survival stack" | method == "Global survival stack")){
    output <- list(pred_data_long_all_pred = pred_data_long_all,
                   out_mod_0 = mod_0,
                   out_mod_1 = mod_1)
  }
  else if (model == "Outcome" & learner == "T-learner" & method == "pcox"){
    output <- list(S_k_pred_long_all_0 = S_k_pred_long_all_0,
                   S_k_pred_long_all_1 = S_k_pred_long_all_1)
  }
  else if (model == "Truncation weight"){
    output <- list(trunc_weight_pred = train_data$trunc_weight_pred)
  }
  else if (model == "Propensity score"){
    output <- list(e_pred = pred,
                   e_mod = mod)
  }
  else if (model == "Censoring" & method == "Parametric"){
    output <- list(G_k_pred_long_all = G_k_pred_long_all,
                   G_mod = mod)
  }
  else if (model == "Censoring" & (method == "Super learner" | method == "Local survival stack" | method == "Global survival stack")){
    output <- list(G_k_pred_long_all = pred_data_long_all,
                   G_mod = mod)
  }
  else if (model == "Truncation" & method == "Parametric" & learner != "ltrc"){
    output <- list(trunc_k_pred_long_all = trunc_k_pred_long_all,
                   trunc_mod = mod)
  }
  else if (model == "Truncation" & method == "Local survival stack" & learner != "ltrc"){
    output <- list(trunc_k_pred_long_all = pred_data_long_all,
                   trunc_mod = mod)
  }
  else if (model == "Truncation" & method == "Local survival stack" & learner == "ltrc"){
    output <- trunc_k_pred_wide_all
  }
  else if (model == "Pseudo outcome - Pooled - Factor" | model == "Pseudo outcome - Pooled - Continuous"){
    output <- list(po_mod = mod,
                   pred = pred_data_long_all)
  }
  else if (model == "Pseudo outcome - Pooled - CI"){
    output <- list(pred = pred_data_long_all$pred)
  }
  else if (model == "Pseudo outcome - One - DR" | model == "Pseudo outcome - One - R"){
    output <- pred
  }

  return(output)
}


