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
                          pred_data_long,
                          pred_data_long_all,
                          evt_times_uni
                          #SL_lib,
){
  
  #------------------------------#
  #--- Creating training data ---#
  #------------------------------#
  
  tryCatch(
    {
      if (model == "Outcome"){
        train_data0 <- subset(data,data$A==0)
        train_data0 <- subset(train_data0,select = c("time","Y",covariates))
        train_data1 <- subset(data,data$A==1)
        train_data1 <- subset(train_data1,select = c("time","Y",covariates))
      } 
      else if (model == "Propensity score"){
        train_data <- as.data.frame(subset(data,select = c(covariates,"A","s")))
      }
      else if (model == "Censoring"){
        train_data <- as.data.frame(subset(data,select = c("time","C",covariates,"A","s")))
      }
      # else if (model == "Pseudo outcome"){
      #   train_data <- data
      # }
    },
    error=function(e) {
      stop('An error occured when creating analysis data')
      print(e)
    }
  )
  
  #----------------------#
  #--- Running models ---#
  #----------------------#
  
  if (model == "Outcome"){
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          mod_0 <- coxph(Surv(time, Y) ~ ., data = train_data0)
          mod_1 <- coxph(Surv(time, Y) ~ ., data = train_data1)
        }
        # if (method == "Random forest"){
        #   X_S_0 <- as.matrix(subset(train_data0, select = covariates))
        #   mod_0 <- regression_forest(X_S_0, train_data0$Y, honesty = FALSE,tune.parameters = "all")
        #
        #   X_S_1 <- as.matrix(subset(train_data1, select = covariates))
        #   mod_1 <- regression_forest(X_S_1, train_data1$Y, honesty = FALSE,tune.parameters = "all")
        # }
        # else
        # else if (method == "Super learner"){
        #   if (Y_bin == 1){
        #     mod_0 <- SuperLearner(Y = train_data0$Y, X = data.frame(subset(train_data0, select = covariates)),
        #                           method = "method.NNLS",
        #                           family = binomial(),
        #                           cvControl = list(V = 10, stratifyCV=TRUE),
        #                           SL.library = SL_lib)
        #     mod_1 <- SuperLearner(Y = train_data1$Y, X = data.frame(subset(train_data1, select = covariates)),
        #                           method = "method.NNLS",
        #                           family = binomial(),
        #                           cvControl = list(V = 10, stratifyCV=TRUE),
        #                           SL.library = SL_lib)
        #   }
      },
      error=function(e) {
        stop('An error occured when running outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Censoring"){
    tryCatch(
      {
        if (method == "Parametric"){
          #Running first outcome models
          mod <- coxph(Surv(time, C) ~ ., data = train_data)
        }
        # if (method == "Random forest"){
        #   X_S_0 <- as.matrix(subset(train_data0, select = covariates))
        #   mod_0 <- regression_forest(X_S_0, train_data0$Y, honesty = FALSE,tune.parameters = "all")
        #
        #   X_S_1 <- as.matrix(subset(train_data1, select = covariates))
        #   mod_1 <- regression_forest(X_S_1, train_data1$Y, honesty = FALSE,tune.parameters = "all")
        # }
        # else
        # else if (method == "Super learner"){
        #   if (Y_bin == 1){
        #     mod_0 <- SuperLearner(Y = train_data0$Y, X = data.frame(subset(train_data0, select = covariates)),
        #                           method = "method.NNLS",
        #                           family = binomial(),
        #                           cvControl = list(V = 10, stratifyCV=TRUE),
        #                           SL.library = SL_lib)
        #     mod_1 <- SuperLearner(Y = train_data1$Y, X = data.frame(subset(train_data1, select = covariates)),
        #                           method = "method.NNLS",
        #                           family = binomial(),
        #                           cvControl = list(V = 10, stratifyCV=TRUE),
        #                           SL.library = SL_lib)
        #   }
      },
      error=function(e) {
        stop('An error occured when running censoring model')
        print(e)
      }
    )
  }
  
  if (model == "Propensity score"){  # | model == "Censoring" | model == "Pseudo outcome" | model == "IPCW"){
    if (method == "Random forest"){
      if (model == "Propensity score"){
        X <- as.matrix(subset(train_data, select = covariates))
        mod <- regression_forest(X, train_data$A, honesty = FALSE,tune.parameters = "all")
      }
      # else if (model == "Pseudo outcome"){
      #   X <- as.matrix(subset(train_data, select = covariates))
      #   mod <- regression_forest(X, train_data$pse_Y)
      # }
    }
    else if (method == "Parametric"){
      if (model== "Propensity score"){
        fit_data <- subset(train_data,select = -c(s))
        mod <- glm(A ~ . , data = fit_data, family = binomial())
      }
      # else if (model == "Pseudo outcome"){
      #   fit_data <- subset(train_data,select = c("pse_Y",covariates))
      #   mod <- lm(pse_Y ~ . , data = fit_data)
      # }
    }
    else if (method == "Super learner"){
      if (model == "Propensity score"){
        sums <- table(train_data$A)
        cv_folds <- min(10,sums[1],sums[2])
        mod <- SuperLearner(Y = train_data$A, X = data.frame(subset(train_data, select = covariates)),
                            method = "method.NNLS",
                            family = binomial(),
                            cvControl = list(V = cv_folds, stratifyCV=TRUE),
                            SL.library = SL_lib)
      }
      # else if (model == "Pseudo outcome"){
      #   mod <- SuperLearner(Y = train_data$pse_Y, X = data.frame(subset(train_data, select = covariates)),
      #                       method = "method.NNLS",
      #                       family = gaussian(),
      #                       cvControl = list(V = 10, stratifyCV=FALSE),
      #                       SL.library = SL_lib)
      # }
    }
  }
  
  


  #-----------------------------------#
  #--- Obtaining model predictions ---#
  #-----------------------------------#

  if (model == "Outcome"){
    tryCatch(
      {
        #--- Obtaining preds from outcome models ---#
        if (method == "Parametric"){
          #Survival estimates at the time point (Training and all data points)
          S_k_pred_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="survival")
          S_k_pred_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="survival")

          #Hazard estimates at the time point (All time point data)
          pred_data_long_all$H_k_long_all_0 <- predict(mod_0,newdata=pred_data_long_all,type="expected")
          pred_data_long_all <- pred_data_long_all %>%
            arrange(ID, tstart) %>%
            group_by(ID) %>%
            mutate(h_k_long_all_0 = c(0,diff(H_k_long_all_0)))

          H_k_pred_long_all_0 <- pred_data_long_all$H_k_long_all_0
          h_k_pred_long_all_0 <- pred_data_long_all$h_k_long_all_0

          pred_data_long_all$H_k_long_all_1 <- predict(mod_1,newdata=pred_data_long_all,type="expected")
          pred_data_long_all <- pred_data_long_all %>%
            arrange(ID, tstart) %>%
            group_by(ID) %>%
            mutate(h_k_long_all_1 = c(0,diff(H_k_long_all_1)))

          H_k_pred_long_all_1 <- pred_data_long_all$H_k_long_all_1
          h_k_pred_long_all_1 <- pred_data_long_all$h_k_long_all_1


          #Survival predictions for time point of interest
          temp_data_list <- vector("list", length(evt_times_uni))  #2)
          for (t in seq_along(evt_times_uni)){  #[2]
            temp_data <- pred_data_long_all
            temp_data$time2 <- temp_data$time      #Need time to be set to a number but don't want to lose original value
            temp_data$time <- evt_times_uni[t]

            temp_data_list[[t]] <- temp_data
          }

          # S_t_pred_long_all_0_list <- vector("list", 2)
          # S_t_pred_long_all_1_list <- vector("list", 2)
          # S_t_pred_long_all_0_list[[2]] <- predict(mod_0,newdata=temp_data_list[[2]],type="survival")
          # S_t_pred_long_all_1_list[[2]] <- predict(mod_1,newdata=temp_data_list[[2]],type="survival")
          # #UNCOMMENT AND REINSTATE WHEN MAKING PREDICTIONS FOR ALL TIME POINTS
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
        # if (method == "Random forest"){
        #   mod_pred_0 <- predict(mod_0, as.matrix(subset(pred_data,select=c(covariates))))$pred
        #   mod_pred_1 <- predict(mod_1, as.matrix(subset(pred_data,select=c(covariates))))$pred
        # }
        # else
        # else if (method == "Super learner"){
        #   pred_data <- subset(pred_data,select = c(covariates))
        #   mod_pred_0 <- predict(mod_0, pred_data)$pred
        #   mod_pred_1 <- predict(mod_1, pred_data)$pred
        # }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }

  if (model == "Censoring"){
    tryCatch(
      {
        #--- Obtaining preds from censoring model ---#
        if (method == "Parametric"){
          #Censoring estimates (Training and all data points)
          G_k_pred_long_all <- predict(mod,newdata=pred_data_long_all,type="survival")
        }
      },
      #if an error occurs, tell me the error
      error=function(e) {
        stop('An error occured when drawning predictions from outcome models')
        print(e)
      }
    )
  }
  
  if (model == "Propensity score"){   # | model == "Censoring"){
    if (method == "Random forest"){
      pred <- predict(mod, pred_data)$predictions
    }
    else if (method == "Parametric"){
      pred <- predict(mod, as.data.frame(pred_data), type = "response")
    }
    else if (method == "Super learner"){
      pred <- predict(mod, as.data.frame(pred_data))$pred
    }
  }

  # if (model == "Pseudo outcome"){
  #   if (method == "Random forest"){
  #     pred_data_matrix <- as.matrix(subset(pred_data, select = c(covariates)))
  #     pred <- predict(mod, pred_data_matrix)
  #   }
  #   else if (method == "Parametric"){
  #     pred_data <- subset(pred_data, select = c(covariates))
  #     pred <- as.data.frame(predict(mod, pred_data, type = "response"))
  #   }
  #   else if (method == "Super learner"){
  #     pred_data <- subset(pred_data, select = c(covariates))
  #     pred <- predict(mod, pred_data)$pred
  #   }
  # }



  #-----------------------------#
  #--- Returning information ---#
  #-----------------------------#

  if (model == "Imputation" | model == "IPCW"){
    output <- analysis_data
  }
  else if (model == "Outcome"){
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
  else if (model == "Propensity score"){
    output <- list(e_pred = pred,
                   e_mod = mod)
  }
  else if (model == "Censoring"){
    output <- list(G_k_pred_long_all = G_k_pred_long_all,
                   G_mod = mod)
  }
  else if (model == "Pseudo outcome"){
    output <- list(po_pred = pred,
                   po_mod = mod)
  }

  return(output)
}


