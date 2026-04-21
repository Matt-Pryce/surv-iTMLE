######################################################################
# Script: Simulations script for Specification 1
# Author: Matt Pryce
# Date: 12/08/24
# Notes:
######################################################################

# #--- Loading libraries needed ---# 

library(base)
library(survival)
library(timereg)
library(boot)
library(grf)
library(SuperLearner)
library(tidyverse)
library(riskRegression)
library(pec)
library(timeROC)
library(xgboost)
library(ranger)
library(KernelKnn)
library(nnet)
library(e1071)
library(randomForest)
library(Sieve)
library(survSuperLearner)
library(pch)
library(glmnet)
library(dplyr)
library(survML)
library(splines)
library(gbm)
library(cvTools)  # used for creating folds for cross-fitting
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("fast_integrals.cpp")


#--- Loading parameter info ---#
load("parameters_log.RData")

#--- Defining parameter values for the specification ---#
scenario <- 1
setting <- 1
n  <- parameters_log$sample_size[setting]
covs  <- parameters_log$num_covs[setting]
ps_spec  <- parameters_log$ps_list[setting]
cen_spec  <- parameters_log$cen_list[setting]
cen_type_spec  <- parameters_log$cen_type_list[setting]
cen_H0_lambda_spec  <- parameters_log$cen_H0_lambda_list[setting]
cen_H0_gamma_spec  <- parameters_log$cen_H0_gamma_list[setting]
cen_unif_spec <- parameters_log$cen_unif_list[setting]
Y_spec  <- parameters_log$Y_list[setting]
Y_type_spec  <- parameters_log$Y_type_list[setting]
Y_H0_lambda_spec  <- parameters_log$Y_H0_lambda_list[setting]
Y_H0_gamma_spec  <- parameters_log$Y_H0_gamma_list[setting]
Y_unif_spec <- parameters_log$Y_unif_list[setting]
Q_type_spec <- parameters_log$Q_type_list[setting]
Q_alpha_base_spec <- parameters_log$Q_alpha_base_list[setting]
Q_alpha_func_spec <- parameters_log$Q_alpha_func_list[setting]
Q_beta_base_spec <- parameters_log$Q_beta_base_list[setting]
Q_beta_func_spec <- parameters_log$Q_beta_func_list[setting]
Q_max_spec <- parameters_log$Q_max_list[setting]
right_cen_spec <- parameters_log$right_cen_list[setting]

if (scenario == 1){
  time_seq <- c(0.1,seq(from=0.25,to=2,by=0.25))
}
if (scenario == 2){
  time_seq <- seq(from=0.1,to=1,by=0.1)
}
if (scenario == 3){
  time_seq <- seq(from=0.2,to=2,by=0.2)
}

#--- Loading and defining functions to generate data ---#
ps_list <- readRDS("ps_funcs.RData")
ps_func <- ps_list[[ps_spec]]

cen_types <- readRDS("cen_types.RData")
cen_type <- cen_types[[cen_type_spec]]

cen_H0_lambda_list <- readRDS("cen_H0_lambda_list.RData")
cen_H0_lambda <- cen_H0_lambda_list[[cen_H0_lambda_spec]]

cen_H0_gamma_list <- readRDS("cen_H0_gamma_list.RData")
cen_H0_gamma <- cen_H0_gamma_list[[cen_H0_gamma_spec]]

cen_funcs <- readRDS("cen_funcs.RData")
cen_func <- cen_funcs[[cen_spec]]

Y_funcs <- readRDS("Y_funcs.RData")
Y_func <- Y_funcs[[Y_spec]]

Y_types <- readRDS("Y_types.RData")
Y_type <- Y_types[[Y_type_spec]]

Y_H0_lambda_list <- readRDS("Y_H0_lambda_list.RData")
Y_H0_lambda <- Y_H0_lambda_list[[Y_H0_lambda_spec]]

Y_H0_gamma_list <- readRDS("Y_H0_gamma_list.RData")
Y_H0_gamma <- Y_H0_gamma_list[[Y_H0_gamma_spec]]

Q_types <- readRDS("Q_types.RData")
Q_type <- Q_types[[Q_type_spec]]

Q_alpha_base_num_list <- readRDS("Q_alpha_base_list.RData")
Q_alpha_base <- Q_alpha_base_num_list[[Q_alpha_base_spec]]

Q_beta_base_num_list <- readRDS("Q_beta_base_list.RData")
Q_beta_base <- Q_beta_base_num_list[[Q_beta_base_spec]]

Q_alpha_func_num_list <- readRDS("Q_alpha_func_list.RData")
Q_alpha_func <- Q_alpha_func_num_list[[Q_alpha_func_spec]]

Q_beta_func_num_list <- readRDS("Q_beta_func_list.RData")
Q_beta_func <- Q_beta_func_num_list[[Q_beta_func_spec]]

Q_max_num_list <- readRDS("Q_max_list.RData")
Q_max <- Q_max_num_list[[Q_max_spec]]

right_cen_list <- readRDS("right_cen_time_list.RData")
right_cen <- right_cen_list[[right_cen_spec]]


#--- Loading test data ---#
test_data_name <- paste("spec_",scenario,"_test_data.RData",sep="")
load(test_data_name)


#--- Creating learners for SL library's ---#
#LASSO & elastic net
nlambda_seq = c(50,100,250)
alpha_seq <- c(0.5,1)
usemin_seq <- c(FALSE,TRUE)
para_learners = create.Learner("SL.glmnet", tune = list(nlambda = nlambda_seq,alpha = alpha_seq,useMin = usemin_seq))

#Random forest
mtry_seq <-  c(3,5)
min_node_seq <- c(10,20)
num_trees_seq <- c(500)
sample_fraction_seq <- c(0.2,0.4,0.6)
rf_learners = create.Learner("SL.ranger", tune = list(mtry = mtry_seq,
                                                      min.node.size = min_node_seq,
                                                      num.trees=num_trees_seq,
                                                      sample.fraction=sample_fraction_seq))



#--- Loading data generation and learner functions ---#
#Data generation
source("Data_gen_surv.R")

# #Scripts needed for learners
source("Data_management_surv.R")
source("Nuisance_models_surv.R")
source("ltrc_all_extra_funcs_2.R")

#Learner
source("T_learner.R")
source("CSF.R")
source("survEP_learner.R")
source("ltrc_estimators.R")


#--- Running loop ---#

# Loop: 1) Set seed
#       2) Generate data
#       3) Run models & store information


start <- Sys.time()

#Setting up temporary list to store model info for the simulation in the spec
model_info_list <- list()

#Number of simulations
sims <- 1

seed <- taskID  #Relates to HPC use, would need defining if run on a computer

#Starting loop
for (i in 1:sims){
  
  #Setting up temporary list to store model info for the simulation in the script
  model_info_sim_list <- list()
  
  #Identifying scenario
  scen <- setting
  
  #-----------------------#
  #--- Simulating data ---#
  #-----------------------#
  #Generate dataset
  sim_data_train <- data_gen_surv(n = n,
                                  x = covs,
                                  x_form = "Unif",
                                  e_func = ps_func,
                                  T_type = Y_type,
                                  cox_func_T = Y_func,
                                  aft_func_T = Y_func,
                                  pois_func_T = Y_func,
                                  H_0_lambda_T = Y_H0_lambda,
                                  H_0_gamma_T = Y_H0_gamma,
                                  C_type = cen_type,
                                  cox_func_C = cen_func,
                                  aft_func_C = cen_func,
                                  pois_func_C = cen_func,
                                  H_0_lambda_C = cen_H0_lambda,
                                  H_0_gamma_C = cen_H0_gamma,
                                  Q_type = Q_type,
                                  Q_alpha_base = Q_alpha_base,
                                  Q_alpha_func = Q_alpha_func,
                                  Q_beta_base = Q_beta_base,
                                  Q_beta_func = Q_beta_func,
                                  max_Q = Q_max,
                                  right_cen = right_cen,
                                  LT = 1)
  
  model_info_sim_list <- append(model_info_sim_list,list(sim_data_train = sim_data_train))
  model_info_sim_list <- append(model_info_sim_list,list(sim_data_test = sim_data_test))
  
  #Creating covariate lists
  cov_list <- NULL
  for (cov in 1:covs){
    cov_list <- append(cov_list,paste("X",cov,sep=""))
  }
  
  
  #-----------------------------------------------------------------------#
  #--- Running T learner - Local survival stacking - Accounting for LT ---#
  #-----------------------------------------------------------------------#
  tryCatch(
    {
      event.SL.library <- c("SL.mean","SL.glm",
                            "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                            "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                            "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                            "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                            "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                            "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")
      
      #Trimming data to contain only observed individuals
      sim_data_train <- subset(sim_data_train,sim_data_train$Observed == 1)
      
      T_learner_SS_LT <- T_learner(data = sim_data_train,
                                   estimand = "Difference",
                                   id = "ID",
                                   time = "T_tilde",
                                   outcome = "delta",
                                   censor = "cen_ind",
                                   exposure = "A",
                                   truncation = "Q",
                                   time_cuts = time_seq,
                                   out_covariates = cov_list,
                                   out_method = "Local survival stack",
                                   out_SL_lib = event.SL.library,
                                   newdata = sim_data_test)
      
      model_info_sim_list <- append(model_info_sim_list,list(T_SS_LT = T_learner_SS_LT))
    },
    error=function(e) {
      message(paste("An error occured when fitting the T-learner SS2 for scenario ",scen," with sample size ",n," in simulation ",seed,sep=""))
      print(e)
    }
  )
  
  
  #-----------------------------------#
  #--- Running CSF - Complete case ---#
  #-----------------------------------#
  tryCatch(
    {
      #Trimming data to contain only observed individuals
      sim_data_train <- subset(sim_data_train,sim_data_train$Observed == 1)

      CSF_mod <- CSF(data = sim_data_train,
                     id = "ID",
                     time = "T_tilde",
                     outcome = "delta",
                     censor = "cen_ind",
                     exposure = "A",
                     time_cuts = time_seq,
                     splits = 10,
                     prop_est = FALSE,
                     covariates = cov_list,
                     e_method = "Super learner",
                     newdata = sim_data_test)

      model_info_sim_list <- append(model_info_sim_list,list(CSF = CSF_mod))
    },
    error=function(e) {
      message(paste("An error occured when fitting the CSF for scenario ",scen," with sample size ",n," in simulation ",seed,sep=""))
      print(e)
    }
  )


  #----------------------------------------------------------------------------------------#
  #--- Running survEP - Target 3 - No Iso - Local survival stacking - Accounting for LT ---#
  #----------------------------------------------------------------------------------------#
  tryCatch(
    {
      pse_lib <- c("SL.mean","SL.lm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")

      bin_lib <- c("SL.mean","SL.glm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")

      #Trimming data to contain only observed individuals
      sim_data_train <- subset(sim_data_train,sim_data_train$Observed == 1)

      survEP_mod3_SS <- survEP_learner(data = sim_data_train,
                                       id = "ID",
                                       time = "T_tilde",
                                       outcome = "delta",
                                       censor = "cen_ind",
                                       exposure = "A",
                                       truncation = "Q",
                                       time_cuts = time_seq,
                                       splits = 10,
                                       e_covariates = cov_list,
                                       e_method = "Super learner",
                                       e_SL_lib = bin_lib,
                                       out_covariates = cov_list,
                                       out_method = "Local survival stack",
                                       out_SL_lib = bin_lib,
                                       g_covariates = cov_list,
                                       g_method = "Local survival stack",
                                       g_SL_lib = bin_lib,
                                       h_covariates = cov_list,
                                       h_method = "Local survival stack",
                                       h_SL_lib = bin_lib,
                                       iso_reg = FALSE,
                                       pse_covariates = cov_list,
                                       pse_approach = "Pooled - Factor",
                                       pse_method = "Super learner",
                                       pse_SL_lib = pse_lib,
                                       newdata = sim_data_test,
                                       target_option = "Lasso - Linear - Option 3",
                                       CI = FALSE)


      model_info_sim_list <- append(model_info_sim_list,list(survEP3_SS = survEP_mod3_SS))
    },
    error=function(e) {
      message(paste("An error occured when fitting the survEP3 SS for scenario ",scen," with sample size ",n," in simulation ",seed,sep=""))
      print(e)
    }
  )
  
  
  
  #--------------------------------------------------------#
  #--- Running ltrc estimator - Local survival stacking ---#
  #--------------------------------------------------------#
  tryCatch(
    {
      sim_data_train <- subset(sim_data_train,sim_data_train$Observed == 1)
      
      trim = 0.05
      trim.est = 0
      
      options.F = list(trim = trim.est,
                       ntree = 500, mtry = 2,
                       df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5)
      options.Sd = list(trim = trim.est,
                        ntree = 500, mtry = 2,
                        df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                        nfolds.OOF = 5)
      options.G = list(trim = trim.est,
                       df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                       trunc.weighting = TRUE)
      options.PS = list(trim = trim.est,
                        df = 7,
                        ntree = 500)
      
      num_search_rounds = 5
      ntrees_max = 500
      
      pse_lib <- c("SL.mean","SL.lm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")
      
      bin_lib <- c("SL.mean","SL.glm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")
      
      
      ests <- trunclearner(sim_data_train,
                           nu = nu,
                           X.name = "T_tilde",    #Observed time
                           Q.name = "Q",          #Truncation time
                           event.name = "delta",  #Event indicator
                           A.name = "A",          #Treatment
                           K = 10,                #Number of folds for cross-fitting of nuisance params
                           model.T = "Local survival stack",
                           model.D = "Local survival stack",
                           model.Q = "Local survival stack",
                           model.A = "Super learner", 
                           cov.names.T = cov_list,
                           cov.names.binary.T = NULL,
                           cov.names.D = cov_list,
                           cov.names.binary.D = NULL,
                           cov.names.Q = cov_list,
                           cov.names.binary.Q = NULL,
                           cov.names.A = cov_list,
                           cov.names.binary.A = NULL,
                           cov.names.CATE = cov_list,
                           cov.names.CATE.binary = NULL,
                           est_approach_G = "truncIPW.F",
                           est_approach_PS = "truncIPW.F",
                           options.F = options.F,
                           options.Sd = options.Sd,
                           options.G = options.G,
                           options.PS = options.PS,
                           num_search_rounds=num_search_rounds, ntrees_max=ntrees_max,
                           e_SL_lib = bin_lib,
                           out_SL_lib = bin_lib,
                           h_SL_lib = bin_lib,
                           pse_SL_lib = pse_lib,
                           pse_option = "boost",
                           newdata = sim_data_test,
                           trim = 0.05,
                           time_grid = time_seq)
      
      DR_ests <- ests$est_DR %>% as.data.frame()
      colnames(DR_ests) <- time_seq
      DR_ests$ID <- c(1:10000)
      DR_ests <- pivot_longer(
        DR_ests,
        cols = -ID,
        names_to = "time",
        values_to = "est"
      )
      
      R_ests <- ests$est_R %>% as.data.frame()
      colnames(R_ests) <- time_seq
      R_ests$ID <- c(1:10000)
      R_ests <- pivot_longer(
        R_ests,
        cols = -ID,
        names_to = "time",
        values_to = "est"
      )
      
      model_info_sim_list <- append(model_info_sim_list,list(ltrc_DR_SS = DR_ests))
      model_info_sim_list <- append(model_info_sim_list,list(ltrc_R_SS = R_ests))
    },
    error=function(e) {
      message(paste("An error occured when fitting the ltrc SS for scenario ",scen," with sample size ",n," in simulation ",seed,sep=""))
      print(e)
    }
  )
  
  
  
  #------------------------------------------------------------------------------------------#
  #--- Running ltrc estimator - Local survival stacking - Their nuisance function choices ---#
  #------------------------------------------------------------------------------------------#
  tryCatch(
    {
      sim_data_train <- subset(sim_data_train,sim_data_train$Observed == 1)
      
      trim = 0.05
      trim.est = 0
      
      options.F = list(trim = trim.est,
                       ntree = 500, mtry = 2,
                       df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5)
      options.Sd = list(trim = trim.est,
                        ntree = 500, mtry = 2,
                        df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                        nfolds.OOF = 5)
      options.G = list(trim = trim.est,
                       df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                       trunc.weighting = TRUE)
      options.PS = list(trim = trim.est,
                        df = 7,
                        ntree = 500)
      
      num_search_rounds = 5
      ntrees_max = 500
      
      pse_lib <- c("SL.mean","SL.lm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")
      
      bin_lib <- c("SL.mean","SL.glm",
                   "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                   "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                   "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                   "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
                   "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                   "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")
      
      
      ests <- trunclearner(sim_data_train,
                           nu = nu,
                           X.name = "T_tilde",    #Observed time
                           Q.name = "Q",          #Truncation time
                           event.name = "delta",  #Event indicator
                           A.name = "A",          #Treatment
                           K = 10,                #Number of folds for cross-fitting of nuisance params
                           model.T = "pCox",
                           model.D = "pCox",
                           model.Q = "pCox", 
                           model.A = "gbm", 
                           cov.names.T = cov_list,
                           cov.names.binary.T = NULL,
                           cov.names.D = cov_list,
                           cov.names.binary.D = NULL,
                           cov.names.Q = cov_list,
                           cov.names.binary.Q = NULL,
                           cov.names.A = cov_list,
                           cov.names.binary.A = NULL,
                           cov.names.CATE = cov_list,
                           cov.names.CATE.binary = NULL,
                           est_approach_G = "truncIPW.F",
                           est_approach_PS = "truncIPW.F",
                           options.F = options.F,
                           options.Sd = options.Sd,
                           options.G = options.G,
                           options.PS = options.PS,
                           num_search_rounds=num_search_rounds, ntrees_max=ntrees_max,
                           e_SL_lib = bin_lib,
                           out_SL_lib = bin_lib,
                           h_SL_lib = bin_lib,
                           pse_SL_lib = pse_lib,
                           pse_option = "boost",
                           newdata = sim_data_test,
                           trim = 0.2,
                           time_grid = time_seq)
      
      DR_ests_pcox <- ests$est_DR %>% as.data.frame()
      colnames(DR_ests) <- time_seq
      DR_ests$ID <- c(1:10000)
      DR_ests <- pivot_longer(
        DR_ests,
        cols = -ID,
        names_to = "time",
        values_to = "est"
      )
      
      R_ests_pcox <- ests$est_R %>% as.data.frame()
      colnames(R_ests) <- time_seq
      R_ests$ID <- c(1:10000)
      R_ests <- pivot_longer(
        R_ests,
        cols = -ID,
        names_to = "time",
        values_to = "est"
      )
      
      model_info_sim_list <- append(model_info_sim_list,list(ltrc_DR_pcox = DR_ests_pcox))
      model_info_sim_list <- append(model_info_sim_list,list(ltrc_R_pcox = R_ests_pcox))
    },
    error=function(e) {
      message(paste("An error occured when fitting the ltrc for scenario ",scen," with sample size ",n," in simulation ",seed,sep=""))
      print(e)
    }
  )
  
  model_info_list <- append(model_info_list,list(i = model_info_sim_list))
}


end <- Sys.time()
end-start



