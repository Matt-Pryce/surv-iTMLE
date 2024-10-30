#######################################################################################
# Script: Function to simulate survival data. 
# Date: 18/10/24
# Author: Matt Pryce 
# Notes:  
#######################################################################################

# library(tidyverse)
# library(stringr)
# library(panelr)
# library(lava)
# library(MASS)
# library(reshape2)
# library(data.table)


#########################################################################
#--- Uni-variate outcome with missingness - Binary treatment setting ---#  
#########################################################################

#' @param n Number of observations to be simulated
#' @param x The number of covariates
#' @param x_form The distribution of the covariates; choose from either ("Unif","Binary")
#' @param e The function which generates the propensity score. Must take arguments from X1 - Xk, where k is number of covariates
#' @param Y_0 The function which generates the baseline outcomes for the covariate set w
#' @param Tau The function which generates the CATE for the covariate set w
#' @param G The function which generates the probability of not having a missing outcome. 
#' @param error The natural variation introduced in the unexposed outcome function. 
#' @param outcome_form Whether the outcome is binary or continuous 
# 
#' @return A dataset containing n observations 


data_gen_surv <- function(n = 1000,
                          x = 10,
                          x_form = c("Unif","Norm"),
                          e_func,
                          Y_type,
                          H_0_lambda_Y,
                          H_0_gamma_Y,
                          cox_func_Y,
                          aft_func_Y,
                          pois_func_Y,
                          unif_end_Y,
                          C_type,
                          H_0_lambda_C,
                          H_0_gamma_C,
                          cox_func_C,
                          aft_func_C,
                          pois_func_C,
                          unif_end_C,
                          right_cen = NULL)
  {
  ID <- rep(1:n)
  
  #Generate baseline covariates
  if (x_form == "Unif"){
    if (x>1){
      X <- data.frame(matrix(runif(n*x,0,1), n, x))
      X <- data.frame(cbind(ID,X))
    }
    if (x==1){
      X1 <- runif(n,0,1)
      X <- data.frame(cbind(ID,X1))
    }
  }
  if (x_form == "Norm"){
    tot_num_covs <- x
    Sigma <- diag(x = 1, nrow = tot_num_covs, ncol = tot_num_covs)
    offdiag(Sigma, type = 0) <- 0
    X <- data.frame(mvrnorm(n = n, mu = rep(x = 0, times = tot_num_covs), Sigma,empirical = T))
    X <- data.frame(cbind(ID,X))
  }
  
  
  #Generate the propensity score
  prop_score <- e_func(X)
  
  #Generate the treatment
  A <- rbinom(n, size=1, prob= prop_score)
  X <- cbind(X,A)
  
  if (Y_type == "Cox"){
    #Cox needs to specify baseline hazard and combination
    U <- runif(n)
    T <- (-log(U) / (H_0_lambda_Y * exp(cox_func_Y(X))))^(1 / H_0_gamma_Y)
  }
  if (Y_type == "AFT"){
    epsilon <- rnorm(dim(X)[1])
    log_T <- aft_func_Y(X) + epsilon
    T <- exp(log_T)
  }
  if (Y_type == "Poisson"){
    T <- rep(-99,n)
    event_counts <- rpois(n, lambda = pois_func_Y(X))
    for (i in 1:n){
      if (event_counts[i] > 0){
        T[i] <- rexp(1, rate = event_counts[i])
      }
      else {
        T[i] <- Inf
      }
    }
  }
  if (Y_type == "Unif"){
    T <- runif(n,0,unif_end_Y)
  }

  if (C_type == "Cox"){
    #Cox needs to specify baseline hazard and combination
    U <- runif(n)
    C <- (-log(U) / (H_0_lambda_C * exp(cox_func_C(X))))^(1 / H_0_gamma_C)
  }
  if (C_type == "AFT"){
    epsilon <- rnorm(dim(X)[1])
    log_C <- aft_func_C(X) + epsilon
    C <- exp(log_C)
  }
  if (C_type == "Poisson"){
    C <- rep(-99,n)
    cen_counts <- rpois(n, lambda = pois_func_C(X))
    for (i in 1:n){
      if (cen_counts[i] > 0){
        C[i] <- rexp(1, rate = cen_counts[i])
      }
      else {
        C[i] <- Inf
      }
    }
  }
  if (C_type == "Unif"){
    C <- runif(n,0,unif_end_C)
  }

  #Determine observed times (minimum of survival and censoring times)
  T_tilde <- pmin(T, C)

  if (is.null(right_cen) == 0){
    T_tilde <- pmin(T_tilde,right_cen)
  }

  # Indicator for censoring (1 if the event occurred, 0 if censored)
  Y <- as.numeric((T <= C) & (T_tilde != right_cen))
  delta <- as.numeric((T <= C))
  cen_ind <- 1-Y


  
  # #Collating dataset
  sim_data <- data.frame(X,prop_score,T,C,T_tilde,Y,delta,cen_ind)
  return(sim_data)
}

#####################################################################
# 
# e_func <- function(X){
#   X2 <- X$X2; 
#   X3 <- X$X3;
#   e_test <- plogis(pmax(X2,X3))
#   e_test
# }
# 
# cox_func <- function(X){
#   X2 <- X$X2; 
#   X3 <- X$X3;
#   cox_test <- pmax(X2,X3)
#   cox_test
# }
# 
# aft_func <- function(X){
#   X2 <- X$X2; 
#   X3 <- X$X3;
#   epsilon <- rnorm(dim(X)[1])
#   aft_test <- pmax(X2,X3) + epsilon 
#   aft_test
# }
# 
# pois_func <- function(X){
#   X2 <- X$X2; 
#   X3 <- X$X3;
#   pois_test <- X2+X3 
#   pois_test
# }
# 
# H_0_lambda <- 1
# H_0_gamma <- 2
# 
# 
# n=400
# 
# covs <- 6
# 
# sim_data <- data_gen_surv(n = n,
#                           x = covs,
#                           x_form = "Unif",
#                           e_func = e_func,
#                           Y_type = "Poisson",
#                           cox_func_Y = cox_func,
#                           aft_func_Y = aft_func,
#                           pois_func_Y = pois_func,
#                           unif_end_Y = 3,
#                           H_0_lambda_Y = H_0_lambda,
#                           H_0_gamma_Y = H_0_gamma,
#                           C_type = "Unif",
#                           cox_func_C = cox_func,
#                           aft_func_C = aft_func,
#                           pois_func_C = pois_func,
#                           H_0_lambda_C = H_0_lambda,
#                           H_0_gamma_C = H_0_gamma,
#                           unif_end_C = 3,
#                           right_cen = 1)
# 

# #Finding truth 
# 
# X <- data.frame(X2=0.4,X3=0.7)
# 
# set.seed(1)
# store <- rep(-99,1000000)
# for (i in 1:1000000){
#   U <- runif(1)
#   T <- (-log(U) / (H_0_lambda * exp(cox_func(X))))^(1 / H_0_gamma)
#   store[i] <- T
# }
# time <- 0.5
# store_time <- (time > store)
# prob <- 1 - sum(store_time)/1000000

#Smarten up and make more effcient 

#Cox with weibull baseline 
#lambda*t^{gamma} is the form of the cumulative baseline hazard
#t^2 has baseline hazard of 2t

#To find true survival prob, generate a survival time for an individual loads of times, what proportion of them times are higher than ... . , however , generate extremely large dataset with 
#How many times did they have an event before that time 
