library(glmnet)
library(survival)
library(xgboost)
library(SuperLearner)
library(dplyr)
library(splines)
library(gbm)
library(cvTools)  # used for creating folds for cross-fitting
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("C:/Users/MatthewPryce/OneDrive - London School of Hygiene and Tropical Medicine/Documents/PhD/iTMLE-surv/Github_repo/SurvEP-learner/src/fast_integrals.cpp")



trunclearner <- function(dat, nu, cov.names.CATE, cov.names.CATE.binary = NULL, K,
                         X.name, Q.name, event.name, A.name, trim = 1e-7,
                         model.T = NULL, model.Q = NULL, model.A = NULL, model.D = NULL, 
                         cov.names.T = NULL, cov.names.Q = NULL, cov.names.A = NULL, cov.names.D = NULL, 
                         cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                         options.F = NULL, options.G = NULL, 
                         options.PS = NULL, options.Sd = NULL,
                         est_approach_G = "truncIPW.F", est_approach_PS = "truncIPW.F", 
                         Fuz.mx = NULL, Fuz_A1.mx = NULL, Fuz_A0.mx = NULL,
                         Gvz.mx = NULL, Sdz.mx = NULL, PS = NULL,
                         df = NULL, nfolds = 10, # alpha = NULL, lambda = NULL,
                         metrics = "rmse",
                         booster = "gbtree",
                         k_folds=NULL,
                         objective= "reg:squarederror",
                         ntrees_max=500,
                         num_search_rounds=20,
                         print_every_n=100,
                         early_stopping_rounds=10,
                         nthread=NULL,
                         verbose=FALSE, 
                         simulation = TRUE,
                         trim1 = 0.05,
                         trim2 = 0,
                         save_cvfit = FALSE,
                         e_SL_lib,
                         out_SL_lib,
                         h_SL_lib,
                         pse_SL_lib,
                         pse_option,
                         newdata,
                         time_grid){
  # Remove the df and alpha option in the above imputs of the function
  # lambda also seems to be un-used
  
  
  dat_A1 = dat
  dat_A1[,A.name] <- 1
  dat_A0 = dat
  dat_A0[,A.name] <- 0
  
  
  # if(is.null(Gvz.mx)){
  #   names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
  #   if(sum(names == "delta.1")){
  #     stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
  #   }
  #   dat$delta.1 = rep(1, nrow(dat))
  #   event.name.Q = "delta.1"
  # }


  n = nrow(dat)
  event.name.T = event.name

  jumps.X = sort(dat[,X.name])                   #Event time jumps
  jumps.Q = sort(dat[,Q.name])                   #Truncation jumps
  jumps.Y = sort(dat[,X.name] - dat[,Q.name])    #Difference jumps 

  tau1 = min(jumps.X)             #Finding min event/censoring time
  tau2 = max(jumps.Q)             #Finding max truncation time
  tau3 = max(jumps.Y)             #Finding max difference time 
  tau.max = max(jumps.X)          #Finding max event/censoring time
  
  # Compute \V{\nu} and \V{1} using out of fold estimate of the nuisance parameters
  truncC_AIPW_1 = replicate(length(time_grid), rep(NA, n), simplify = FALSE)  # \V{1}
  truncC_AIPW_nuT = replicate(length(time_grid), rep(NA, n), simplify = FALSE)  # \V{\nu(T)}
  # y_R = rep(NA,n)  # [\V{\nu(T)\} - m(Z)] / [\V{1}(A-\pi(Z))]
  # PS_weights = rep(NA,n)  # {A-\pi(Z)}^2

  # Store some intermediate results that are useful for constructing IPW and naive estimators
  mZ_hat = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  mu0_hat = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  mu1_hat = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  GXZ_hat = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  SdXZ_hat = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  PS_hat = rep(NA, n)

  folds <- cvFolds(n, K)
  
  #dat$s <- rep(1:length(dat[,X.name]),1) %% K
  for(k in 1:K){
  #for(k in 0:(K-1)){
    id.est = folds$subsets[folds$which == k]
    id.fit = folds$subsets[folds$which != k]
    dat.est = dat[id.est, ]
    dat.fit = dat[id.fit, ]
    dat_A1.est = dat_A1[id.est, ]
    dat_A0.est = dat_A0[id.est, ]
    
    ### Estimate the nuisance parameters

    # Estimate F and F_A1 F_A0
    if (model.T == "pCox"){
      if(is.null(Fuz.mx)){
        u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)   #All times

        Fuz.mx.est = est_F(dat.fit, dat.est, u, model.T, X.name, Q.name, event.name,
                           cov.names = cov.names.T, cov.names.binary = cov.names.binary.T,
                           mtry = options.F$mtry, ntree = options.F$ntree,
                           nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                           lambda = options.F$lambda, df = options.F$df)
      }

      if(is.null(Fuz_A1.mx)){
       Fuz_A1.mx.est = est_F(dat.fit, dat_A1.est, u, model.T, X.name, Q.name, event.name,
                             cov.names = cov.names.T, cov.names.binary = cov.names.binary.T,
                             trim = options.F$trim,
                             mtry = options.F$mtry, ntree = options.F$ntree,
                             nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                             lambda = options.F$lambda, df = options.F$df)
      }

      if(is.null(Fuz_A0.mx)){
       Fuz_A0.mx.est = est_F(dat.fit, dat_A0.est, u, model.T, X.name, Q.name, event.name,
                             cov.names = cov.names.T, cov.names.binary = cov.names.binary.T,
                             mtry = options.F$mtry, ntree = options.F$ntree,
                             nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                             lambda = options.F$lambda, df = options.F$df)
      }
    }
    else if (model.T == "Local survival stack"){
      names(dat.fit)[names(dat.fit) == A.name] <- "A"
      names(dat.est)[names(dat.est) == A.name] <- "A"
      names(dat.fit)[names(dat.fit) == X.name] <- "time"
      names(dat.est)[names(dat.est) == X.name] <- "time"
      names(dat.fit)[names(dat.fit) == event.name] <- "Y"
      names(dat.est)[names(dat.est) == event.name] <- "Y"
      names(dat.fit)[names(dat.fit) == Q.name] <- "Q"
      names(dat.est)[names(dat.est) == Q.name] <- "Q"

      #dat.est <- dat.est %>% arrange(ID)

      newdata_long_all <- dat.est

      u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)   #All times

      outcome_models <- nuis_mod_surv(model = "Outcome",
                                      data = dat.fit,
                                      method = model.T,
                                      covariates = cov.names.T,
                                      pred_data_long_all = newdata_long_all,
                                      evt_times_uni = u,
                                      SL_lib = out_SL_lib,
                                      learner = "ltrc",
                                      LT = 1)

      Fuz_A0.mx.est = 1 - outcome_models$S_0_pred
      colnames(Fuz_A0.mx.est) <- u
      Fuz_A0.mx.est <- as.matrix(Fuz_A0.mx.est)

      Fuz_A1.mx.est = 1 - outcome_models$S_1_pred
      colnames(Fuz_A1.mx.est) <- u
      Fuz_A1.mx.est <- as.matrix(Fuz_A1.mx.est)

      #Creating Fuz.mx.est
      for (id_order in 1:dim(Fuz_A0.mx.est)[1]){
        if (id_order == 1 & dat.est$A[id_order] == 0){
          Fuz.mx.est <- Fuz_A0.mx.est[1,]
        }
        else if (id_order == 1 & dat.est$A[id_order] == 1){
          Fuz.mx.est <- Fuz_A1.mx.est[1,]
        }
        else if (id_order > 1 & dat.est$A[id_order] == 0){
          Fuz.mx.est <- rbind(Fuz.mx.est,Fuz_A0.mx.est[id_order,])
        }
        else if (id_order > 1 & dat.est$A[id_order] == 1){
          Fuz.mx.est <- rbind(Fuz.mx.est,Fuz_A1.mx.est[id_order,])
        }
      }
      Fuz.mx.est <- as.matrix(Fuz.mx.est)
    }


    ## Estimate S_D
    if(model.D == "pCox"){
      if(is.null(Sdz.mx)){
        d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)

        Sdz.mx.est = est_Sd(dat.fit, dat.est, d, model.D, X.name, Q.name, event.name,
                            cov.names = cov.names.D, cov.names.binary = cov.names.binary.D,
                            trim = options.Sd$trim,
                            mtry = options.Sd$mtry, ntree = options.Sd$ntree,
                            nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                            lambda = options.Sd$lambda, df = options.Sd$df)
      }
    }
    else if (model.D == "Local survival stack"){
      #dat.est <- dat.est %>% arrange(ID)

      newdata_long_all <- dat.est

      d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)

      Sdz.mx.est <- nuis_mod_surv(model = "Outcome - Diff",
                                      data = dat.fit,
                                      method = model.D,
                                      covariates = cov.names.D,
                                      pred_data_long_all = newdata_long_all,
                                      evt_times_uni = d,
                                      SL_lib = out_SL_lib,
                                      learner = "ltrc",
                                      LT = 1)

      Sdz.mx.est <- as.matrix(Sdz.mx.est)
      # Sdz.mx.est <- 1 - Sdz.mx.est
      colnames(Sdz.mx.est) <- d
      Sdz.mx.est <- as.matrix(Sdz.mx.est)
    }



    ## Estimate G - using truncation weights 1/{1-F(Q|A,Z)} estimated using fit.si, and use fit.sj to estimate G, i\neq j \in \{1,2\}
    v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
    if(model.Q == "pCox"){
     if(est_approach_G == "truncIPW.F"){
       # dat.fit2 <- dat.fit
       # dat.est2 <- dat.est
       # names(dat.fit2)[names(dat.fit2) == A.name] <- "A"
       # names(dat.est2)[names(dat.est2) == A.name] <- "A"
       # names(dat.fit2)[names(dat.fit2) == X.name] <- "time"
       # names(dat.est2)[names(dat.est2) == X.name] <- "time"
       # names(dat.fit2)[names(dat.fit2) == event.name] <- "Y"
       # names(dat.est2)[names(dat.est2) == event.name] <- "Y"
       # names(dat.fit2)[names(dat.fit2) == Q.name] <- "Q"
       # names(dat.est2)[names(dat.est2) == Q.name] <- "Q"
       #
       # # Compute the truncation weights 1/{1-F(Q|A,Z)}
       # Fq.vec.fit <- nuis_mod_surv(model = "Truncation weight",
       #                             data = dat.fit2,
       #                             method = "Local survival stack",
       #                             covariates = cov.names.T,
       #                             pred_data_long_all = dat.fit2,
       #                             evt_times_uni = dat.fit[,Q.name],
       #                             SL_lib = out_SL_lib,
       #                             learner = "ltrc",
       #                             LT = 1)
       #
       # w_truncF.fit = 1/pmax(Fq.vec.fit$trunc_weight_pred, trim)

       if(is.null(Fuz.mx)){
         Fq.vec.fit = diag(est_F(dat.fit, dat.fit, dat.fit[,Q.name],
                                 model.T, X.name, Q.name, event.name,
                                 cov.names = cov.names.T, cov.names.binary = cov.names.binary.T,
                                 mtry = options.F$mtry, ntree = options.F$ntree,
                                 nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                 lambda = options.F$lambda, df = options.F$df))
       }

       w_truncF.fit = 1/pmax(1-Fq.vec.fit, trim)

       dat.fit$delta.1 = rep(1, nrow(dat.fit))
       event.name.Q = "delta.1"

       # Estimate G using fit.s2
       Gvz.mx.est = est_G(dat.fit, dat.est, v, model.Q, X.name, Q.name, event.name.Q,
                          cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                          weights = w_truncF.fit, trunc = FALSE,
                          tau = tau.max, trim = options.G$trim,
                          nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                          lambda = options.G$lambda, df = options.G$df)
     }
    }
    else if (model.Q == "Local survival stack"){
      Fq.vec.fit = NULL

      # Compute the truncation weights 1/{1-F(Q|A,Z)} estimated using fit.s1
      Fq.vec.fit <- nuis_mod_surv(model = "Truncation weight",
                                  data = dat.fit,
                                  method = model.T,
                                  covariates = cov.names.T,
                                  pred_data_long_all = dat.fit,
                                  evt_times_uni = dat.fit[,Q.name],
                                  SL_lib = out_SL_lib,
                                  learner = "ltrc",
                                  LT = 1)

      dat.fit$trunc_weight_pred <- 1/pmax(Fq.vec.fit$trunc_weight_pred, trim)

      dat.fit$delta.1 = rep(1, nrow(dat.fit))

      Gvz.mx.est <- nuis_mod_surv(model = "Truncation",
                                  data = dat.fit,
                                  method = model.Q,
                                  covariates = cov.names.Q,
                                  evt_times_uni = v,
                                  SL_lib = h_SL_lib,
                                  pred_data_long_all = dat.est,
                                  learner = "ltrc",
                                  LT = 1)
      Gvz.mx.est <- as.matrix(Gvz.mx.est)
      Gvz.mx.est <- 1 - Gvz.mx.est
      colnames(Gvz.mx.est) <- v
      Gvz.mx.est <- as.matrix(Gvz.mx.est)
    }



    ## Estimate PS
    if(is.null(PS)){

     if(est_approach_PS == "truncIPW.F" & model.A == "gbm"){

       if(is.null(Fq.vec.fit)){ # The truncation weights from F haven't been computed yet

         # Compute the truncation weights 1/{1-F(Q|A,Z)} estimated using fit.s1
         if(is.null(Fuz.mx)){
           Fq.vec.fit = diag(est_F(dat.fit, dat.fit, dat.fit[,Q.name],
                                   model.T, X.name, Q.name, event.name.T,
                                   cov.names = cov.names.T, cov.names.binary = cov.names.binary.T,
                                   trim = options.F$trim,
                                   mtry = options.F$mtry, ntree = options.F$ntree,
                                   formula.survPen = options.F$formula.survPen,
                                   nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                   lambda = options.F$lambda, df = options.F$df))
         }

         w_truncF.fit = 1/pmax(1-Fq.vec.fit, trim)

       }

       # Estimate PS using truncation weights 1/(1-F(Q|A,Z))
       PS.est = est_PS(dat.fit, dat.est, model.A, A.name,
                       cov.names = cov.names.A, cov.names.binary = cov.names.binary.A,
                       weights = w_truncF.fit,
                       trim = options.PS$trim, df = options.PS$df,
                       ntree = options.PS$ntree)   # a vector of estimated PS for each individual
     }

     if(est_approach_PS == "truncIPW.F" & model.A == "Super learner"){

       dat.fit$s <- 0   #Added to prevent error in nuis_mod_surv script

       # Compute the truncation weights 1/{1-F(Q|A,Z)} estimated using fit.s1
       Fq.vec.fit <- nuis_mod_surv(model = "Truncation weight",
                                   data = dat.fit,
                                   method = model.T,
                                   covariates = cov.names.T,
                                   pred_data_long_all = dat.fit,
                                   evt_times_uni = dat.fit[,Q.name],
                                   SL_lib = out_SL_lib,
                                   learner = "ltrc",
                                   LT = 1)

       dat.fit$trunc_weight_pred <- 1/pmax(Fq.vec.fit$trunc_weight_pred, trim)

       PS_model <- nuis_mod_surv(model = "Propensity score",
                                 data = dat.fit,
                                 method = model.A,
                                 covariates = cov.names.A,
                                 SL_lib = e_SL_lib,
                                 pred_data = dat.est,
                                 LT = 1)

       PS.est <- PS_model$e_pred
     }
    }


    nuuT_A1.mx_list <- vector("list", length(time_grid))
    nuuT_A0.mx_list <- vector("list", length(time_grid))
    mu_A1.est_list <-  vector("list", length(time_grid))
    mu_A0.est_list <-  vector("list", length(time_grid))
    mu_Z.est_list <-  vector("list", length(time_grid))
    truncC_AIPW_result_list <- vector("list", length(time_grid))
    tau.Tmax_A1 = max(u) + 1e-5
    tau.Tmax_A0 = max(u) + 1e-5
    for (time in 1:length(time_grid)){

      t00 = time_grid[time]

      nu <- function(t,t0=t00){
        # indicator function
        result = as.numeric(t>t0)

        return(result)
      }

      ## Compute the mu_1(Z), mu_0(Z) for the treatment augmentation, and m(Z) invoved in the R-learner
      nuuT_A1.mx_list[[time]] = matrix(rep(nu(u), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
      nuuT_A0.mx_list[[time]] = matrix(rep(nu(u), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
      mu_A1.est_list[[time]] = as.vector(int_fmx_dF_cpp(tau.Tmax_A1, nuuT_A1.mx_list[[time]], Fuz_A1.mx.est, u))
      mu_A0.est_list[[time]] = as.vector(int_fmx_dF_cpp(tau.Tmax_A0, nuuT_A0.mx_list[[time]], Fuz_A0.mx.est, u))
      mu_Z.est_list[[time]] = mu_A1.est_list[[time]] * PS.est + mu_A0.est_list[[time]] * (1-PS.est)

      ## Compute \V{1}, y_R, and (A-\pi(Z))^2
      if (model.D == "Local survival stack"){
        truncC_AIPW_result_list[[time]] <- truncC_AIPW_transMean(dat.est, nu, Fuz.mx.est, Gvz.mx.est, Sdz.mx.est,
                                                    "time", "Q", "Y", trim)
      }
      if (model.D == "pCox"){
        truncC_AIPW_result_list[[time]] <- truncC_AIPW_transMean(dat.est, nu, Fuz.mx.est, Gvz.mx.est, Sdz.mx.est,
                                                    X.name, Q.name, event.name, trim)
      }

      truncC_AIPW_1[[time]][id.est] = truncC_AIPW_result_list[[time]]$truncC_const1
      truncC_AIPW_nuT[[time]][id.est] = truncC_AIPW_result_list[[time]]$truncC_nu

      # Record some intermediate quantities
      mZ_hat[[time]][id.est] = mu_Z.est_list[[time]]
      mu1_hat[[time]][id.est] = mu_A1.est_list[[time]]
      mu0_hat[[time]][id.est] = mu_A0.est_list[[time]]
      PS_hat[id.est] = PS.est
      if (model.Q == "Local survival stack"){
        GXZ_hat[[time]][id.est] = CDF_eval(dat.est[,"time"], Gvz.mx.est)
      }
      if (model.Q == "pCox"){
        GXZ_hat[[time]][id.est] = CDF_eval(dat.est[,X.name], Gvz.mx.est)
      }
      if (model.D == "Local survival stack"){
        SdXZ_hat[[time]][id.est] = 1 - CDF_eval(dat.est[,"time"] - dat.est[,"Q"], 1-Sdz.mx.est)
      }
      if (model.D == "pCox"){
        SdXZ_hat[[time]][id.est] = 1 - CDF_eval(dat.est[,X.name] - dat.est[,Q.name], 1-Sdz.mx.est)
      }
    }
  }

  A = dat[,A.name]

  nuT_tilde = replicate(length(time_grid), rep(NA, n), simplify = FALSE)
  for (time in 1:length(time_grid)){
    nuT_tilde[[time]] = truncC_AIPW_nuT[[time]]/bound_away_zero(truncC_AIPW_1[[time]], trim1)
  }


  ############### Estimate CATE ################

  #CATE_est_R <- replicate(length(time_grid), rep(NA, dim(newdata)[1]), simplify = FALSE)
  #CATE_est_DR <- replicate(length(time_grid), rep(NA, dim(newdata)[1]), simplify = FALSE)
  #for (time in 1:length(time_grid)){
  #  cov.names.CATE = c(cov.names.CATE, cov.names.CATE.binary)
  #  Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])
  #  Z_CATE_new = as.matrix(newdata[,cov.names.CATE, drop = FALSE])
  #
  #  if (pse_option == "boost"){
  #    # truncR-learner with xgboost
  #    PS_weights = (bound_away_zero(A-PS_hat, 0))^2  # used to be trim1
  #    weights_R = bound_away_zero(bound_away_zero(truncC_AIPW_1[[time]], 0)*PS_weights, trim2)  # used to be trim1
  #    y_R = (nuT_tilde[[time]] - mZ_hat[[time]]) / bound_away_zero(A-PS_hat, trim1)
  #
  #    cvfit_wsq_R  = cvboost_wsq(Z_CATE, y_R, weights = weights_R,
  #                               # k_folds = k_folds,
  #                               ntrees_max = ntrees_max,
  #                               num_search_rounds = num_search_rounds)#,
  #    # print_every_n = print_every_n,
  #    # nthread  = nthread,
  #    # verbose = verbose)
  #    CATE_est_R[[time]] = predict(cvfit_wsq_R, newx = Z_CATE_new)
  #
  #    # truncDR-learner with xgboost
  #    mu_A.vec = A*mu1_hat[[time]] + (1-A)*mu0_hat[[time]]
  #    y_DR = (A-PS_hat)/(pmax(PS_hat, trim)*pmax(1-PS_hat, trim)) * (nuT_tilde[[time]] - mu_A.vec) + mu1_hat[[time]] - mu0_hat[[time]]
  #    weights_DR = bound_away_zero(truncC_AIPW_1[[time]], 0)  # used to be trim1
  #    cvfit_wsq_DR  = cvboost_wsq(Z_CATE, y_DR, weights = weights_DR,
  #                               # k_folds = k_folds,
  #                               ntrees_max = ntrees_max,
  #                               num_search_rounds = num_search_rounds)#,
  #                               # print_every_n = print_every_n,
  #                               # nthread  = nthread,
  #                               # verbose = verbose)
  #    CATE_est_DR[[time]] = predict(cvfit_wsq_DR, newx = Z_CATE_new)
  #  }
  #  else if (pse_option == "Super learner"){
  #    #--- DR ---#
  #    mu_A.vec = A*mu1_hat + (1-A)*mu0_hat
  #    y_DR = (A-PS_hat)/(pmax(PS_hat, trim)*pmax(1-PS_hat, trim)) * (nuT_tilde - mu_A.vec) + mu1_hat - mu0_hat
  #    weights_DR = bound_away_zero(truncC_AIPW_1, 0)
  #    DR_data <- data.frame(Z_CATE,
  #                          mu_A.vec=mu_A.vec,
  #                          pse_Y=y_DR,
  #                          weights_DR=weights_DR)
  #
  #    CATE_est_DR <- nuis_mod_surv(model = "Pseudo outcome - One - DR",
  #                                  data = DR_data,
  #                                  method = pse_option,
  #                                  covariates = cov.names.CATE,
  #                                  SL_lib = pse_SL_lib,
  #                                  pred_data = newdata)
  #
  #    #--- R ---#
  #    PS_weights = (bound_away_zero(A-PS_hat, 0))^2  # used to be trim1
  #    weights_R = bound_away_zero(bound_away_zero(truncC_AIPW_1, 0)*PS_weights, trim2)  # used to be trim1
  #    y_R = (nuT_tilde - mZ_hat) / bound_away_zero(A-PS_hat, trim1)
  #    R_data <- data.frame(Z_CATE,
  #                         pse_Y=y_R,
  #                         weights_R=weights_R)
  #
  #    CATE_est_R <- nuis_mod_surv(model = "Pseudo outcome - One - R",
  #                                  data = R_data,
  #                                  method = pse_option,
  #                                  covariates = cov.names.CATE,
  #                                  SL_lib = pse_SL_lib,
  #                                  pred_data = newdata)
  #  }
  #
  #}
  #
  #result_list = list(est_DR = CATE_est_DR,
  #                   est_R = CATE_est_R)
  result_list = list(truncC_AIPW_result_list=truncC_AIPW_result_list,
                     GXZ_hat=GXZ_hat,
                     SdXZ_hat=SdXZ_hat)
  return(result_list)
}


load("~/PhD/iTMLE-surv/Simulations/LTRC_simulations/Output/Model_results/Full_sims/Setting_1/setting_1_all_output.RData")
test_data <- data_list[[1]]$sim_data_train
test_data_all <- subset(test_data,test_data$Observed == 1)

#test_data <- sim_data_test
#test_data_all <- subset(test_data,test_data$Observed == 1)


t00_all = c(0.1,seq(from=0.25,to=2,by=0.25))

t00 = t00_all[9]

nu <- function(t,t0=t00){
  # indicator function
  result = as.numeric(t>t0)

  # result = pmin(t, t0)

  # result = log(t)

  return(result)
}


bound_away_zero <- function(x, trim){
  sign_x = sign(x)
  sign_x[x == 0] <- 1
  y = sign_x * pmax(abs(x), trim)

  return(y)
}


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

num_search_rounds = 1
ntrees_max = 100

event.SL.library2 <- c("SL.mean","SL.glm",
                       "SL.glmnet_1","SL.glmnet_2","SL.glmnet_3", "SL.glmnet_4",
                       # "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
                       # "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
                       "SL.ranger_1","SL.ranger_2")#,"SL.ranger_3","SL.ranger_4",
                       # "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
                       # "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")

pse_lib <- c("SL.mean","SL.lm")#,
             # "SL.glmnet_1","SL.glmnet_2")#,"SL.glmnet_3", "SL.glmnet_4",
             # "SL.glmnet_5","SL.glmnet_6","SL.glmnet_7", "SL.glmnet_8",
             # "SL.glmnet_8","SL.glmnet_9","SL.glmnet_11", "SL.glmnet_12",
             # "SL.ranger_1","SL.ranger_2","SL.ranger_3","SL.ranger_4",
             # "SL.ranger_5","SL.ranger_6","SL.ranger_7","SL.ranger_8",
             # "SL.ranger_9","SL.ranger_10","SL.ranger_11","SL.ranger_12")



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


#time_seq <- c(0.1,seq(from=0.25,to=2,by=0.25)) 
#for (i in 2:length(time_seq)){
#  time_seq_use <- time_seq[1:i]
#}

check <- trunclearner(test_data_all,
                      nu = nu,
                      X.name = "T_tilde",    #Observed time
                      Q.name = "Q",          #Truncation time
                      event.name = "delta",  #Event indicator
                      A.name = "A",          #Treatment
                      K = 10,                 #Number of folds for cross-fitting of nuisance params
                      model.T = "Local survival stack",#"pCox",#
                      model.D = "Local survival stack",#"pCox",#
                      model.Q = "Local survival stack",#"pCox", #
                      model.A = "Super learner", #"gbm",  #
                      cov.names.T = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                                      "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),
                      cov.names.binary.T = NULL,
                      cov.names.D = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                                      "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),
                      cov.names.binary.D = NULL,
                      cov.names.Q = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                                      "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),
                      cov.names.binary.Q = NULL,
                      cov.names.A = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                                      "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),
                      cov.names.binary.A = NULL,
                      cov.names.CATE = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                                         "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),
                      cov.names.CATE.binary = NULL,
                      est_approach_G = "truncIPW.F",
                      est_approach_PS = "truncIPW.F",
                      options.F = options.F,
                      options.Sd = options.Sd,
                      options.G = options.G,
                      options.PS = options.PS,
                      num_search_rounds=num_search_rounds, ntrees_max=ntrees_max,
                      e_SL_lib = event.SL.library2,
                      out_SL_lib = event.SL.library2,
                      h_SL_lib = event.SL.library2,
                      pse_SL_lib = pse_lib,
                      pse_option = "boost",
                      newdata = test_data,
                      trim = 0.05,
                      time_grid = c(0.1,seq(from=0.25,to=2,by=0.25)))



# #--- Trail for running for multiple times ---#
# time_seq <- c(0.05,0.1,seq(from=0.2,to=0.2,by=0.2))  #Need an extra time before the one you want
# DR_ests <- data.frame(matrix(NA, nrow = dim(test_data)[1], ncol = (length(time_seq)-1)))
# names(DR_ests) <- time_seq[2:length(time_seq)]
# R_ests <- data.frame(matrix(NA, nrow = dim(test_data)[1], ncol = (length(time_seq)-1)))
# names(R_ests) <- time_seq[2:length(time_seq)]
# for (i in 2:length(time_seq)){
#   time_seq_use <- time_seq[1:i]
# 
#   t00 = max(time_seq_use)
# 
#   nu <- function(t,t0=t00){
#     # indicator function
#     result = as.numeric(t>t0)
# 
#     # result = pmin(t, t0)
# 
#     # result = log(t)
# 
#     return(result)
#   }
# 
#   ests <- trunclearner(test_data_all2,
#                        nu = nu,
#                        X.name = "T_tilde",    #Observed time
#                        Q.name = "Q",          #Truncation time
#                        event.name = "delta",  #Event indicator
#                        A.name = "A",          #Treatment
#                        K = 2,                 #Number of folds for cross-fitting of nuisance params
#                        model.T = "Local survival stack",#"pCox", #"pCox",
#                        model.D = "Local survival stack",#"pCox", #"pCox",
#                        model.Q = "Local survival stack",#"pCox", ##"pCox",
#                        model.A = "Super learner", #"gbm",  #
#                        cov.names.T = c("X1","X2","X3","X4","X5"),
#                        cov.names.binary.T = NULL,
#                        cov.names.D = c("X1","X2","X3","X4","X5"),
#                        cov.names.binary.D = NULL,
#                        cov.names.Q = c("X1","X2","X3","X4","X5"),
#                        cov.names.binary.Q = NULL,
#                        cov.names.A = c("X1","X2","X3","X4","X5"),
#                        cov.names.binary.A = NULL,
#                        cov.names.CATE = c("X1","X2","X3","X4","X5"),
#                        cov.names.CATE.binary = NULL,
#                        est_approach_G = "truncIPW.F",
#                        est_approach_PS = "truncIPW.F",
#                        options.F = options.F,
#                        options.Sd = options.Sd,
#                        options.G = options.G,
#                        options.PS = options.PS,
#                        num_search_rounds=num_search_rounds, ntrees_max=ntrees_max,
#                        e_SL_lib = event.SL.library2,
#                        out_SL_lib = event.SL.library2,
#                        h_SL_lib = event.SL.library2,
#                        pse_SL_lib = pse_lib,
#                        pse_option = "boost",
#                        newdata = test_data,
#                        trim = 0.2)
# 
#   DR_ests[,i-1] <- ests$est_DR#[,1]
#   R_ests[,i-1] <- ests$est_R#[,1]
# 
# }


