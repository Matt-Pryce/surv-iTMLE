

###############################
###  Nuisance function est  ###
###############################

est_Sd <- function(dat.fit, dat.est = dat.fit, time.eval, model, 
                   time.name, Q.name, event.name, cov.names = NULL, trim = 0, 
                   mtry = NA, ntree = NA, formula.survPen = NA, 
                   x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                   nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, 
                   df = 5, OOF = FALSE, nfolds.OOF = 10){
  
  names = c(time.name, Q.name, event.name, cov.names)
  if(is.na(match("Y", names)) & is.na(match("delta.D", names))){
    dat.fit$Y = dat.fit[,time.name] - dat.fit[,Q.name]    # residual censored event time
    dat.fit$delta.D = 1- dat.fit[,event.name]    # censoring indicator
  }else{
    stop("The names of the variables cannot be 'Y' or 'delta.D'. Conflict with intermediate variables in the function.")
  }
  
  
 if(model == "pCox"){
    
    if(is.null(cov.names)){
      if(is.null(x.fit)|is.null(x.est)){
        stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
      }
    }else{
      dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])
      XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
      x.fit = XX[1:nrow(dat.fit), ]
      x.est = XX[-(1:nrow(dat.fit)), ]
    }
    
    yss = Surv(dat.fit$Y, dat.fit$delta.D)
    cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, alpha = alpha, lambda = lambda)
    fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est)
    Sdz.mx = t(fit.pred$surv)
    colnames(Sdz.mx) <- fit.pred$time
    
    Fdz.mx0 = 1 - Sdz.mx
    Sdz.mx = 1 - CDF_eval.mx(time.eval, Fdz.mx0)
    # # Visualize the estimates --------------------------------------------
    # nplot = 8
    # plot(d, Sdz.mx[1,], type = "l", col = 1)
    # for(i in 2:nplot){
    #     lines(d, Sdz.mx[i,], col = 4+i)
    # }
    # # End visualize -------------------------------------------------------
    
  }
  
  Sdz.mx = pmax(Sdz.mx, trim)
  return(Sdz.mx)
}






est_F <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name = NULL, event.name, cov.names, trim = 0,
                  weights = NULL, 
                  mtry = NULL, ntree = NULL, 
                  formula.survPen = NULL, 
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, 
                  df = 5){
  
  # if(is.null(trim)){ trim = 1e-7 }
  if(is.null(nfolds)){ nfolds = 10 }
  if(is.null(s)){ s = "lambda.min" }
  if(is.null(alpha)){ alpha = 1 }
  if(is.null(df)){ df = 5 }
  
  if((!is.null(weights)) & model == "RF"){
    stop("The function LTRCforests::ltrcrrf() does not take in case weights.")
  }
  if((!is.null(weights)) & model == "survPen"){
    stop("The function survPen::survPen() does not take in case weights.")
  }
  
  
  u = time.eval
  
  
  if(model == "pCox"){
    
    if(is.null(cov.names)){
      if(is.null(x.fit)|is.null(x.est)){
        stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
      }
    }else{
      dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])            
      XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
      x.fit = XX[1:nrow(dat.fit), ]
      x.est = XX[-(1:nrow(dat.fit)), ]
    }
    
    if(is.null(Q.name)){
      yss = Surv(dat.fit[,time.name], dat.fit[,event.name])
    }else{
      yss = Surv(dat.fit[,Q.name], dat.fit[,time.name], dat.fit[,event.name])
    }
    
    if(is.null(weights)){
      cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, 
                          alpha = alpha, lambda = lambda)
      fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est,
                                   alpha = alpha, lambda = lambda)
    }else{
      cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, weights = weights, 
                          alpha = alpha, lambda = lambda)
      fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est,
                                   weights = weights, alpha = alpha, lambda = lambda)
    }
    
    Fuz.mx0 = 1 - t(fit.pred$surv)
    colnames(Fuz.mx0) <- fit.pred$time
    
    Fuz.mx = CDF_eval.mx(time.eval, Fuz.mx0)
    # # Visualize the estimates --------------------------------------------
    # nplot = 8
    # plot(time.eval, Fuz.mx[1,], type = "l", col = 1)
    # for(i in 2:nplot){
    #     lines(time.eval, Fuz.mx[i,], col = i)
    # }
    # # End visualize ------------------------------------------------------
    
  }else{
    stop("This T model is not implemented in this function!")
  }
  Fuz.mx = pmin(Fuz.mx, 1-trim)

  return(Fuz.mx)
  # return(list(Fuz.mx = Fuz.mx, beta.T = beta.T))
}







est_G <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name, event.name, cov.names, trim = 0, weights = NULL,
                  formula.survPen = NULL, tau = NULL, trunc = TRUE,
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, df = 5){
  
  if((!is.null(weights)) & model == "RF"){
    stop("The function LTRCforests::ltrcrrf() does not take in case weights.")
  }
  if((!is.null(weights)) & model == "survPen"){
    stop("The function survPen::survPen() does not take in case weights.")
  }
  
  if(mean(dat.fit[,event.name])<1){
    stop("The truncation time is always observed, so dat.fit[,event.name] should be a vector of one's.")
  }
  
  names = c(time.name, Q.name, event.name, cov.names)
  if(sum(names == "Q2")){
    stop("The names of the variables cannot be 'Q2'. It is used in the middle of the computation.")
  }
  if(sum(names == "T2")){
    stop("The names of the variables cannot be 'T2'. It is used in the middle of the computation.")
  }
  
  v = time.eval
  
  if(is.null(tau)){
    tau = max(c(dat.fit[,time.name], dat.est[,time.name])) + 1
  }
  
  dat.fit$Q2 = tau - dat.fit[ ,Q.name]
  dat.fit$T2 = tau - dat.fit[ ,time.name]
  dat.fit$delta.1 = rep(1, nrow(dat.fit))
  
  dat.est$Q2 = tau - dat.est[ ,Q.name]
  dat.est$T2 = tau - dat.est[ ,time.name]
  dat.est$delta.1 = rep(1, nrow(dat.est))
  
  if(model == "pCox"){
    
    if(is.null(cov.names)){
      if(is.null(x.fit)|is.null(x.est)){
        stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
      }
    }else{
      dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])            
      XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
      x.fit = XX[1:nrow(dat.fit), ]
      x.est = XX[-(1:nrow(dat.fit)), ]
    }
    
    if(trunc){
      yss = Surv(dat.fit$T2, dat.fit$Q2, dat.fit$delta.1)
    }else{
      yss = Surv(dat.fit$Q2, dat.fit$delta.1)
    }
    
    
    if(is.null(weights)){
      cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, alpha = alpha, lambda = lambda)
    }else{
      cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, weights = weights, alpha = alpha, lambda = lambda)
    }
    
    fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est, 
                                 weights = weights, alpha = alpha, lambda = lambda)
    Gvz.mx0 = t(fit.pred$surv)
    colnames(Gvz.mx0) <- tau - fit.pred$time
    id0 = order(tau-fit.pred$time)
    Gvz.mx0 = Gvz.mx0[,id0]
    
    Gvz.mx = CDF_eval.mx(time.eval, Gvz.mx0)
  }
  
  colnames(Gvz.mx) = v
  Gvz.mx = pmax(Gvz.mx, trim)
  
  return(Gvz.mx)

}



est_PS <- function(dat.fit, dat.est, model, A.name, cov.names, cov.names.binary = NULL, 
                   weights = rep(1, nrow(dat.fit)), 
                   trim = 0, df = 7, ntree = 2000){
  if(model =="gbm"){
    if(is.null(ntree)){ ntree = 2000 }
    cov.names = c(cov.names, cov.names.binary)
    formula.A = formula(paste(A.name, "~", paste(cov.names, collapse = " + ")))
    psfit = gbm(formula.A, data = dat.fit, distribution = 'bernoulli', 
                weights = weights, n.trees = ntree)
    
    PS <- predict(psfit, newdata = dat.est, n.trees = ntree, type = 'response')
    
  }
  else{
    stop("This model for A is not implemented yet.")
  }
  
  PS = pmin(pmax(PS, trim), 1-trim)
  return(PS)
}




############################
###  AIPW and aug funcs  ###
############################


truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx,
                                   T.name, Q.name, trim = 1e-7, u = NULL, v = NULL) {
  
  if (nrow(dat) != nrow(Fuz.mx) || nrow(dat) != nrow(Gvz.mx)) {
    stop("The number of rows of dat, Fuz.mx, and Gvz.mx must be the same.")
  }
  
  nn <- nrow(dat)
  Ttime <- as.numeric(dat[, T.name])
  Q <- as.numeric(dat[, Q.name])
  
  if(is.null(u)){
    u <- as.numeric(colnames(Fuz.mx))  # jumps.T
  }
  if(is.null(v)){
    v <- as.numeric(colnames(Gvz.mx))  # jumps.Q
  }
  
  nuu <- nu(u)            # used in matrix form
  nu_time <- nu(Ttime)     # pointwise evaluation
  
  result <- truncAIPW_transMean_EF_cpp(
    nn = nn,
    time = Ttime,
    Q = Q,
    Fuz_mx = Fuz.mx,
    Gvz_mx = Gvz.mx,
    u = u,
    v = v,
    nu_time = nu_time,
    nuu = nuu,
    trim = trim
  )
  
  result <- lapply(result, as.vector)
  
  return(result)
}

  






### truncC AIPW ----------------------------------------------------------------------------

# Return applying truncC_AIPW to the function \nu(T) and the constant function 1
#### ! Need testing - already tested?
#### Double check plus or minus augmentation term - should be correct, but be aware if the results are not as expected.
truncC_AIPW_transMean <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                  X.name, Q.name, event.name, trim = 1e-7,
                                  u = NULL, v = NULL, d = NULL){
  
  ## compute the augmentation term
  nn <- nrow(dat)
  X <- as.numeric(dat[, X.name])
  Q <- as.numeric(dat[, Q.name])
  Delta <- as.integer(dat[, event.name])
  
  if(is.null(u)){
    u <- as.numeric(colnames(Fuz.mx))
  }
  if(is.null(v)){
    v <- as.numeric(colnames(Gvz.mx))
  }
  if(is.null(d)){
    d <- as.numeric(colnames(Sdz.mx))
  }
  
  # compute the IPCW weights
  X.res = X - Q
  Sdyz.vec <- diag(1 - CDF_eval_mx_cpp(X.res, 1 - Sdz.mx, d))
  id1 = (dat[, event.name] == 1)
  w_IPCW = rep(0, nrow(dat))
  w_IPCW[id1] = 1/pmax(Sdyz.vec[id1], trim)
  
  nuu <- nu(u)
  
  # apply truncAIPW
  efs = truncAIPW_transMean_EF(dat = dat, nu = nu, Fuz.mx = Fuz.mx, Gvz.mx = Gvz.mx, 
                               T.name = X.name, Q.name = Q.name, trim = trim)
  Num_AIPW = as.numeric(efs$Num_AIPW)
  Den_AIPW = as.numeric(efs$Den_AIPW)
  
  Num_IPW.Q = as.numeric(efs$Num_IPW.Q)
  Den_IPW.Q = as.numeric(efs$Den_IPW.Q)
  
  result_Aug_QD = aug_QD(nn, nuu, X, Q, Delta, Fuz.mx, Gvz.mx, Sdz.mx,
                         u, v, d, Sdyz.vec, trim)
  
  Aug_QD_nu = result_Aug_QD$Aug_QD_nu
  Aug_QD_const1 = result_Aug_QD$Aug_QD_const1
  
  truncC_nu = w_IPCW * Num_AIPW + Aug_QD_nu 
  truncC_const1 = w_IPCW * Den_AIPW + Aug_QD_const1
  
  
  return(list(truncC_nu = truncC_nu, 
              truncC_const1 = truncC_const1))
}






###########################
###  XGboost functions  ###
###########################

xgb_cv_wsq <- function(x,
                       y,
                       weights,
                       k_folds = NULL,
                       params,
                       ntrees_max = 500,
                       print_every_n = 100,
                       nthread = NULL,
                       verbose = FALSE){
  
  if(length(y) != nrow(x)){
    stop("The length of y and the row number of x are not the same!")
  }
  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10,length(y)/4)))
  }
  if (is.null(weights)) {
    weights = rep(1, length(y))
  }
  
  if (is.null(nthread)){
    nthread = parallel::detectCores()
  }
  
  n = length(y)
  folds <- cvFolds(n, k_folds)
  cv_errors <- matrix(NA, nrow = ntrees_max, ncol = k_folds)
  
  for(k in 1:k_folds){
    id.test = folds$subsets[folds$which == k]
    id.train = folds$subsets[folds$which != k]
    x.train = x[id.train, , drop = FALSE]
    x.test= x[id.test, , drop = FALSE]
    y.train = y[id.train]
    y.test= y[id.test]
    weights.test = weights[id.test]
    weights.train = weights[id.train]
    
    dtrain <- xgboost::xgb.DMatrix(data = x.train, label = y.train)
    dtest <- xgboost::xgb.DMatrix(data = x.test, label = y.test)
    
    
    # Define the custom weighted squared loss objective function
    wsq_loss_train <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      residuals <- preds - labels
      grad <- 2 * residuals * weights.train
      hess <- 2 * weights.train
      return(list(grad = grad, hess = hess))
    }
    
    
    # fit the model and compute the error on testing set
    xgb_fit <- xgboost::xgb.train(params = params, data = dtrain,
                                  nrounds = ntrees_max, 
                                  maximize = FALSE,
                                  obj = wsq_loss_train,
                                  verbose = verbose,
                                  print_every_n = print_every_n, 
                                  nthread = nthread)
    
    for (i in 1:ntrees_max) {
      preds <- predict(xgb_fit, newdata = dtest, iterationrange = c(1, i+1))
      weighted_errors <- weights.test * (preds - y.test)^2
      cv_errors[i, k] <- mean(weighted_errors)
    }
    
  }
  
  
  mean_cv_errors <- rowMeans(cv_errors)
  best_nrounds <- which.min(mean_cv_errors)
  best_cv_error <- mean_cv_errors[best_nrounds]
  
  return(list(best_nrounds = best_nrounds, 
              best_cv_error = best_cv_error, 
              cv_errors = cv_errors))
}






cvboost_wsq <- function(x,
                        y,
                        weights = NULL,
                        k_folds = NULL,
                        ntrees_max = 500,
                        num_search_rounds = 20,
                        print_every_n = 100,
                        nthread = NULL,
                        verbose = FALSE) {
  
  if (is.null(k_folds)) {
    k_folds = floor(max(3, min(10, length(y) / 4)))
  }
  if (is.null(weights)) {
    weights = rep(1, length(y))
  }
  
  best_param = list()
  best_seednumber = 1234
  best_loss = Inf
  best_nrounds = 0
  
  if (is.null(nthread)) {
    nthread = parallel::detectCores()
  }
  
  for (iter in 1:num_search_rounds) {
    # the parameters used for simulation and previous application
    param <- list(subsample = sample(c(0.5, 0.75, 1), 1),
                  colsample_bytree = sample(c(0.6, 0.8, 1), 1),
                  eta = sample(c(5e-3, 1e-2, 0.015, 0.025, 5e-2, 8e-2, 1e-1, 2e-1), 1),
                  max_depth = sample(c(3:20), 1),
                  gamma = runif(1, 0.0, 0.2),
                  min_child_weight = sample(1:20, 1),
                  max_delta_step = sample(1:10, 1))
    
    # # The tuning parameters to improve the smoothness of CATE 3D plots for HAAS application
    # param <- list(subsample = sample(c(0.5, 0.7, 0.9), 1),
    #               colsample_bytree = sample(c(0.6, 0.8, 1), 1),
    #               eta = sample(c(1e-3, 0.005, 0.01, 0.02, 5e-2, 1e-1), 1),
    #               max_depth = sample(c(2:6), 1),
    #               gamma = sample(c(0, 0.5, 1, 2, 3, 5), 1),
    #               min_child_weight = sample(seq(10, 50, by = 5), 1),
    #               max_delta_step = sample(seq(0, 10, by = 2), 1))
    
    seed_number = sample.int(100000, 1)[[1]]
    set.seed(seed_number)
    
    xgb_cvfit <- xgb_cv_wsq(x = x,
                            y = y,
                            weights = weights,
                            k_folds = k_folds,
                            params = param,
                            ntrees_max = ntrees_max,
                            print_every_n = print_every_n,
                            nthread = nthread,
                            verbose = verbose)
    
    min_loss = xgb_cvfit$best_cv_error
    
    if (min_loss < best_loss) {
      best_loss = min_loss
      best_seednumber = seed_number
      best_param = param
      best_xgb_cvfit = xgb_cvfit
      best_nrounds = xgb_cvfit$best_nrounds
    }
  }
  
  set.seed(best_seednumber)
  
  dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
  
  wsq_loss_final <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    residuals <- preds - labels
    grad <- 2 * residuals * weights
    hess <- 2 * weights
    return(list(grad = grad, hess = hess))
  }
  
  xgb_train_args = list(data = dtrain,
                        params = best_param,
                        nthread = nthread,
                        nrounds = best_nrounds,
                        obj = wsq_loss_final)
  
  xgb_fit <- do.call(xgboost::xgb.train, xgb_train_args)
  
  ret = list(xgb_fit = xgb_fit,
             best_xgb_cvfit = best_xgb_cvfit,
             best_seednumber = best_seednumber,
             best_param = best_param,
             best_loss = best_loss,
             best_nrounds = best_nrounds,
             best_xgb_train_args = xgb_train_args)
  
  class(ret) <- "cvboost_wsq"
  
  return(ret)
}

#' predict for cvboost_wsq
#'
#' @param object a cvboost_wsq object
#' @param newx covariate matrix to make predictions on. If null, return the predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @return vector of predictions
#' @export
predict.cvboost_wsq<- function(object,
                               newx = NULL,
                               ...) {
  
  if (is.null(newx)) {
    return(object$best_xgb_cvfit$pred)
  } else {
    dtest <- xgboost::xgb.DMatrix(data = newx)
    return(predict(object$xgb_fit, newdata = dtest))
  }
}




#########################
###  Extra functions  ###
#########################



CDF_eval <- function(time.eval, CDF.mx){
  if(length(time.eval) != nrow(CDF.mx)){
    stop("The number of time points does not equal the number of subjects!")
  }
  
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  
  probs = rep(NA, length(time.eval))
  for(i in 1:length(time.eval)){
    id = findInterval(time.eval[i], c(0, jumps, Inf))
    probs[i] = CDF.mx[i,id]
  }
  
  return(probs)
}



# Return a matrix of CDF(time.eval_j|Z_i) - (i,j)-th element
CDF_eval.mx <- function(time.eval, CDF.mx){
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  id = findInterval(time.eval, c(0, jumps, Inf))
  probs = CDF.mx[,id]
  
  if(length(time.eval)>1 & is.null(dim(probs))){
    names(probs) = time.eval
  }else if(length(time.eval)>1 & (!is.null(dim(probs)))){
    colnames(probs) = time.eval
  }
  
  return(probs)
}




# functions for computing L-S integrals ---------------------------------------------

int_fmx_dF <- function(v, f.mx, F.mx){
  
  if(mean(dim(f.mx) == dim(F.mx))<1){
    stop("The dimensions of f.mx and F.mx are not the same!")
  }
  
  jumps = as.numeric(colnames(F.mx))
  
  dF.mx = F.mx - cbind(0, F.mx[,-ncol(F.mx)]) 
  # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)]) 
  id = findInterval(v, c(0, jumps, Inf))
  vals = cbind(0, f.mx*dF.mx)
  
  inte = int_sum(id, vals)
  
  rownames(inte) = rownames(F.mx)
  colnames(inte) = v
  
  return(inte)
}


# # compute the integral for one subject at a single time point, i.e., f.mx and F.mx are just vectors, and v is a given time
# int_fmx_dF.single <- function(v, f.vec, F.vec, jumps){
#     if(length(f.vec) != length(F.vec) | length(f.vec) != length(jumps)){
#         stop("The lengths of f.mx, F.mx, and jumps are not all the same!")
#     }
#     
#     dF.vec = F.vec - c(0, F.vec[-length(F.vec)])
#     id = findInterval(v, c(0, jumps, Inf))
#     vals = c(0, f.vec*dF.vec)
#     
#     inte = sum(vals[1:id])
#     
#     return(inte)
# }



# int_fmx_dS <- function(v, f.mx, F.mx){
#     
#     if(mean(dim(f.mx) == dim(F.mx))<1){
#         stop("The dimensions of f.mx and F.mx are not the same!")
#     }
#     
#     jumps = as.numeric(colnames(F.mx))
#     dF.mx = F.mx - cbind(1, F.mx[,-ncol(F.mx)]) 
#     # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)]) 
#     id = findInterval(v, c(0, jumps, Inf))
#     vals = cbind(0, f.mx*dF.mx)
#     
#     inte = int_sum(id, vals)
#     
#     rownames(inte) = rownames(F.mx)
#     colnames(inte) = v
#     
#     return(inte)
# }






# function for returning a matric for \int_v^\infty fmx dF
int_infty_fmx_dF <- function(v, f.mx, F.mx){   
  if(mean(dim(f.mx) == dim(F.mx))<1){
    stop("The dimensions of f.mx and F.mx are not the same!")
  }
  jumps = as.numeric(colnames(F.mx))
  if(nrow(f.mx)==1){
    dF.mx = unlist(F.mx) - unlist(cbind(0, matrix(F.mx[,-ncol(F.mx)], nrow=1)))     #CHANGED TO ADD UNLIST HERE
    dF.mx = matrix(dF.mx, nrow=1)
  }else{
    dF.mx = unlist(F.mx) - unlist(cbind(0, F.mx[,-ncol(F.mx)]))
    dF.mx = matrix(dF.mx, nrow=1)
  }
  
  id = findInterval(v, c(0, jumps, Inf))
  vals = cbind(0, f.mx*dF.mx)
  
  inte = int_sum_infty(id, vals)
  
  # rownames(inte) = rownames(F.mx)
  # colnames(inte) = v
  
  return(inte)
}



# id: a single index denoting the 
int_sum.single <- function(id, vals){
  if(id==1){
    inte = vals[,1]
  }else{
    inte = rowSums(vals[,(1:id)])
  }
  return(inte)
}

int_sum <- function(id, vals){
  inte = sapply(id, int_sum.single, vals = vals)
  return(inte)
}


# id: a single index denoting the 
int_sum_infty.single <- function(id, vals){
  if(id == ncol(vals)){
    inte = rep(0, nrow(vals))
  }else if(id == ncol(vals)-1){
    inte = vals[,ncol(vals)]
  }else{
    if(nrow(vals) == 1){
      inte = sum(vals[,((id+1):ncol(vals))])
    }else{
      inte = rowSums(vals[,((id+1):ncol(vals))])
    }
  }
  return(inte)
}

int_sum_infty <- function(id, vals){
  inte = sapply(id, int_sum_infty.single, vals = vals)
  return(inte)
}



CDF_eval.mx <- function(time.eval, CDF.mx){
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  id = findInterval(time.eval, c(0, jumps, Inf))
  probs = CDF.mx[,id]
  
  if(length(time.eval)>1 & is.null(dim(probs))){
    names(probs) = time.eval
  }else if(length(time.eval)>1 & (!is.null(dim(probs)))){
    colnames(probs) = time.eval
  }
  
  return(probs)
}


ns_mx <- function(data, cov.names, cov.names.binary, df){
  
  X_ns = do.call(cbind, lapply(cov.names, function(cov){splines::ns(data[, cov], df = df)}))
  X_ns_wbinary = cbind(X_ns, data[, cov.names.binary])  # combine with the binary covariates
  dim_ns = dim(X_ns)[2]
  X_ns_linear_interaction = stats::model.matrix(~.*.-1, data.frame(X_ns_wbinary)) # pairwise interaction (not including squared term for each column), including the first order terms
  X_ns_sq = do.call(cbind, lapply(1:dim_ns, function(col){matrix(X_ns[,col]^2)})) # squared term for each column
  X_ns = cbind(X_ns_linear_interaction, X_ns_sq)
  
  return(X_ns)
}


int_fmx_dF <- function(v, f.mx, F.mx){
  
  if(mean(dim(f.mx) == dim(F.mx))<1){
    stop("The dimensions of f.mx and F.mx are not the same!")
  }
  
  jumps = as.numeric(colnames(F.mx))
  dF.mx = F.mx - cbind(0, F.mx[,-ncol(F.mx)]) 
  # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)]) 
  id = findInterval(v, c(0, jumps, Inf))
  
  vals = cbind(0, f.mx*dF.mx)
  
  inte = int_sum(id, vals)
  
  rownames(inte) = rownames(F.mx)
  colnames(inte) = v
  
  return(inte)
}



CDF_eval <- function(time.eval, CDF.mx){
  if(length(time.eval) != nrow(CDF.mx)){
    stop("The number of time points does not equal the number of subjects!")
  }
  
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  
  probs = rep(NA, length(time.eval))
  id = rep(NA, length(time.eval))
  for(i in 1:length(time.eval)){
    id = findInterval(time.eval[i], c(0, jumps, Inf))
    if (id == 1){
      id=2
    }
    probs[i] = CDF.mx[i,id]
  }
  
  return(probs)
}



# Return a matrix of CDF(time.eval_j|Z_i) - (i,j)-th element
CDF_eval.mx <- function(time.eval, CDF.mx){
  jumps = as.numeric(colnames(CDF.mx))
  CDF.mx = cbind(0, CDF.mx)
  id = findInterval(time.eval, c(0, jumps, Inf))
  probs = CDF.mx[,id]
  
  if(length(time.eval)>1 & is.null(dim(probs))){
    names(probs) = time.eval
  }else if(length(time.eval)>1 & (!is.null(dim(probs)))){
    colnames(probs) = time.eval
  }
  
  return(probs)
}

# x: a vector
# return max(x_i, trim) if x_i>0, and min(x_i, -trim) if x_i<0
bound_away_zero <- function(x, trim){
  sign_x = sign(x)
  sign_x[x == 0] <- 1
  y = sign_x * pmax(abs(x), trim)
  
  return(y)
}
