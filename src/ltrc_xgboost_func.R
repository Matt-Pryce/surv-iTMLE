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
    # # the parameters used for simulation and previous application
    # param <- list(subsample = sample(c(0.5, 0.75, 1), 1),
    #               colsample_bytree = sample(c(0.6, 0.8, 1), 1),
    #               eta = sample(c(5e-3, 1e-2, 0.015, 0.025, 5e-2, 8e-2, 1e-1, 2e-1), 1),
    #               max_depth = sample(c(3:20), 1),
    #               gamma = runif(1, 0.0, 0.2),
    #               min_child_weight = sample(1:20, 1),
    #               max_delta_step = sample(1:10, 1))
    
    # The tuning parameters to improve the smoothness of CATE 3D plots for HAAS application
    param <- list(subsample = sample(c(0.5, 0.7, 0.9), 1),
                  colsample_bytree = sample(c(0.6, 0.8, 1), 1),
                  eta = sample(c(1e-3, 0.005, 0.01, 0.02, 5e-2, 1e-1), 1),
                  max_depth = sample(c(2:6), 1),
                  gamma = sample(c(0, 0.5, 1, 2, 3, 5), 1),
                  min_child_weight = sample(seq(10, 50, by = 5), 1),
                  max_delta_step = sample(seq(0, 10, by = 2), 1))
    
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