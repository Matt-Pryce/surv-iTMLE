

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
    
  }else{
    stop("This T model is not implemented in this function!")
  }
  
  Fuz.mx = pmin(Fuz.mx, 1-trim)
  
  return(Fuz.mx)
  # return(list(fit.pred = fit.pred,
  #             Fuz.mx0 = Fuz.mx0,
  #             Fuz.mx = Fuz.mx))
  # return(cv.fit)
}


