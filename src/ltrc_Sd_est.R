est_Sd <- function(dat.fit, dat.est = dat.fit, time.eval, model, 
                   time.name, Q.name, event.name, cov.names = NULL, trim = 0, 
                   mtry = NA, ntree = NA, formula.survPen = NA, 
                   x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                   nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, 
                   df = 5, OOF = FALSE, nfolds.OOF = 10){
  
  # u = time.eval
  
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
      
      # x.fit <- as.matrix(dat.fit[,c(cov.names, cov.names.binary)])
      # x.est <- as.matrix(dat.est[,c(cov.names, cov.names.binary)])
    }

   yss = Surv(dat.fit$Y, dat.fit$delta.D)
   cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, alpha = alpha, lambda = lambda)
   fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est)
   Sdz.mx = t(fit.pred$surv)
   colnames(Sdz.mx) <- fit.pred$time

   Fdz.mx0 = 1 - Sdz.mx
   Sdz.mx = 1 - CDF_eval.mx(time.eval, Fdz.mx0)

 }
  
  Sdz.mx = pmax(Sdz.mx, trim)
  return(Sdz.mx)
  # return(dat.fit)
}