est_G <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name, event.name, cov.names, trim = 0, weights = NULL,
                  formula.survPen = NULL, tau = NULL, trunc = TRUE,
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, df = 5){
  
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
  
  #Consider from here
  if(is.null(tau)){
    tau = max(c(dat.fit[,time.name], dat.est[,time.name])) + 1
  }
  
  #Adding in so we can look at earlier times
  dat.fit[ ,time.name] <- pmin(dat.fit[ ,time.name],(tau-0.01))
  dat.fit[ ,Q.name] <- pmin(dat.fit[ ,Q.name],(tau-0.01))
  dat.est[ ,time.name] <- pmin(dat.est[ ,time.name],(tau-0.01))
  dat.est[ ,Q.name] <- pmin(dat.est[ ,Q.name],(tau-0.01))
  
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