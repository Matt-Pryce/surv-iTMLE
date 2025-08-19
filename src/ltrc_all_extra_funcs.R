
###############################
###  Nuisance function est  ###
###############################

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




############################
###  AIPW and aug funcs  ###
############################


truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx,
                                   T.name, Q.name, trim = 0){
  
  if(nrow(dat) != nrow(Fuz.mx) | nrow(dat) != nrow(Gvz.mx)){
    stop("The number of rows of dat, Fuz.mx and Gvz.mx are not the same. ")
  }
  
  time = dat[, T.name]
  Q = dat[, Q.name]
  
  u = as.numeric(colnames(Fuz.mx))  # jumps.T
  v = as.numeric(colnames(Gvz.mx))  # jumps.Q
  
  tau2 = max(v)+1e-10
  
  Gtz = CDF_eval(time, Gvz.mx)
  Gqz = CDF_eval(Q, Gvz.mx)
  Fqz = CDF_eval(Q, Fuz.mx)
  
  DDen1 = 1/pmax(Gtz, trim)
  DDen2 = Fqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
  
  Fvz.mx = CDF_eval.mx(v, Fuz.mx)
  
  nn = nrow(dat)
  atrisk.mx = matrix(nrow = nn, ncol = length(v))
  for(i in 1:nn){
    atrisk.mx[i,] = (Q[i] <= v & v < time[i])    # at risk indicator for subject i at times in jumps.Q
  }
  f.mx = atrisk.mx*Fvz.mx/(pmax(Gvz.mx, trim)^2*pmax(as.matrix(1-Fvz.mx), trim))
  DDen3 = as.vector(int_fmx_dF(tau2, f.mx, Gvz.mx))
  
  NNum1 = nu(time)/pmax(Gtz, trim)
  
  nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
  mqz = diag(int_fmx_dF(Q, nuu.mx, Fuz.mx))
  NNum2 = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
  
  mvz.mx = int_fmx_dF(v, nuu.mx, Fuz.mx)
  fnu.mx = atrisk.mx*mvz.mx/(pmax(Gvz.mx, trim)^2*pmax(as.matrix(1-Fvz.mx), trim))
  NNum3 = as.vector(int_fmx_dF(tau2, fnu.mx, Gvz.mx))
  
  Num_AIPW = NNum1 + NNum2 - NNum3
  Num_AIPW = pmax(Num_AIPW,0)       #Added for positive weights
  Den_AIPW = DDen1 + DDen2 - DDen3
  Den_AIPW = pmax(Den_AIPW,0)       #Added for positive weights
  
  # For the estimating equation of the IPW and Reg.T1, Reg.T2 estimators
  DDen4 = 1/pmax(1-Fqz, trim)
  NNum4 = nu(time) + mqz/pmax(1-Fqz, trim)
  
  tau.Tmax = max(u) + 1
  Enutz = as.vector(int_fmx_dF(tau.Tmax, nuu.mx, Fuz.mx))
  NNum5 = Enutz/pmax(1-Fqz, trim)
  
  
  
  return(list(Num_AIPW = Num_AIPW, Den_AIPW = Den_AIPW,
              Num_IPW.Q = NNum1, Den_IPW.Q = DDen1,
              Num_Reg.T1 = NNum4, Num_Reg.T2 = NNum5, Den_Reg = DDen4))
}




truncC_AIPW_transMean <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                  X.name, Q.name, event.name, trim = 1e-7){
  
  X.res = dat[,X.name] - dat[,Q.name]
  
  # apply truncAIPW
  efs = truncAIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, X.name, Q.name, trim)
  Num_AIPW = efs$Num_AIPW
  Den_AIPW = efs$Den_AIPW
  
  Num_IPW.Q = efs$Num_IPW.Q
  Den_IPW.Q = efs$Den_IPW.Q
  
  # compute the IPCW weights
  Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))
  id1 = (dat[, event.name] == 1)
  w_IPCW = rep(0, nrow(dat))
  w_IPCW[id1] = 1/pmax(Sdyz.vec[id1], trim)
  
  result_Aug_QD = aug_QD(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, Sdyz.vec,
                         X.name, Q.name, event.name, trim)
  Aug_QD_nu = result_Aug_QD$Aug_QD_nu
  Aug_QD_const1 = result_Aug_QD$Aug_QD_const1
  
  truncC_nu = w_IPCW * Num_AIPW + pmax(Aug_QD_nu,0)
  truncC_const1 = w_IPCW * Den_AIPW + pmax(Aug_QD_const1,0)
  
  
  return(list(truncC_nu = truncC_nu,
              truncC_const1 = truncC_const1))
}


aug_QD <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, Sdyz.vec = NA,
                   X.name, Q.name, event.name, trim = 1e-7){
  if(sum(is.na(Sdyz.vec))>0){
    X.res = dat[,X.name] - dat[,Q.name]
    Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))
  }
  
  aug_result_1 = aug_QD_1(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, X.name, Q.name, event.name, trim)
  aug_result_2 = aug_QD_2(dat, nu, Fuz.mx, Gvz.mx, Sdyz.vec, Q.name, event.name, trim)
  aug_result_3 = aug_QD_3(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, X.name, Q.name, event.name, trim)
  
  Aug_QD_nu_1 = aug_result_1$aug_1
  Aug_QD_nu_2 = aug_result_2$aug_2
  Aug_QD_nu_3 = aug_result_3$aug_3
  Aug_QD_nu = Aug_QD_nu_1 + Aug_QD_nu_2 - Aug_QD_nu_3
  
  Aug_QD_const1_1 = aug_result_1$aug_1_const1
  Aug_QD_const1_2 = aug_result_2$aug_2_const1
  Aug_QD_const1_3 = aug_result_3$aug_3_const1
  Aug_QD_const1 = Aug_QD_const1_1 + Aug_QD_const1_2 - Aug_QD_const1_3
  
  # Aug_11 = aug_result_1$aug_11
  # Aug_21 = aug_result_1$aug_21
  # Aug_31 = aug_result_1$aug_31
  
  return(list(Aug_QD_nu = Aug_QD_nu,
              Aug_QD_const1 = Aug_QD_const1))
}



aug_QD_1 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                     X.name, Q.name, event.name, trim = 1e-7){
  
  nn = nrow(dat)
  X = dat[,X.name]
  Q = dat[,Q.name]
  X.res = dat[,X.name] - dat[,Q.name]
  Delta = dat[, event.name]
  
  u = as.numeric(colnames(Fuz.mx))
  v = as.numeric(colnames(Gvz.mx))
  d = as.numeric(colnames(Sdz.mx))
  
  tau.Tmax = max(u) + 1e-5
  tau.Dmax = max(d) + 1e-5
  
  # compute the IPCW weights
  Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))   # S_D(X-Q|Q,A,Z)
  id0 = (dat[, event.name] == 0)
  w_D = rep(0, nrow(dat))
  w_D[id0] = 1/pmax(Sdyz.vec[id0], trim)
  
  ## compute aug_1
  # compute aug_11
  Fxz.vec = diag(as.matrix(CDF_eval.mx(X, Fuz.mx)))  # F(X|A,Z)
  Guz.mx = CDF_eval.mx(u, Gvz.mx)    # G(u|A,Z).mx
  
  
  nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
  
  ind_Xu.mx = matrix(nrow = nn, ncol = length(u))     #  I(X <= u)
  for(i in 1:nn){
    ind_Xu.mx[i,] = as.numeric(X[i] <= u)
  }
  fu1.mx = ind_Xu.mx * nuu.mx / pmax(Guz.mx, trim)
  aug_11 = w_D * as.vector(int_fmx_dF(tau.Tmax, fu1.mx, Fuz.mx)) / pmax(1-Fxz.vec, trim)   # this is faster than diag(int_infty_fmx_dF(X,...))
  
  fu1_const1.mx = ind_Xu.mx / pmax(Guz.mx, trim)
  aug_11_const1 = w_D * as.vector(int_fmx_dF(tau.Tmax, fu1_const1.mx, Fuz.mx)) / pmax(1-Fxz.vec, trim)   # this is faster than diag(int_infty_fmx_dF(X,...))
  
  
  # compute aug_12
  # start_time <- Sys.time()
  # system.time({
  wnu.mx = nuu.mx / pmax(Guz.mx, trim)
  int_wnu_dF = matrix(nrow = nn, ncol = length(d))  # int \ind(Q+d<=u)\nu(u)/G(u|A,Z) dF(u|A,Z)
  F_Qu_z.mx = matrix(nrow = nn, ncol = length(d)) # F(Q+u|A,Z)
  
  wnu_const1.mx = 1 / pmax(Guz.mx, trim)
  int_wnu_const1_dF = matrix(nrow = nn, ncol = length(d))  # int \ind(Q+d<=u) 1/G(u|A,Z) dF(u|A,Z)
  
  for(i in 1:nn){
    wnu.mx.i = matrix(wnu.mx[i,], nrow=1)
    wnu_const1.mx.i = matrix(wnu_const1.mx[i,], nrow=1)
    Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
    colnames(Fuz.mx.i) = colnames(Fuz.mx)
    
    int_wnu_dF[i,] = int_infty_fmx_dF(Q[i]+d, wnu.mx.i, Fuz.mx.i)
    int_wnu_const1_dF[i,] = int_infty_fmx_dF(Q[i]+d, wnu_const1.mx.i, Fuz.mx.i)
    
    F_Qu_z.mx[i,] = CDF_eval.mx(Q[i]+d, matrix(unlist(Fuz.mx.i), nrow=1))   #cHANGED TO ADD UNLIST
  }
  # })
  # # now_time <- Sys.time()
  # # now_time - start_time      # ~ 10.66324 secs
  # 
  FDdz.mx = 1-Sdz.mx
  ind_dXQ.mx = matrix(nrow = nn, ncol = length(d))
  for(i in 1:nn){
    ind_dXQ.mx[i,] = as.numeric(d <= X[i]-Q[i])
  }
  
  fd2.mx = ind_dXQ.mx * int_wnu_dF/ (pmax(1-F_Qu_z.mx, trim) * (pmax(Sdz.mx, trim))^2)
  aug_12 = as.vector(int_fmx_dF(tau.Dmax, fd2.mx, FDdz.mx))
  
  fd2_const1.mx = ind_dXQ.mx * int_wnu_const1_dF/ (pmax(1-F_Qu_z.mx, trim) * (pmax(Sdz.mx, trim))^2)
  aug_12_const1 = as.vector(int_fmx_dF(tau.Dmax, fd2_const1.mx, FDdz.mx))
  
  aug_1 = aug_11 - aug_12
  
  aug_1_const1 = aug_11_const1 - aug_12_const1
  
  
  return(list(aug_1 = aug_1, # aug_11 = aug_11, aug_12 = aug_12,
              aug_1_const1 = aug_1_const1))
}




aug_QD_2 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdyz.vec,
                     Q.name, event.name, trim = 1e-7){
  
  nn = nrow(dat)
  Q = dat[,Q.name]
  Delta = dat[, event.name]
  
  u = as.numeric(colnames(Fuz.mx))
  nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
  
  mqz = diag(int_fmx_dF(Q, nuu.mx, Fuz.mx))
  Gqz = CDF_eval(Q, Gvz.mx)
  Fqz = CDF_eval(Q, Fuz.mx)
  
  # aug_21 = (1-Delta)/pmax(Sdyz.vec, trim) * mqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))  # for testing purposes
  aug_2 = (1 - Delta/pmax(Sdyz.vec, trim)) * mqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))  
  
  aug_2_const1 = (1 - Delta/pmax(Sdyz.vec, trim)) * Fqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))
  
  
  return(list(aug_2 = aug_2, 
              aug_2_const1 = aug_2_const1))
}




aug_QD_3 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                     X.name, Q.name, event.name, trim = 1e-7){
  
  nn = nrow(dat)
  X = dat[,X.name]
  Q = dat[,Q.name]
  X.res = dat[,X.name] - dat[,Q.name]
  Delta = dat[, event.name]
  
  u = as.numeric(colnames(Fuz.mx))
  v = as.numeric(colnames(Gvz.mx))
  d = as.numeric(colnames(Sdz.mx))
  
  nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
  
  tau.Tmax = max(u) + 1e-5
  tau.Qmax = max(v) + 1e-5
  tau.Dmax = max(d) + 1e-5
  
  # compute the IPCW weights
  Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))   # S_D(X-Q|Q,A,Z)
  id0 = (dat[, event.name] == 0)
  w_D = rep(0, nrow(dat))
  w_D[id0] = 1/pmax(Sdyz.vec[id0], trim)
  
  # compute aug_31
  mvz.mx = int_fmx_dF(v, nuu.mx, Fuz.mx)
  Fvz.mx = CDF_eval.mx(v, Fuz.mx)
  
  F_Xv_z.mx = matrix(nrow = nn, ncol = length(v))  # used in the computation of aug_31
  ind_Qv.mx = matrix(nrow = nn, ncol = length(v))  # used in both aug_31 and aug_32
  for(i in 1:nn){
    Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
    colnames(Fuz.mx.i) = colnames(Fuz.mx)
    
    F_Xv_z.mx[i,] = CDF_eval.mx(pmin(v, X[i]), matrix(unlist(Fuz.mx.i), nrow=1))    #CHANGED TO UNLIST
    
    ind_Qv.mx[i,] = (Q[i] <= v)
  }
  
  fv31.mx = mvz.mx * ind_Qv.mx / ( pmax(1-F_Xv_z.mx, trim) * (pmax(Gvz.mx, trim))^2 )
  aug_31 = w_D * int_fmx_dF(tau.Qmax, fv31.mx, Gvz.mx)
  
  fv31_const1.mx = Fvz.mx * ind_Qv.mx / ( pmax(1-F_Xv_z.mx, trim) * (pmax(Gvz.mx, trim))^2 )
  aug_31_const1 = w_D * int_fmx_dF(tau.Qmax, fv31_const1.mx, Gvz.mx)
  
  
  ## Compute aug_32
  # start_time <- Sys.time()
  # system.time({
  aug32_int_dG = matrix(nrow = nn, ncol = length(d))   # int ... dG(v|A,Z) in aug_32
  aug32_int_dG_const1 = matrix(nrow = nn, ncol = length(d))   # int F(v|A,Z)/... dG(v|A,Z) in aug_32_const1
  v.vec = rep(v, length(d))
  d.vec = rep(d, each = length(v))
  for(i in 1:nn){
    Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
    colnames(Fuz.mx.i) = colnames(Fuz.mx)
    
    F_Qdv_z.mx = matrix(CDF_eval.mx(pmin(Q[i]+d.vec, v.vec), Fuz.mx.i), ncol = length(v), byrow = TRUE)
    F_Qdv_z.mx <- matrix(unlist(F_Qdv_z.mx),nrow=dim(F_Qdv_z.mx)[1],ncol = dim(F_Qdv_z.mx)[2])             #CHANGED HERE TO UNLIST
    
    fv32_i.mx = t(replicate(length(d), mvz.mx[i,]*ind_Qv.mx[i,])) / ( pmax(1-F_Qdv_z.mx, trim) * t( replicate(length(d), (pmax(Gvz.mx[i,], trim))^2) ) )
    aug32_int_dG[i,] = int_fmx_dF(tau.Qmax, fv32_i.mx, t(replicate(length(d), Gvz.mx[i,])) )
    
    fv32_i_const1.mx = t(replicate(length(d), unlist(Fvz.mx[i,]*ind_Qv.mx[i,]))) / ( pmax(1-F_Qdv_z.mx, trim) * t( replicate(length(d), (pmax(Gvz.mx[i,], trim))^2) ) )
    aug32_int_dG_const1[i,] = int_fmx_dF(tau.Qmax, fv32_i_const1.mx, t(replicate(length(d), Gvz.mx[i,])) )    #ADDED UNLIST
  }
  # })
  # now_time <- Sys.time()
  # now_time - start_time   # ~ 3min
  
  atrisk_XQ.mx = matrix(nrow = nn, ncol = length(d))
  for(i in 1:nn){
    atrisk_XQ.mx[i,] = (X[i]-Q[i] >= d)
  }
  
  FDdz.mx = 1-Sdz.mx
  fd32.mx = aug32_int_dG * atrisk_XQ.mx / ((pmax(Sdz.mx, trim))^2)
  aug_32 = int_fmx_dF(tau.Dmax, fd32.mx, FDdz.mx)
  
  fd32_const1.mx = aug32_int_dG_const1 * atrisk_XQ.mx / ((pmax(Sdz.mx, trim))^2)
  aug_32_const1 = int_fmx_dF(tau.Dmax, fd32_const1.mx, FDdz.mx)
  
  
  aug_3 = aug_31 - aug_32
  
  aug_3_const1 = aug_31_const1 - aug_32_const1
  
  return(list(aug_3 = aug_3, # aug_31 = aug_31, aug_32 = aug_32,
              aug_3_const1 = aug_3_const1))
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