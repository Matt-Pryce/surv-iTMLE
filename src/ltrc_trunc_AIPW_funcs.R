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