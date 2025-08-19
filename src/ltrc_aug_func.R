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