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