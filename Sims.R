library("FAdist")
library("statmod")

source('weight.R')
source('hazard.R')
source('bootstrapp.R')
source('estima.R')
source('data.R')
source('rmse.R')


##--------------Simulate Data Set------------------------------------------------------------

  #--------------Sample size------------------------------------#
  n = 50
  
  #--------------Number of simulations at each point------------------------------------#
  p = 10
  
  #--------------Number of bootstrapped resamples------------------------------------#
  B = 10
  
  #--------------Length of the for loop in MLE------------------------------------#
  m = 1000
  
  #--------------Fixed mu, change value of sig------------------------------------#
  mu = 3
  
  sig = 1 
  
  #--------------Interevent times for simulation------------------------------------#
  t <- seq(1, 100)
  
  
  #--------------Date for simulations------------------------------------#
  pts <- gen.data(n,mu,sig,p)
  
  #--------------Results from simulations------------------------------------#
  est = estima(pts[,1],B,m)
  for (i in 2:p){
    
    est <- mapply(c,est, estima(pts[,i],B,m),SIMPLIFY = F, USE.NAMES = T)
  }
# lapply(1:3, function(i) seq(df[1,i], df[2,i], df[3,i]))
  #--------------Mdoel Averaged hazard rates with AIC weights------------------------------------#
  ha.aic <- matrix(est$aic,nrow = length(t),ncol = p*length(sig))
  
  #--------------Model Averaged hazard rates with bootstrapped propotion weights------------------------------------#
  ha.bstrp <- matrix(est$bstrp,nrow = length(t),ncol = p*length(sig))
  
  #--------------Model Averaged hazard rates from median of bootstrapped best estimates------------------------------------#
  ha.best <- matrix(est$best,nrow = length(t),ncol = p*length(sig))
  
  #-------------- Calculating RMSE------------------------------------#
  rmse.aic <- rmse(ha.aic,mu,sig,t)
  
  rmse.best <- rmse(ha.best,mu,sig,t)
  
  rmse.bstrp <- rmse(ha.bstrp,mu,sig,t)

  rmse = list(out = est, aic=ha.aic, bstrp = ha.bstrp, best = ha.best, rmse_aic = rmse.aic, rmse_best=rmse.best,rmse_bstrp=rmse.bstrp)

dput(rmse,"simulations/result.txt")

  

