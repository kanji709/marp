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

  rmse = list(esr = est, ha.aic=ha.aic, ha.bstrp=ha.bstrp, ha.best=ha.best, rmse.aic=rmse.aic, rmse.best=rmse.best,rmse.bstrp=rmse.bstrp)

dput(rmse,"simulations/result.txt")

  

