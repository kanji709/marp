estima <- function(sims, B, m){
  
  
  ##-------------Estimated parameters and hazard rates for each model-------------------------
  
  
  #Estimated Paras, logLik, AIC, BIC and Weights with MLE
  est <- weights(sims, m)
  
  
  
  #values of estimated parameters
  par.hat <- matrix(c(est$Par_1, est$Par_2), nrow = 6, ncol = 2)
  
  
  
  #Hazard rate for each model with estimated paramters
  hazard.hat <- hazard(par.hat, t)
  
  ##-----------------Hazard rates from log normal model------------------------------
  ha.lnorm <- hazard.hat[5,]
  
  
  ##-----------------Hazard rates from log normal model------------------------------
  id <-  as.numeric(which.min(est$AIC))
  
  ha.sel <- hazard.hat[id,]
  
  ##-----------------Model averaged hazard rates with AIC weights------------------------------
  
  
  #AIC Weights for each canidate model
  weight.aic <- matrix(c(est$Weights_AIC), nrow = 1, ncol = 6)
  
  
  
  #Model Averaged Hazard Rate with AIC weights
  hazard.aic <- weight.aic %*% hazard.hat
  
  
  
  ##---------------Percentile bootstrap-------------------------------
  #Results from boostrap
  
  #start.time <- Sys.time()
  result <- bootstrap(sims, t, B, m)
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #time.taken
  
  
  #Call index for the best model in each bootstrapped resample
  id.best <- result$id
  
  
  #Compute bootstrap propotion weights
  weights.bstrp <-
    matrix(c(
      length(which(id.best == 1)),
      length(which(id.best == 2)),
      length(which(id.best == 3)),
      length(which(id.best == 4)),
      length(which(id.best == 5)),
      length(which(id.best == 6))
    ) / B,ncol = 6)
  
  
  
  #Model averaged hazard rate with bootstrapped propotion weights 
  hazard.bstrp <- weights.bstrp %*% hazard.hat
  
  
  
  
  ##------------Model averaged hazard rates as the median of best estimat from bootstrap---------------
  
  
  #Get the median of 50% quantile of bootstrapped hazard rates
  median <-apply(result$ma, 2, function(x)
      quantile(x, 0.5))
  hazard.best <- matrix(median, nrow = 1)
  
  
  
 out <- list(lnorm = ha.lnorm,
             sel = ha.sel,
             w_aic = weight.aic,
             aic = hazard.aic,
             w_bstrp = weights.bstrp,
             bstrp = hazard.bstrp,
             best = hazard.best) 
 return(out)
}

