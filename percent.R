percentile <-  function(data,B,m,t){

#------------------------------------------ Model-Averaged Hazard Rates with AIC Weights -----------------------------------------------------------------------#  
  #Estimate parameters and obtain AIC weights
  prelim <-  weights(data,m)
  
  par.hat <- c(prelim$Par_1,prelim$Par_2)
  
  aic.weights <- prelim$Weights_AIC
  
  haz.hat <- hazard(par.hat,t)
  
  haz.aic <- haz.hat%*%aic.weights
  
#------------------------------------------ Model-Averaged Hazard Rates with Bootstarpped Proportion Weights -----------------------------------------------------------------------#  
  #Parametric or Non-parametric bootstrap
  resamples <- bstrp.np(data,B)
  
  #All the estimated values for each model of each bootstrapped resamples
  temp <- apply(resamples,2,function(x) weights(x,m))
  
  #Store the index for a single best model of each bootsrtapped resamples  
  id <- sapply(temp,'[[',8,USE.NAMES = F,simplify = T)
  
  #Proportion weights
  prop.weights <- sapply(1:6,function(i) length(which(i==id))/B)
  
  haz.prop <- haz.hat%*%prop.weights


#------------------------------------------------------- Percentile Model-Averaged Confidence Intervals --------------------------------------------------------------------------#
  
  #Extract the paramters from the results above  
  par <- rbind(sapply(temp,'[[',1),sapply(temp,'[[',2))
  
  #Calculate the hazard rates   
  haz <- apply(par,2,function(x) as.data.frame(hazard(x,t))) 
  
  #Hazard rates from the data generated model   
  logn <- sapply(1:B,function(i) haz[[i]][,5])
  
  #Median of the Hazard rates from the data generated model
  haz.gen <- apply(logn,1,function(x) quantile(x,probs = 0.5))
  
  #Hazard rates from the single best model of each bootstrapped resamples  
  best <- sapply(1:B,function(i) haz[[i]][,id[i]])
  
  #Percentile CI
  percent.ci <- apply(best,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))

#------------------------------------------ Model-Averaged Hazard Rates from median of bootstrapped estimates -----------------------------------------------------------------------#  
  haz.med <- percent.ci[2,]
  
#------------------------------------------ output -----------------------------------------------------------------------#  
  
  out <- list('mahr.aic' = t(haz.aic),'mahr.prop' = t(haz.prop), 'mahr.med' = haz.med, 'hr.gen' = haz.gen, 'aic.weights' = aic.weights, 'id' = id, 'prop.weights' = prop.weights, 'lower' = percent.ci[1,],'upper' = percent.ci[3,])
  
  return(out)
  
}
  

