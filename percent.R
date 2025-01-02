percentile <-  function(data,B,t,m){

#------------------------------------------ Model-Averaged Hazard Rates with Bootstarpped Proportion Weights -----------------------------------------------------------------------#  
  #Parametric or Non-parametric bootstrap
  resamples <- replicate(B,sample(data,size = length(data),replace = TRUE))
  
  #All the estimated values for each model of each bootstrapped resamples
  temp <- apply(resamples,2,function(x) prelim(x,t,m))
  
  #Store the index for a single best model of each bootsrtapped resamples  
  id <- sapply(temp,'[[',9,USE.NAMES = F,simplify = T)
  
  #Proportion weights
  prop.weights <- sapply(1:6,function(i) length(which(i==id))/B)
  

#------------------------------------------------------- Percentile Model-Averaged Confidence Intervals --------------------------------------------------------------------------#
  
  #Extract the paramters from the results above  
  par <- rbind(sapply(temp,'[[',1),sapply(temp,'[[',2))
  
  #Calculate the hazard rates   
  haz <- sapply(temp,'[[',6) 
  
  #Hazard rates from the data generated model   
  logn <- sapply(1:B,function(i) matrix(haz[,i],length(t),6)[,5])
  
  #Median of the Hazard rates from the data generated model
  haz.gen <- apply(logn,1,function(x) quantile(x,probs = 0.5))
  
  #Hazard rates from the single best model of each bootstrapped resamples  
  best <- sapply(1:B,function(i) matrix(haz[,i],length(t),6)[,id[i]])
  
  #Percentile CI
  percent.ci <- apply(best,1,function(x) quantile(x,probs = c(0.025,0.5,0.975)))

#------------------------------------------ Model-Averaged Hazard Rates from median of bootstrapped estimates -----------------------------------------------------------------------#  
  haz.med <- percent.ci[2,]
  
#------------------------------------------ output -----------------------------------------------------------------------#  
  
  out <- list('id' = id, 'prop.weights' = prop.weights, 'haz.med' = haz.med, 'haz.gen' = haz.gen, 'lower' = percent.ci[1,],'upper' = percent.ci[3,])
  
  return(out)
  
}
  

