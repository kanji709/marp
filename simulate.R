simulate <- function(n,mu,sig,p,m,t,B){

#------------------------------------------ Simulate p datasets(size n) from log-normal(mu,sig) -------------------------------------------------------------------#  
   pts <- replicate(p,rlnorm(n,mu,sig))
  
#------------------------------------------ Preliminary analysis --------------------------------------------------------------------------------------------------#       
   intro <- sapply(1:p,function(i) prelim(pts[,i],t,m),simplify = F)
   
   par.hat <- rbind(sapply(intro,'[[',1),sapply(intro,'[[',2))

#------------------------------------------ Model-averaged hazard rates using AIC weights -------------------------------------------------------------------------#  
   ### AIC weights 
   aic.weights <- sapply(intro,'[[',7)
   
   ### Estimated hazard rates
   haz.hat <- sapply(intro,'[[',6)
   
   ### Model-averaged hazard rates using AIC weights
   haz.aic <- sapply(1:p,function(i) matrix(haz.hat[,i],length(t),6)%*%aic.weights[,i])
   
#------------------------------------------ Percentile bootstrapped confidence intervals and other types of model-averhaed estimates-------------------------------#   
   percent <- apply(pts,2,function(x) percentile(x,B,t,m))
   
   ### Bootstrapped propotion weights
   prop.weights <- sapply(percent,'[[',2)
   
   ### Model-averaged hazard rates using bootstrapped propotion weights
   haz.prop <- sapply(1:p,function(i) matrix(haz.hat[,i],length(t),6)%*%prop.weights[,i])
   
   ### Median of the Hazard rates from the data generated model
   haz.gen <- sapply(percent,'[[',4)
   
   ### Model-Averaged Hazard Rates using the median of best models from each bootstrapped resample
   haz.med <- sapply(percent,'[[',3)
   
   ### Bootstrapped percentile confidence intervals
   percent.ci <- sapply(1:p,function(i) cbind(percent[[i]]$lower,percent[[i]]$upper),simplify = F)
#------------------------------------------ Calculate the variance and bias for each estimates -----------------------------------------------------------------------#   
   ### Model-averaged hazard rates with AIC weights
   rmse.aic <- rmse(haz.aic,mu,sig,t)
   
   ### Model-averaged hazard rates with propotion weights
   rmse.prop <- rmse(haz.prop,mu,sig,t)
   
   ### Model-averaged hazard rates as median of bootstrapped estimates 
   rmse.med <- rmse(haz.med,mu,sig,t)
   
   ### Estimatd hazard rates from data generation model
   rmse.gen <- rmse(haz.gen,mu,sig,t)

   return(list('rmse.aic' = rmse.aic, 'rmse.prop' = rmse.prop, 'rmse.med' = rmse.med, 'rmse.gen' = rmse.gen))
   
}
