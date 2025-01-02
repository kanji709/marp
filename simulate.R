simulate <- function(n,mu,sig,p,B,m,t){

#------------------------------------------ Simulate p datasets(size n) from log-normal(mu,sig) ------------------------------------------------------------------#  
   pts <- replicate(p,rlnorm(n,mu,sig))
  
#------------------------------------------ Preliminary Analysis -----------------------------------------------------------------------------------------------------#       
   prelim <- apply(pts,2,function(x) weight(x,m))

### Estimated parameters for each generated data sets
   par.hat <- rbind(sapply(prelim,'[[',1),sapply(prelim,'[[',2))

### Estimated hazard rates from each model    
   haz.hat <- 
#------------------------------------------ Model-averaged hazard rates and Percentile bootstrapped confidence intervals ---------------------------------------------#   
   out <- apply(pts,2,function(x) percentile(x,B,m,t))
   
#------------------------------------------ Calculate the variance and bias for each estimates -----------------------------------------------------------------------#   
   
   #Estimated rates from all different approache
   haz <- sapply(1:4,function(i) sapply(out,'[[',i,simplify = T),simplify = F)
   
   #Model-averaged hazard rates with AIC weights
   rmse.aic <- rmse(haz[[1]],mu,sig,t)
   
   #Model-averaged hazard rates with propotion weights
   rmse.prop <- rmse(haz[[2]],mu,sig,t)
   
   #Model-averaged hazard rates as median of bootstrapped estimates 
   rmse.med <- rmse(haz[[3]],mu,sig,t)
   
   #Estimatd hazard rates from data generation model
   rmse.gen <- rmse(haz[[4]],mu,sig,t)

   return(list('results' = out,'rmse.aic' = rmse.aic, 'rmse.prop' = rmse.prop, 'rmse.med' = rmse.med, 'rmse.gen' = rmse.gen))
   
}
