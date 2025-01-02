library("statmod")
library("FAdist")


sim <- function(n, mu, sig, P, m, t, B, R, alpha){

  # True mean inter-event times, Pr(X< mu) and hazrd rates ------------------
  
  # Mean inter-event times
  mu.true <- exp(mu + sig^2/2)
   
  # Pr(X < mu) 
  pr.true <- plnorm(mu.true,mu,sig)
   
  # Hazard rates
  haz.true <- dlnorm(t,mu,sig)/plnorm(t,mu,sig,lower.tail = F)
  

# Generate P samples(size n) from log-normal(mu,sig) ----------------------

  pts <- replicate(P,rlnorm(n,mu,sig))
  

# Preliminary analysis ----------------------------------------------------

  # Analyze each simulated sample
  hat <- apply(pts, 2, function(x) prelim(x,t,m, mu.true))
   
  # Estimated parameter values
  par.hat <- array(rbind(sapply(hat,'[[',1),sapply(hat,'[[',2)),c(12,1,P))
  
  # Estimated mean inter-event times
  mu.hat <- sapply(hat,'[[',7)
  
  # Estimated Pr(X < mu)
  pr.hat <- sapply(hat,'[[',8)
  
  # Estimated hazard rates
  haz.hat <- array(sapply(hat,'[[',6),c(length(t),6,P))
  
# Estimates from a single best model (Model Selection) --------------------
  
  # ID of a single best model
  id <- sapply(hat,'[[',11)
  
  # Estimated mean inter-event times
  mu.hat.best <- sapply(hat,'[[',12)
  
  # Estimated Pr(X < mu)
  pr.hat.best <- sapply(hat,'[[',13)
  
  # Estimated hazard rates
  haz.hat.best <- sapply(hat,'[[',14)
  
  
# Estimates from the generating model -------------------------------------
  
  # Mean inter-event times
  mu.hat.gen <- sapply(hat,'[[',15)
  
  # Pr(X < mu)
  pr.hat.gen <- sapply(hat,'[[',16)
  
  # Hazard rates
  haz.hat.gen <- sapply(hat,'[[',17)
  
   
# Model-averaging using AIC weights ---------------------------------------

  # AIC weights 
  weights.aic <- sapply(hat,'[[',9)
  
  # Model-averaved mean inter-event times
  
  mu.aic <- as.matrix(sapply(1:P,function(i) mu.hat[,i]%*%weights.aic[,i]))
  
  # Model-averaved probablity of Pr(X < mu)
  pr.aic <- as.matrix(sapply(1:P,function(i) pr.hat[,i]%*%weights.aic[,i]))
  
  # Model-averaged hazard rates
  haz.aic <- sapply(1:P,function(i) haz.hat[,,i]%*%weights.aic[,i])
   
# Generate and analyze bootstarpped resamples for each sample i=1,.,P ---- 
  
  star <- apply(pts,2,function(x) percent(x,B,t,m,mu.true))
  
# Estimates from the generating model in bootstrapped resamples ----------

  # Median of estimated mean inter-event times
  mu.med.gen <- sapply(star,'[[',2)
  
  # Median of estimated Pr(X < mu)
  pr.med.gen <- sapply(star,'[[',8)
  
  # Median of estimated Hazard rates
  haz.med.gen <- sapply(star,'[[',14)
  
# Percentile bootstrapped C.I for estimates from generating model --------
  
  # Percentile bootstrapped C.I for mean inter-event times 
  mu.percent.gen.lower <- sapply(star,"[[",3)
  
  mu.percent.gen.upper <- sapply(star,"[[",4)
  
  # Percentile bootstrapped C.I for Pr(X < mu)
  pr.percent.gen.lower <- sapply(star,"[[",9)
  
  pr.percent.gen.upper <- sapply(star,"[[",10)
  
  # Percentile bootstrapped C.I for hazard rates
  haz.percent.gen.lower <- sapply(star,"[[",15) 
  
  haz.percent.gen.upper <- sapply(star,"[[",16)
  
# Model-averaging using bootstrapped propotion weights --------------------
  
  # Bootstrapped propotion weights
  weights.prop <- sapply(star,'[[',1)
  
  # Model-averaved mean inter-event times
  mu.prop <- as.matrix(sapply(1:P,function(i) mu.hat[,i]%*%weights.prop[,i]))
  
  # Model-averaved mean Pr(X < mu) 
  pr.prop <- as.matrix(sapply(1:P,function(i) pr.hat[,i]%*%weights.prop[,i]))
  
  # Model-averaved hazard rates  
  haz.prop <- sapply(1:P,function(i) haz.hat[,,i]%*%weights.prop[,i])
   
# Model-averaging from the median of best estimates ------------------------
  
  # Model-averaved mean inter-event times
  mu.med.best <- sapply(star,'[[',5)
  
  # Model-averaved Pr(X < mu)
  pr.med.best <- sapply(star,'[[',11)
  
  # Model-averaved mean inter-event times
  haz.med.best <- sapply(star,'[[',17)

# Model-averaged percentile bootstrapped confidence intervals --------------
  
  # Percentile bootstrapped C.I for mean inter-event times 
  mu.percent.best.lower <- sapply(star,"[[",6)
  
  mu.percent.best.upper <- sapply(star,"[[",7)
  
  # Percentile bootstrapped C.I for Pr(X < mu)
  pr.percent.best.lower <- sapply(star,"[[",12)
  
  pr.percent.best.upper <- sapply(star,"[[",13) 
  
  # Percentile bootstrapped C.I for hazard rates
  haz.percent.best.lower <- sapply(star,"[[",18)
  
  haz.percent.best.upper <- sapply(star,"[[",19)  

# Studentized bootstrapped confidence intervals ---------------------------
  
  # Studentized bootstrapped C.I from generating model
  student.gen.index <- sapply(1:P, function(i) replace(matrix(0,6,P)[,i], 5, 1))
  
  student.gen.ci <- sapply(1:P,function(i) student(par.hat[,,i],n,B,t,m,R,mu.hat[,i],pr.hat[,i],haz.hat[,,i],student.gen.index[,i], alpha, mu.true))
  
  # Mean inter-event times
  mu.student.gen.lower <- matrix(unlist(student.gen.ci[1,]),1,P)
  
  mu.student.gen.upper <- matrix(unlist(student.gen.ci[2,]),1,P)
  
  # Pr(X < mu)
  pr.student.gen.lower <- matrix(unlist(student.gen.ci[3,]),1,P)
  
  pr.student.gen.upper <- matrix(unlist(student.gen.ci[4,]),1,P)
  
  # Hazard rates
  haz.student.gen.lower <- matrix(unlist(student.gen.ci[5,]),length(t),P)
  
  haz.student.gen.upper <- matrix(unlist(student.gen.ci[6,]),length(t),P)
  
  
 
  # Studentized bootstrapped C.I from best models
  student.best.index <- sapply(1:P, function(i) replace(matrix(0,6,P)[,i], id[i], 1))
  
  student.best.ci <- sapply(1:P,function(i) student(par.hat[,,i],n,B,t,m,R,mu.hat[,i],pr.hat[,i],haz.hat[,,i],student.best.index[,i], alpha, mu.true))
  
  # Mean inter-event times
  mu.student.best.lower <- matrix(unlist(student.best.ci[1,]),1,P)
  
  mu.student.best.upper <- matrix(unlist(student.best.ci[2,]),1,P)
  
  # Pr(X < mu)
  pr.student.best.lower <- matrix(unlist(student.best.ci[3,]),1,P)
  
  pr.student.best.upper <- matrix(unlist(student.best.ci[4,]),1,P)
  
  # Hazard rates
  haz.student.best.lower <- matrix(unlist(student.best.ci[5,]),length(t),P)
  
  haz.student.best.upper <- matrix(unlist(student.best.ci[6,]),length(t),P)
 
  # Model-averaged studentized bootstrapped C.I
  student.ma.ci <- sapply(1:P,function(i) student(par.hat[,,i],n,B,t,m,R,mu.hat[,i],pr.hat[,i],haz.hat[,,i],weights.aic[,i], alpha, mu.true))
  
  # Mean inter-event times
  mu.student.ma.lower <- matrix(unlist(student.ma.ci[1,]),1,P)
  
  mu.student.ma.upper <- matrix(unlist(student.ma.ci[2,]),1,P)
  
  # Pr(X < mu)
  pr.student.ma.lower <- matrix(unlist(student.ma.ci[3,]),1,P)
  
  pr.student.ma.upper <- matrix(unlist(student.ma.ci[4,]),1,P)
  
  # Hazard rates
  haz.student.ma.lower <- matrix(unlist(student.ma.ci[5,]),length(t),P)
  
  haz.student.ma.upper <- matrix(unlist(student.ma.ci[6,]),length(t),P)
  
  
# Var, bias and RMSE ------------------------------------------------------
  
  ### Mdel Selcetion (Best model)
  
  # Mean inter-event times
  
  rmse.mu.best <- rmse(t(mu.hat.best),mu.true)
  
  # Pr(X < mu)
  rmse.pr.best <- rmse(t(pr.hat.best),pr.true)
  
  # Hazard rates
  rmse.haz.best <- rmse(haz.hat.best,haz.true)
  
  ### Estimates from generating model
  
  # Mean inter-event times
  
  rmse.mu.gen <- rmse(t(mu.hat.gen),mu.true)
  
  # Pr(X < mu)
  rmse.pr.gen <- rmse(t(pr.hat.gen),pr.true)
  
  # Hazard rates
  rmse.haz.gen <- rmse(haz.hat.gen,haz.true)
  
  
  ### Model-averaged estimates using AIC weights
  
  # Mean inter-event times
  
  rmse.mu.aic <- rmse(t(mu.aic),mu.true)
  
  # Pr(X < mu)
  rmse.pr.aic <- rmse(t(pr.aic),pr.true)
  
  # Hazard rates
  rmse.haz.aic <- rmse(haz.aic,haz.true)
  
  ### Model-averaged hazard rates with propotion weights
  
  # Mean inter-event times
  
  rmse.mu.prop <- rmse(t(mu.prop),mu.true)
  
  # Pr(X < mu)
  rmse.pr.prop <- rmse(t(pr.prop),pr.true)
  
  # Hazard rates
  rmse.haz.prop <- rmse(haz.prop,haz.true)
  
  ### Median of estimates in generating model in bootstarpps
  
  # Mean inter-event times
  
  rmse.mu.med.gen <- rmse(matrix(mu.med.gen,1,P),mu.true)
  
  # Pr(X < mu)
  rmse.pr.med.gen <- rmse(matrix(pr.med.gen,1,P),pr.true)
  
  # Hazard rates
  rmse.haz.med.gen <- rmse(haz.med.gen,haz.true)
  
  
  ### Model-averaged hazard rates as median of bootstrapped estimates 
  
  # Mean inter-event times
  
  rmse.mu.med.best <- rmse(matrix(mu.med.best,1,P),mu.true)
  
  # Pr(X < mu)
  rmse.pr.med.best <- rmse(matrix(pr.med.best,1,P),pr.true)
  
  # Hazard rates
  rmse.haz.med.best <- rmse(haz.med.best,haz.true)
  
# Coverage rates ----------------------------------------------------------
  
  ### Percentile bootstrapped C.I for estimates from generating model --------
  
  # Percentile bootstrapped C.I for mean inter-event times 
  coverage.mu.percent.gen <- cover(mu.percent.gen.lower[i],mu.percent.gen.upper[i],mu.true,P)
  
  # Percentile bootstrapped C.I for Pr(X < mu)
  coverage.pr.percent.gen <- cover(pr.percent.gen.lower[i],pr.percent.gen.upper[i],pr.true,P)
  
  # Percentile bootstrapped C.I for hazard rates
  coverage.haz.percent.gen <- sapply(1:length(t),function(i) cover(haz.percent.gen.lower[i,],haz.percent.gen.upper[i,],haz.true,P))
  
  ### Model-averaged percentile bootstrapped confidence intervals
  
  # Percentile bootstrapped C.I for mean inter-event times 
  coverage.mu.percent.best <- cover(mu.percent.best.lower[i],mu.percent.best.upper[i],mu.true,P)
  
  # Percentile bootstrapped C.I for Pr(X < mu)
  coverage.pr.percent.best <- cover(pr.percent.best.lower[i],pr.percent.best.upper[i],pr.true,P)

  # Percentile bootstrapped C.I for hazard rates
  coverage.haz.percent.best <- sapply(1:length(t),function(i) cover(haz.percent.best.lower[i,],haz.percent.best.upper[i,],haz.true,P))

  
  
  ### Studentized bootstrapped C.I for estimates from generating model
  # Studentized bootstrapped C.I for mean inter-event times
  coverage.mu.student.gen <- cover(mu.student.gen.lower[i],mu.student.gen.upper[i],mu.true,P)
    
  # Studentized bootstrapped C.I for Pr(X < mu)
  coverage.pr.student.gen <- cover(pr.student.gen.lower[i],pr.student.gen.upper[i],pr.true,P)
  
  # Studentized bootstrapped C.I for hazard rates
  coverage.haz.student.gen <- sapply(1:length(t),function(i) cover(haz.student.gen.lower[i,],haz.student.gen.upper[i,],haz.true,P))
  
  ### Studentized bootstrapped C.I for estimates from best models
  # Studentized bootstrapped C.I for mean inter-event times
  coverage.mu.student.best <- sapply(1:length(t),function(i) cover(mu.student.best.lower[i],mu.student.best.upper[i],mu.true,P))
  
  # Studentized bootstrapped C.I for Pr(X < mu)
  coverage.pr.student.best <- sapply(1:length(t),function(i) cover(pr.student.best.lower[i],pr.student.best.upper[i],pr.true,P))
  
  # Studentized bootstrapped C.I for hazard rates
  coverage.haz.student.best <- sapply(1:P,function(i) cover(haz.student.best.lower[,i],haz.student.best.upper[,i],haz.true,length(t)))
  
  ### Model-averaged studentized bootstrapped C.I
  
  # Studentized bootstrapped C.I for mean inter-event times
  coverage.mu.student.ma <- cover(mu.student.ma.lower[i],mu.student.ma.upper[i],mu.true,P)
  
  # Studentized bootstrapped C.I for Pr(X < mu)
  coverage.pr.student.ma <- cover(pr.student.ma.lower[i],pr.student.ma.upper[i],pr.true,P)
  
  # Studentized bootstrapped C.I for hazard rates
  coverage.haz.student.ma <- sapply(1:length(t),function(i) cover(haz.student.ma.lower[,i],haz.student.ma.upper[,i],haz.true,P))
  
  
  return(
    list(
      'rmse.mu.gen' = rmse.mu.gen,
      'rmse.pr.gen' = rmse.pr.gen,
      'rmse.haz.gen' = rmse.haz.gen,
      'rmse.mu.best' = rmse.mu.best,
      'rmse.pr.best' = rmse.pr.best,
      'rmse.haz.best' = rmse.haz.best,
      'rmse.mu.aic' = rmse.mu.aic,
      'rmse.pr.aic' = rmse.pr.aic,
      'rmse.haz.aic' = rmse.haz.aic,
      'rmse.mu.prop' = rmse.mu.prop,
      'rmse.pr.prop' = rmse.pr.prop,
      'rmse.haz.prop' = rmse.haz.prop,
      'rmse.mu.med.gen' = rmse.mu.med.gen,
      'rmse.pr.med.gen' = rmse.pr.med.gen,
      'rmse.haz.med.gen' = rmse.haz.med.gen,
      'rmse.mu.med.best' = rmse.mu.med.best,
      'rmse.pr.med.best' = rmse.pr.med.best,
      'rmse.haz.med.best' = rmse.haz.med.best,
      'coverage.mu.percent.gen' = coverage.mu.percent.gen,
      'coverage.pr.percent.gen' = coverage.pr.percent.gen,
      'coverage.haz.percent.gen' = coverage.haz.percent.gen,
      'coverage.mu.percent.best' = coverage.mu.percent.best,
      'coverage.pr.percent.best' = coverage.pr.percent.best,
      'coverage.haz.percent.best' = coverage.haz.percent.best,
      'coverage.mu.student.gen' = coverage.mu.student.gen,
      'coverage.pr.student.gen' = coverage.pr.student.gen,
      'coverage.haz.student.gen' = coverage.haz.student.gen,
      'coverage.mu.student.best' = coverage.mu.student.best,
      'coverage.pr.student.best' = coverage.pr.student.best,
      'coverage.haz.student.best' = coverage.haz.student.best,
      'coverage.mu.student.ma' = coverage.mu.student.ma,
      'coverage.pr.student.ma' = coverage.pr.student.ma,
      'coverage.haz.student.ma' = coverage.haz.student.ma
    )
  )
  
}





