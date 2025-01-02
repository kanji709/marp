student <- function(data,mu,sig,m,B,R,t){
#---------------------------------------- Estimated hazard rates from each model ------------------------------------------------------#
  prelim <-  weights(data,m)
  
  par.hat <- c(prelim$Par_1,prelim$Par_2)
  
  haz.hat <- hazard(par.hat,t)
#---------------------------------------- Parametric bootstrapping and Double boostrapping ------------------------------------------------------#
  #---------------------------------------- Exp ------------------------------------------------------#
  exp <- bstrp.exp(data,B)
  
  double.exp <- sapply(1:B,function(i) bstrp.exp(exp[,i],R),simplify = F)
  
  #---------------------------------------- Log-Logistic ------------------------------------------------------#
  llog <- bstrp.llog(data,B,m)
  
  double.llog <- sapply(1:B,function(i) bstrp.llog(llog[,i],R,m),simplify = F)
  
  #---------------------------------------- Gamma ------------------------------------------------------#
  gamma <- bstrp.gamma(data,B,m)
  
  double.gamma <- sapply(1:B,function(i) bstrp.gamma(gamma[,i],R,m),simplify = F)
  
  #---------------------------------------- Weibull ------------------------------------------------------#
  weibull <- bstrp.weibull(data,B,m)
  
  double.weibull <- sapply(1:B,function(i) bstrp.weibull(weibull[,i],R,m),simplify = F)
  
  #---------------------------------------- Log-Normal ------------------------------------------------------#
  logn <- bstrp.logn(data,B)
  
  double.logn <- sapply(1:B,function(i) bstrp.logn(logn[,i],R),simplify = F)
  
  #---------------------------------------- BPT ------------------------------------------------------#
  bpt <- bstrp.bpt(data,B,m)
  
  double.bpt <- sapply(1:B,function(i) bstrp.bpt(bpt[,i],R,m),simplify = F)
  
#---------------------------------- Compute the the estimated hazard rates, variance and Tstar  for each model of each resamples------------------------------------------------------#
  #---------------------------------------- Exp ------------------------------------------------------#
  haz.star.exp <- apply(exp,2,function(x) hazard_exp(t,exp.ic(x)$Par_1))
  
  var.exp <- apply(haz.star.exp,1,function(x) var(x))  
  
  var.double.exp <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.exp[[i]],2,function(x) hazard_exp(t,exp.ic(x)$Par_1)),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.exp <- apply(haz.star.exp,2,function(x) x- as.matrix(haz.hat[,1]))/sqrt(var.double.exp)
  #---------------------------------------- Log-Logistic ------------------------------------------------------#
  haz.star.llog <- apply(llog,2,function(x) hazard_llog(t,c(llog.ic(x,m)$Par_1,llog.ic(x,m)$Par_2)))
  
  var.llog <- apply(haz.star.llog,1,function(x) var(x))  
  
  var.double.llog <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.llog[[i]],2,function(x) hazard_llog(t,c(llog.ic(x,m)$Par_1,llog.ic(x,m)$Par_2))),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.llog <- apply(haz.star.llog,2,function(x) x- as.matrix(haz.hat[,3]))/sqrt(var.double.llog)
  #---------------------------------------- Gamma ------------------------------------------------------#
  haz.star.gamma <- apply(gamma,2,function(x) hazard_gamma(t,c(gamma.ic(x,m)$Par_1,gamma.ic(x,m)$Par_2)))
  
  var.gamma <- apply(haz.star.gamma,1,function(x) var(x))  
  
  var.double.gamma <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.gamma[[i]],2,function(x) hazard_gamma(t,c(gamma.ic(x,m)$Par_1,gamma.ic(x,m)$Par_2))),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.gamma <- apply(haz.star.gamma,2,function(x) x- as.matrix(haz.hat[,2]))/sqrt(var.double.gamma)
  #---------------------------------------- Weibull ------------------------------------------------------#
  haz.star.weibull <- apply(weibull,2,function(x) hazard_weibull(t,c(weibull.ic(x,m)$Par_1,weibull.ic(x,m)$Par_2)))
  
  var.weibull <- apply(haz.star.weibull,1,function(x) var(x))  
  
  var.double.weibull <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.weibull[[i]],2,function(x) hazard_weibull(t,c(weibull.ic(x,m)$Par_1,weibull.ic(x,m)$Par_2))),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.weibull <- apply(haz.star.weibull,2,function(x) x- as.matrix(haz.hat[,4]))/sqrt(var.double.weibull)
  #---------------------------------------- Log-Normal ------------------------------------------------------#
  haz.star.logn <- apply(logn,2,function(x) hazard_logn(t,c(logn.ic(x)$Par_1,logn.ic(x)$Par_2)))
  
  var.logn <- apply(haz.star.logn,1,function(x) var(x))  
  
  var.double.logn <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.logn[[i]],2,function(x) hazard_logn(t,c(logn.ic(x)$Par_1,logn.ic(x)$Par_2))),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.logn <- apply(haz.star.logn,2,function(x) x- as.matrix(haz.hat[,5]))/sqrt(var.double.logn)
  #---------------------------------------- BPT ------------------------------------------------------#
  haz.star.bpt <- apply(bpt,2,function(x) hazard_bpt(t,c(bpt.ic(x,m)$Par_1,bpt.ic(x,m)$Par_2))) 
  
  var.bpt <- apply(haz.star.bpt,1,function(x) var(x))  
  
  var.double.bpt <- sapply(1:B,function(j) apply(sapply(1:B,function(i) apply(double.bpt[[i]],2,function(x) hazard_bpt(t,c(bpt.ic(x,m)$Par_1,bpt.ic(x,m)$Par_2))),simplify = F)[[j]],1,function(y) var(y)))
  
  Tstar.bpt <- apply(haz.star.bpt,2,function(x) x- as.matrix(haz.hat[,6]))/sqrt(var.double.bpt)
#---------------------------------- Calculate lower and upper limitss of the CI for each hazard rates ------------------------------------------------------#
  Tstar <- list(Tstar.exp,Tstar.llog,Tstar.gamma,Tstar.weibull,Tstar.logn,Tstar.bpt)
  
  var.star <- cbind(var.exp,var.llog,var.gamma,var.weibull,var.logn,var.bpt);colnames(var.star) = NULL
}
