#--------------------------------------- Parametric Bootstrapping -------------------------------------------------------------------#
#----------------------------------------- Exp Distribution---------------------------#
bstrp <- function(par,n,B,t,R){
  
  ### Exp Distribution
  bstrp <- replicate(B,rexp(n,par[1]))
  
  intro <- sapply(1:B,function(i) prelim.exp(bstrp[,i],t))
  
  par.star <- unlist(intro[1,])
  
  haz.star <- unlist(intro[6,])
  
  var <- apply(haz.star,1,var)
  
  double <- sapply(1:B,function(i) replicate(R,rexp(n,par.star[i])),simplify = F)
  
  intro.double <- sapply(1:B,function(i) apply(double[[i]],2,function(x) prelim.exp(x,t)),simplify = F)
  
  haz.double <- sapply(1:B,function(j) sapply(1:R,function(i) intro.double[[j]][[i]]$hazard),simplify = F)
  
  var.double <- sapply(1:B,function(i) apply(haz.double[[i]],1,function(x) var(x)))
  
}  
#------------------------------------------ Log-Logistic Distribution ------------------------------------------------------#
bstrp.llog <- function(data,B,m){
  
  par <- c(llog.ic(data,m)$Par_1,llog.ic(data,m)$Par_2)
  
  samples <- replicate(B,rllog(length(data),par[1],par[2]))
  
  return(samples)

}  

#--------------------------------- Gamma Distribution ---------------------------------------------------#
bstrp.gamma <- function(data,B,m){
  
  par <- c(gamma.ic(data,m)$Par_1,gamma.ic(data,m)$Par_2)
  
  samples <- replicate(B,rgamma(length(data),par[1],par[2]))
  
  return(samples)
  
}  

#---------------------------------------- Weibull Distribution ------------------------------------------------------#
bstrp.weibull <- function(data,B,m){
  
  par <- c(weibull.ic(data,m)$Par_1,weibull.ic(data,m)$Par_2)
  
  samples <- replicate(B,rweibull(length(data),par[2],par[1]))
  
  return(samples)
  
}  

#------------------------------------------ Log-Normal Distribution ------------------------------------------------------#
bstrp.logn <- function(data,B){
  
  par <- c(logn.ic(data)$Par_1,logn.ic(data)$Par_2)
  
  samples <- replicate(B,rlnorm(length(data),par[1],par[2]))
  
  return(samples)
  
}  

#------------------------------------------ BPT Distribution ------------------------------------------------------#
bstrp.bpt <- function(data,B,m){
  
  par <- c(bpt.ic(data,m)$Par_1,bpt.ic(data,m)$Par_2)
  
  samples <- replicate(B,rinvgauss(length(data),par[1],par[1]/par[2]^2))
  
  return(samples)
  
}  
