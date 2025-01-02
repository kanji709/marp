#-------------------------------------------- Non-parametric Bootstrapping ---------------------------------#
bstrp.np <- function(data,B){
  
  samples <- replicate(B,sample(data,size = length(data),replace = TRUE))
    
  return(samples)
    
}

#--------------------------------------- Parametric Bootstrapping -------------------------------------------------------------------#
#----------------------------------------- Exp Distribution---------------------------#
bstrp.exp <- function(data,B){
  
  par <- exp.ic(data)$Par_1
  
  samples <- replicate(B,rexp(length(data),par))
  
  return(samples)
  
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