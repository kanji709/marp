uppert <- function(upper_theta){
  
  ut <- (haz.hat - upper_theta)/sqrt(var.star)
  
  return(ut)
}
upperT <- function(upper_theta){
  
  temp <- prelim$Weights_AIC*sapply(1:6,function(i) sum(Tstar[[i]] <= replicate(B,uppert(upper_theta)[,i]))/B)

  return((sum(temp)-alpha/2)^2)
}

nlm(upperT,runif(1))
