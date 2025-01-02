rmse <- function(haz, mu, sig, t) {
  
#------------------------------------------ Hazard rates given by data generated model with true paramters  -----------------------------------------------------------------------#     
  haz.true <- dlnorm(t,mu,sig)/plnorm(t,mu,sig,lower.tail = F)

#------------------------------------------ Variance  -----------------------------------------------------------------------#     
  var <- apply(haz,1,function(x) var(x))
  
#------------------------------------------ Bias  -----------------------------------------------------------------------#     
  bias <- apply(haz,1,function(x) mean(x)) - haz.true
  
#------------------------------------------ R.M.S.E  -----------------------------------------------------------------------#     
  rmse <- sqrt(var+bias^2)
  
  return(list("Var" = var, "Bias"= bias, "RMSE" = rmse))

}
