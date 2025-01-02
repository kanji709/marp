rmse <- function(est, mu, sig, t) {
  
  var <- rep(NA,length(t))
  
  bias <- rep(NA,length(t))
  
  rmse <- rep(NA, length(t))
  
  for (i in 1:length(t)) {
    
    var[i] <- var(est[i, ])
    
    bias[i] <- mean(est[i, ]) - hazard_logn(t[i], c(mu, sig))
    
    rmse[i] <- sqrt(var[i] + bias[i] ^ 2)
  
    }
  
  return(data.frame("var" = var, "Bias"= bias, "RMSE" = rmse))

}
