#---------------------------- Non-parametric Boostrapping --------------------------------------------#
bstrap <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  for (i in 1:B) {
    #Resample Data
    nsims <- sample(data, size = length(data), replace = TRUE)
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}

#---------------------------- parametric Percentile Boostrapping ----------------------------------------------------#
#------------------------------------------ Exponential Distribution ------------------------------------------------------#
bstrap.exp <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of exp distribution
  par <- exp.ic(data)$Par_1
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rexp(length(data),par)
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}

#------------------------------------------ Log-Logistic Distribution ------------------------------------------------------#
bstrap.llog <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of llog distribution
  par <- c(llog.ic(data,m)$Par_1,llog.ic(data,m)$Par_2)
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rllog(length(data),par[1],par[2])
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}

#------------------------------------------ Gamma Distribution ------------------------------------------------------#
bstrap.llog <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of gamma distribution
  par <- c(gamma.ic(data,m)$Par_1,gamma.ic(data,m)$Par_2)
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rgamma(length(data),par[1],par[2])
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}


#------------------------------------------ Weibull Distribution ------------------------------------------------------#
bstrap.llog <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of weibull distribution
  par <- c(weibull.ic(data,m)$Par_1,weibull.ic(data,m)$Par_2)
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rweibull(length(data),par[2],par[1])
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}

#------------------------------------------ Log-Normal Distribution ------------------------------------------------------#
bstrap.llog <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of log-normal distribution
  par <- c(logn.ic(data)$Par_1,logn.ic(data)$Par_2)
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rlnorm(length(data),par[1],par[2])
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}


#------------------------------------------ BPT Distribution ------------------------------------------------------#
bstrap.llog <- function(data, t, B, m) {
  haz.exp <- NULL
  haz.gamma <- NULL
  haz.llog <- NULL
  haz.weibull <- NULL
  haz.lnorm <- NULL
  haz.bpt <- NULL
  haz.ma <- NULL
  index <- NULL
  
  #Estiamted paramter of log-normal distribution
  par <- c(bpt.ic(data,m)$Par_1,bpt.ic(data,m)$Par_2)
  
  for (i in 1:B) {
    
    #Resample Data
    nsims <- rinvgauss(length(data),par[1],par[1]/par[2]^2)
    
    #Compute estimated paramters, log-Likelihood, AIC, BIc and Weights(both AIC and BIC)
    wei <- weights(nsims, m)
    
    #Extract the parameter(s)
    par.star <- matrix(c(wei$Par_1, wei$Par_2), nrow = 6, ncol = 2)
    
    #Compute the hazard rate of each model with estimated paramters
    hazard.hat <- hazard(par.star, t)
    
    #Store the value of hazard rate for each model
    haz.exp <- rbind(haz.exp, hazard.hat[1, ])
    haz.gamma <- rbind(haz.gamma, hazard.hat[2, ])
    haz.llog <- rbind(haz.llog, hazard.hat[3, ])
    haz.weibull <- rbind(haz.weibull, hazard.hat[4, ])
    haz.lnorm <- rbind(haz.lnorm, hazard.hat[5, ])
    haz.bpt <- rbind(haz.bpt, hazard.hat[6, ])
    
    #Find the single best model with Min. AIC
    ind <-  as.numeric(which.min(wei$AIC))
    
    
    index <- rbind(index, unlist(ind))
    
    #Find the model averaged hazard rate
    haz.ma <- rbind(haz.ma, hazard.hat[ind, ])
    
  }
  output <- list(
    "exp" = haz.exp,
    "gamma" = haz.gamma,
    "llog" = haz.llog,
    "weibull" = haz.weibull,
    "lnorm" = haz.lnorm,
    "bpt" = haz.bpt,
    "ma" = haz.ma,
    "id" = index
  )
  return(output)
}
