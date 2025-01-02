#--------------- Estimated paramters, Log-likelihood, AIC and BIC for each Model -------------------------------#

#------------------------ Exponential Distribution ----------------------------#
exp.ic <- function(data) {
  #Estimate parameter value with MLE
  par1 <- 1 / mean(data)
  par2 <- NA
  
  #Log-likelihood
  ll <- sum(dexp(data, par1, log = TRUE))
  
  #AIC
  aic <- -2 * ll + 2
  
  #BIC
  bic <- -2 * ll + log(length(data))
  
  #Combine results
  ic <-
    list(
      "Par_1" = par1,
      "Par_2" = par2,
      "logL" = ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(ic)
  
}



#-------------------------Log-logistic Distribution-------------------------------------------#

llog.ll <- function(param, x) {
  #Define the parameter(s)
  alpha <- exp(param[1])
  
  beta <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dllog(x, beta, alpha, log = T))
  
  return(-logl)
}


llog.ic <- function(data, m) {
  
  i = 0
  
  inits <- NULL
  loop <- NULL
  
  while (i < m){
    tmp.init <- cbind(runif(1, 0, 2 * log(median(data))),runif(1,0,2 * log(mean(data)) / log(median(data))))
    
    tmp <- tryCatch(nlm(llog.ll, log(tmp.init), x = data),error = function(e){list('code' = 3)})
    if (tmp$code <= 2.5){
      
      loop <- c(loop,tmp$minimum)
      inits <- rbind(inits,tmp.init)
      i = i +1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(llog.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[2])
  par2 <- exp(mle$estimate[1])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  ic <-
    list(
      "Par_1" = par1,
      "Par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  
  return(ic)
  
}

#--------------------------------- Gamma Distribution ---------------------------------------------------#
gamma.ll <- function(param, x) {
  
  #Define the parameter(s)
  alpha <- exp(param[1])
  beta <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dgamma(x, alpha, beta, log = T))
  
  return(-logl)
  
}

gamma.ic <- function(data, m) {
  
  i = 0
  
  inits <- NULL
  loop <- NULL
  
  while (i < m){
    tmp.init <- cbind(
      runif(1,0.8*mean(data) ^ 2 / var(data) ,1.2*mean(data) ^ 2 / var(data) ),
      runif(1, 0.8*mean(data) / var(data), 1.2 * mean(data) / var(data))
    )
    
    tmp <- tryCatch(nlm(gamma.ll, log(tmp.init), x = data),error = function(e){list('code' = 3)})
    if (tmp$code <= 2.5){
      
      loop <- c(loop,tmp$minimum)
      inits <- rbind(inits,tmp.init)
      i = i +1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(gamma.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  ic <-
    list(
      "Par_1" = par1,
      "Par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(ic)
  
  
}




#--------------------------------------------- Weibull Distribution ----------------------------------------------------------------#
weibull.ll <- function(param, x) {
  #Define the parameter(s)
  lambda <- exp(param[1])
  k <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dweibull(x, k, lambda, log = T))
  
  return(-logl)
  
}

weibull.ic <- function(data, m) {
  
  i = 0
  
  inits <- NULL
  loop <- NULL
  
  #Estimate k0 from the weibull plot
  da <- sort(data)
  Fhat <- ppoints(da)
  k0 <- as.numeric(lm(log(-log(1 - Fhat)) ~ log(da))$coefficients[2])
  
  while (i < m){
    tmp.init <- cbind(runif(1, 0.8* mean(data),1.1* mean(data)), runif(1, 0.8*k0,1.1*k0))
    
    
    tmp <- tryCatch(nlm(weibull.ll, log(tmp.init), x = data),error = function(e){list('code' = 3)})
    if (tmp$code <= 2.5){
      
      loop <- c(loop,tmp$minimum)
      inits <- rbind(inits,tmp.init)
      i = i +1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(weibull.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  ic <-
    list(
      "Par_1" = par1,
      #lambda
      "Par_2" = par2,
      #k
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(ic)
  
}


#------------------------------------------- Log-normal Distribution --------------------------------------------------------------#
logn.ic <- function(data) {
  
  # mu
  par1 <- sum(log(data)) / length(data)
  
  # sigma
  par2 <- sqrt(sum((log(data) - par1) ^ 2) / length(data))
  
  
  ll <- sum(dlnorm(data, par1, par2, log = T))
  
  aic <- -2 * ll + 4
  
  bic <- -2 * ll + 2 * log(length(data))
  
  
  ic <-
    list(
      "Par_1" = par1,
      "Par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(ic)
  
}

#-------------------------------------------- BPT Distribution ------------------------------------------------------#
bpt.ll <- function(param, x) {
  #Define the parameter(s)
  mu <- exp(param[1])
  alpha <- exp(param[2])
  
  
  #Log-likelihood
  n <- length(x)
  logl <-
    log(1 - 0) + (n / 2) * (log(mu) - log(2 * pi * alpha)) - sum((x - mu) ^
                                                                   2 / (2 * x * mu * alpha)) - 3 * sum(log(x)) / 2
  return(-logl)
  
}


bpt.ic <- function(data, m) {
  
  i = 0
  
  inits <- NULL
  loop <- NULL
  
  while (i < m){
    tmp.init <-
      cbind(runif(1, 0, 2 * mean(data)),sqrt(runif(1, 0.5 * sum(data) / sum(abs(data - mean(data))), 1.5 * sum(data) / sum(abs(data - mean(data))))))
    
    tmp <- tryCatch(nlm(bpt.ll, log(tmp.init), x = data),error = function(e){list('code' = 3)})
    if (tmp$code <= 2.5){
      
      loop <- c(loop,tmp$minimum)
      inits <- rbind(inits,tmp.init)
      i = i +1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(bpt.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- sqrt(exp(mle$estimate[2]))
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  ic <-
    list(
      "Par_1" = par1,
      "Par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(ic)
  
}

#------------------------------------- AIC and BIC Weights ----------------------------------------------------------# 
weights <- function(data, m) {
  
  weights <- list(NULL)
  
  #Using ICs functions to compute parameters, logLikelihood, AIC and BIC. Also, combine the results
  ic <- mapply(
    c,
    exp.ic(data),
    gamma.ic(data, m),
    llog.ic(data, m),
    weibull.ic(data, m),
    logn.ic(data),
    bpt.ic(data, m),
    USE.NAMES = T,
    SIMPLIFY = F
  )
  
  #Locate the min. AIC
  index.aic <- which.min(ic$AIC)
  
  #Find the min. AIC
  min.aic <- ic$AIC[index.aic]
  
  #Calculate the delta AIC
  del.aic <- ic$AIC - min.aic
  
  #Model weights(AIC)
  aic.weight <- exp(-.5 * del.aic) / sum(exp(-.5 * del.aic))
  
  #Same procedure for BIC
  index.bic <- which.min(ic$BIC)
  
  min.bic <-  ic$BIC[index.bic]
  
  del.bic <- ic$BIC - min.bic
  
  bic.weight <- exp(-.5 * del.bic) / sum(exp(-.5 * del.bic))
  
  #Combine allof the results
  results <-
    list("Weights_AIC" = aic.weight, "Weights_BIC" = bic.weight)
  
  weights <- append(ic, results)
  
  
  return(weights)
  
}
