
#------------------------- Parameter estimation, Log-likelihood, AIC and BIC for each Model -------------------------------#

### Exponential Distribution
exp.ic <- function(data) {
  
  # Estimate the parameter with MLE
  par1 <- 1 / mean(data)
  
  par2 <- NA
  
  # Log-likelihood
  ll <- sum(dexp(data, par1, log = TRUE))
  
  # AIC
  aic <- -2 * ll + 2
  
  # BIC
  bic <- -2 * ll + log(length(data))
  
  # results
  out <-
    list(
      "par_1" = par1,
      "par_2" = par2,
      "logL" = ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(out)
  
}



### Log-logistic Distribution

llog.ll <- function(param, x) {
  
  #Define the parameters
  alpha <- exp(param[1])
  
  beta <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dllog(x, beta, alpha, log = T))
  
  return(-logl)

  }


llog.ic <- function(data, m) {
  
  i <- 0
  
  inits <- NULL
  
  loop <- NULL
  
  while (i < m) {
    tmp.init <-
      cbind(runif(1, 0, 2 * log(median(data))), runif(1, 0, 2 * log(mean(data)) / log(median(data))))
    
    tmp <-
      tryCatch(
        nlm(llog.ll, log(tmp.init), x = data),
        error = function(e) {
          list('code' = 3)
        }
      )
    if (tmp$code <= 2.5) {
      loop <- c(loop, tmp$minimum)
      inits <- rbind(inits, tmp.init)
      i = i + 1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(llog.ll, log(inits[index, ]), x = data)
  
  par1 <- exp(mle$estimate[2])
  par2 <- exp(mle$estimate[1])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  out <-
    list(
      "par_1" = par1,
      "par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  
  return(out)
  
}

### Gamma Distribution
gamma.ll <- function(param, x) {
  
  alpha <- exp(param[1])
 
  beta <- exp(param[2])
  
  logl <- sum(dgamma(x, alpha, beta, log = T))
  
  return(-logl)
  
}

gamma.ic <- function(data, m) {
  
  i = 0
  
  inits <- NULL
  
  loop <- NULL
  
  while (i < m) {
    tmp.init <- cbind(
      runif(
        1,
        0.8 * mean(data) ^ 2 / var(data) ,
        1.2 * mean(data) ^ 2 / var(data)
      ),
      runif(1, 0.8 * mean(data) / var(data), 1.2 * mean(data) / var(data))
    )
    
    tmp <-
      tryCatch(
        nlm(gamma.ll, log(tmp.init), x = data),
        error = function(e) {
          list('code' = 3)
        }
      )
    if (tmp$code <= 2.5) {
      loop <- c(loop, tmp$minimum)
      inits <- rbind(inits, tmp.init)
      i = i + 1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(gamma.ll, log(inits[index, ]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  out <-
    list(
      "par_1" = par1,
      "par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(out)
  
  
}




### Weibull Distribution
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
  
  k0 <-
    as.numeric(lm(log(-log(1 - Fhat)) ~ log(da))$coefficients[2])
  
  while (i < m) {
    tmp.init <-
      cbind(runif(1, 0.8 * mean(data), 1.1 * mean(data)), runif(1, 0.8 * k0, 1.1 *
                                                                  k0))
    
    
    tmp <-
      tryCatch(
        nlm(weibull.ll, log(tmp.init), x = data),
        error = function(e) {
          list('code' = 3)
        }
      )
    if (tmp$code <= 2.5) {
      loop <- c(loop, tmp$minimum)
      inits <- rbind(inits, tmp.init)
      i = i + 1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(weibull.ll, log(inits[index, ]), x = data)
  
  #lambda
  par1 <- exp(mle$estimate[1])
  
  #k
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  out <-
    list(
      "par_1" = par1, 
      "par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(out)
  
}


### Log-normal Distribution
logn.ic <- function(data) {
  
  # mu
  par1 <- sum(log(data)) / length(data)
  
  # sigma
  par2 <- sqrt(sum((log(data) - par1) ^ 2) / length(data))
  
  
  ll <- sum(dlnorm(data, par1, par2, log = T))
  
  aic <- -2 * ll + 4
  
  bic <- -2 * ll + 2 * log(length(data))
  
  
  out <-
    list(
      "par_1" = par1,
      "par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(out)
  
}

### BPT Distribution
bpt.ll <- function(param, x){
  
  mu <- exp(param[1])
  
  alpha <- exp(param[2])
  
  
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
  
  while (i < m) {
    tmp.init <-
      cbind(runif(1, 0, 2 * mean(data)), sqrt(runif(
        1, 0.5 * sum(data) / sum(abs(data - mean(data))), 1.5 * sum(data) / sum(abs(data - mean(data)))
      )))
    
    tmp <-
      tryCatch(
        nlm(bpt.ll, log(tmp.init), x = data),
        error = function(e) {
          list('code' = 3)
        }
      )
    if (tmp$code <= 2.5) {
      loop <- c(loop, tmp$minimum)
      inits <- rbind(inits, tmp.init)
      i = i + 1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(bpt.ll, log(inits[index, ]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- sqrt(exp(mle$estimate[2]))
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  out <-
    list(
      "par_1" = par1,
      "par_2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic
    )
  
  return(out)
  
}

#------------------------- Summary of the values of estimated parameters, AIC, BIC, AIC weights, BIC weights and index of "best" model -------------------------------#

weight <- function(data, m) {

  # Estimated parameters, log-likelihood, AIC and BIC
  results <- mapply(
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
  
  # Locate the min. AIC and min. BIC 
  index.aic <- which.min(results$AIC)
  
  index.bic <- which.min(results$BIC)
  
  # Obtain the min. AIC and min. BIC
  min.aic <- results$AIC[index.aic]
  
  min.bic <-  results$BIC[index.bic]
  
  # Calculate the delta AIC and delta BIC
  del.aic <- results$AIC - min.aic
  
  del.bic <- results$BIC - min.bic
  
  # Coumpute model weights using AIC and BIC
  aic.weight <- exp(-.5 * del.aic) / sum(exp(-.5 * del.aic))
  
  bic.weight <- exp(-.5 * del.bic) / sum(exp(-.5 * del.bic))
  
  #Combine all of the results
  weights <-
    list("weights.AIC" = aic.weight,
         "weights.BIC" = bic.weight,
         'ID' = index.aic)
  
  output <- append(results, weights)
  
  
  return(output)
  
}
