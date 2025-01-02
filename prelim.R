#   For each model, the functions (prelim. & prelim) compute 
#   1. Estimated parameters;
#   2. Log-likelihood, AIC and BIC; 
#   3. Hazard rates, mean inter-event times (u) and probablity of Pr(X > u).
   
#   For each model, the functions (bstrp. & bstrp) compute 
#   1. Parametric bootstraps;
#   2. Double parametric bootstraps; 
#   3. Estimaed hazard rates from first and second layer bootstrapped resamples;
#   4. Variance of hazard rates from original data sets and fisrt layer bootstrapped resamples;
#   5. Tstars.




# Exp Distribution --------------------------------------------------------

prelim.exp <- function(data, t, mu.true) {
  # Estimate the parameter with MLE
  par1 <- 1 / mean(data)
  
  par2 <- NA
  
  # Log-likelihood
  ll <- sum(dexp(data, par1, log = TRUE))
  
  # AIC
  aic <- -2 * ll + 2
  
  # BIC
  bic <- -2 * ll + log(length(data))
  
  # Hazard rates
  haz <- rep(par1, length(t))
  
  # Mean inter-event times 
  mu <- 1/par1 
  
  # Pr(X > u)
  pr <- pexp(mu.true, par1)
  
  # results
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
}


bstrp.exp <- function(par, n, B, t, R, mean, prob, haz, mu.true) {
  
  bstrp <- replicate(B, rexp(n, par[1]))
  
  intro <- apply(bstrp,2,function(x) prelim.exp(x, t, mu.true))
  
  par.star <- sapply(intro,'[[',1)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)
  
  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(pr.star)
  
  double <- sapply(1:B, function(i) replicate(R, rexp(n, par.star[i])), simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.exp(x, t, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[1]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[1]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,1]) / sqrt(haz.var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'mu.star' = mu.star,
      'pr.star' = pr.star,
      'mu.var.hat' =  mu.var.hat,
      'pr.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mu.var.double' = mu.var.double,
      'pr.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mu.Tstar' = mu.Tstar,
      'pr.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}


# Log-Logistic Distribution -----------------------------------------------

dllog <- function (x, shape = 1, scale = 1, log = FALSE) {
  
  fx <- (shape/scale)*(x/scale)^{shape - 1}/(1 + (x/scale)^shape)^2
  
  if (log) 
    
    return(log(fx))
  
  else return(fx)
  
}

pllog <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE){
  
  Fx <- 1/(1+(q/scale)^{-shape})
  
  if (!lower.tail)
    
    Fx <- 1 - Fx
  
  if (log.p) 
    
    Fx <- log(Fx)
  
  return(Fx)
  
}

llog.ll <- function(param, x) {
  #Define the parameters
  alpha <- exp(param[1])
  
  beta <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dllog(x, beta, alpha, log = T))
  
  return(-logl)
  
}


prelim.llog <- function(data, t, m, mu.true) {
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
      i <- i + 1
    }
  }
  
  index <- which.min(loop)
  
  mle <- nlm(llog.ll, log(inits[index,]), x = data)
  
  # shape
  par1 <- 1/exp(mle$estimate[2])
  
  # scale
  par2 <- log(exp(mle$estimate[1]))
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  haz <- dllog(t, par1, par2) / pllog(t, par1, par2, lower.tail = F)
  
  # mean is undefined if shape parameter par1 < 1
  if (par1 > 1){
    
    mu <- (par2*pi/par1)/sin(pi/par1)
    
  }else{
    
    mu <- mean(data)
      
  }
  
  pr <- pllog(mu.true,par1,par2)
  
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
}

bstrp.llog <- function(par, n, B, t, m, R, mean, prob, haz,mu.true) {
  
  bstrp <- replicate(B, rllog(n, par[3], par[9]))
  
  intro <- apply(bstrp,2,function(x) prelim.llog(x, t, m, mu.true))
  
  par1.star <- sapply(intro,'[[',1)
  
  par2.star <- sapply(intro,'[[',2)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)
  
  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(pr.star)
  
  double <- sapply(1:B, function(i) replicate(R, rllog(n, par1.star[i], par2.star[i])), simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.llog(x, t, m, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[3]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[3]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,3]) / sqrt(haz.var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'mu.star' = mu.star,
      'pr.star' = pr.star,
      'mu.var.hat' =  mu.var.hat,
      'pr.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mu.var.double' = mu.var.double,
      'pr.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mu.Tstar' = mu.Tstar,
      'pr.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}





# Gamma Distribution ------------------------------------------------------

gamma.ll <- function(param, x) {
  alpha <- exp(param[1])
  
  beta <- exp(param[2])
  
  logl <- sum(dgamma(x, alpha, beta, log = T))
  
  return(-logl)
  
}


prelim.gamma <- function(data, t, m, mu.true) {
  i <- 0
  
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
  
  mle <- nlm(gamma.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  haz <-
    dgamma(t, par1, par2) / pgamma(t, par1, par2, lower.tail = F)
  
  mu <- par1/par2
  
  pr <- pgamma(mu.true,par1,par2)
  
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
  
}


bstrp.gamma <- function(par, n, B, t, m, R, mean, prob, haz,mu.true) {
  
  bstrp <- replicate(B, rgamma(n, par[2], par[8]))
  
  intro <- apply(bstrp,2,function(x) prelim.gamma(x, t, m, mu.true))
  
  par1.star <- sapply(intro,'[[',1)
  
  par2.star <- sapply(intro,'[[',2)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)
  
  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(pr.star)
  
  double <- sapply(1:B, function(i) replicate(R, rgamma(n, par1.star[i], par2.star[i])), simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.gamma(x, t, m, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[2]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[2]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,2]) / sqrt(haz.var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'mean.star' = mu.star,
      'prob.star' = pr.star,
      'mean.var.hat' =  mu.var.hat,
      'prob.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mean.var.double' = mu.var.double,
      'prob.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mean.Tstar' = mu.Tstar,
      'prob.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}



# Weibull Distribution ----------------------------------------------------

weibull.ll <- function(param, x) {
  #Define the parameter(s)
  lambda <- exp(param[1])
  
  k <- exp(param[2])
  
  #Log-likelihood
  logl <- sum(dweibull(x, k, lambda, log = T))
  
  return(-logl)
  
}

prelim.weibull <- function(data, t, m, mu.true) {
  i <-  0
  
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
  
  mle <- nlm(weibull.ll, log(inits[index,]), x = data)
  
  # scale
  par1 <- exp(mle$estimate[1])
  
  # shape
  par2 <- exp(mle$estimate[2])
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  haz <- dweibull(t, par2, par1) / pweibull(t, par2, par1, lower.tail = F)
  
  mu <- par1*gamma(1+1/par2)
  
  pr <- pweibull(mu.true,par2,par1) 
    
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
}

bstrp.weibull <- function(par, n, B, t, m, R, mean, prob, haz,mu.true) {
  
  bstrp <- replicate(B, rweibull(n, par[10], par[4]))
  
  intro <- apply(bstrp,2,function(x) prelim.weibull(x, t, m, mu.true))
  
  par1.star <- sapply(intro,'[[',1)
  
  par2.star <- sapply(intro,'[[',2)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)
  
  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(mu.star)
  
  double <- sapply(1:B, function(i) replicate(R, rweibull(n, par2.star[i], par1.star[i])), simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.weibull(x, t, m, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[4]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[4]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,4]) / sqrt(haz.var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'mu.star' = mu.star,
      'pr.star' = pr.star,
      'mu.var.hat' =  mu.var.hat,
      'pr.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mu.var.double' = mu.var.double,
      'pr.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mu.Tstar' = mu.Tstar,
      'pr.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}

# Log-Normal Distribution -------------------------------------------------

prelim.lnorm <- function(data, t, mu.true) {
  # mu
  par1 <- sum(log(data)) / length(data)
  
  # sigma
  par2 <- sqrt(sum((log(data) - par1) ^ 2) / length(data))
  
  
  ll <- sum(dlnorm(data, par1, par2, log = T))
  
  aic <- -2 * ll + 4
  
  bic <- -2 * ll + 2 * log(length(data))
  
  haz <- dlnorm(t, par1, par2) / plnorm(t, par1, par2, lower.tail = F)
  
  mu <- exp(par1+par2^2/2)
  
  pr <- plnorm(mu.true,par1,par2)
  
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
}

bstrp.lnorm <- function(par, n, B, t, R, mean, prob, haz,mu.true) {
  
  bstrp <- replicate(B, rlnorm(n, par[5], par[11]))
  
  intro <- apply(bstrp,2,function(x) prelim.lnorm(x, t, mu.true))
  
  par1.star <- sapply(intro,'[[',1)
  
  par2.star <- sapply(intro,'[[',2)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)
  
  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(pr.star)
  
  double <- sapply(1:B, function(i) replicate(R, rlnorm(n, par1.star[i], par2.star[i])),simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.lnorm(x, t, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[5]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[5]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,5]) / sqrt(haz.var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'mu.star' = mu.star,
      'pr.star' = pr.star,
      'mu.var.hat' =  mu.var.hat,
      'pr.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mu.var.double' = mu.var.double,
      'pr.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mu.Tstar' = mu.Tstar,
      'pr.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}


# Brownian passage time (BPT) Distribution --------------------------------

bpt.ll <- function(param, x) {
  mu <- exp(param[1])
  
  alpha <- exp(param[2])
  
  
  n <- length(x)
  logl <-
    log(1 - 0) + (n / 2) * (log(mu) - log(2 * pi * alpha)) - sum((x - mu) ^
                                                                   2 / (2 * x * mu * alpha)) - 3 * sum(log(x)) / 2
  return(-logl)
  
}


prelim.bpt <- function(data, t, m, mu.true) {
  i <- 0
  
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
  
  mle <- nlm(bpt.ll, log(inits[index,]), x = data)
  
  par1 <- exp(mle$estimate[1])
  
  par2 <- sqrt(exp(mle$estimate[2]))
  
  ll <- mle$minimum
  
  aic <- 2 * ll + 4
  
  bic <- 2 * ll + 2 * log(length(data))
  
  haz <-
    dinvgauss(t, par1, par1 / par2 ^ 2) / pinvgauss(t, par1, par1 / par2 ^
                                                      2, lower.tail = F)
  mu <- par1
  
  pr <- pinvgauss(mu.true,par1,par1/par2^2)
  
  out <-
    list(
      "par1" = par1,
      "par2" = par2,
      "logL" = -ll,
      "AIC" = aic,
      "BIC" = bic,
      "hazard" = haz,
      "mean" = mu,
      "prob" = pr
    )
  
  return(out)
  
}

bstrp.bpt <- function(par, n, B, t, m, R, mean, prob, haz, mu.true) {
  
  bstrp <- replicate(B, rinvgauss(n, par[6], par[6] / par[12] ^ 2))
  
  intro <- apply(bstrp,2,function(x) prelim.bpt(x, t, m, mu.true))
  
  par1.star <- sapply(intro,'[[',1)
  
  par2.star <- sapply(intro,'[[',2)
  
  haz.star <- sapply(intro,'[[',6)
  
  mu.star <- sapply(intro,'[[',7)
  
  pr.star <- sapply(intro,'[[',8)
  
  haz.var.hat <- apply(haz.star, 1, var)

  mu.var.hat <- var(mu.star)
  
  pr.var.hat <- var(pr.star)
    
  double <- sapply(1:B, function(i) replicate(R, rinvgauss(n, par1.star[i], par1.star / par2.star ^ 2)),simplify = "array")
  
  intro.double <- apply(double,c(2,3),function(x) prelim.bpt(x, t, m, mu.true))
  
  haz.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',6))
  
  mu.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',7))
  
  pr.double <- apply(intro.double,c(1,2),function(x) sapply(x,'[[',8))
  
  haz.var.double <- apply(haz.double,c(1,3),var)
  
  mu.var.double <- apply(mu.double,2,var)
  
  pr.var.double <- apply(pr.double,2,var)
  
  mu.Tstar <- (mu.star - mean[6]) / sqrt(mu.var.double)
  
  pr.Tstar <- (pr.star - prob[6]) / sqrt(pr.var.double)
  
  haz.Tstar <- (haz.star - haz[,6]) / sqrt(haz.var.double)
  
  
  return(
    list(
      'haz.star' = haz.star,
      'mu.star' = mu.star,
      'pr.star' = pr.star,
      'mu.var.hat' =  mu.var.hat,
      'pr.var.hat' =  pr.var.hat,
      'haz.var.hat' = matrix(haz.var.hat,200,1),
      'mu.var.double' = mu.var.double,
      'pr.var.double' = pr.var.double,
      'haz.var.double' = haz.var.double,
      'mu.Tstar' = mu.Tstar,
      'pr.Tstar' = pr.Tstar,
      'haz.Tstar' = haz.Tstar
    )
  )
  
}
  


# Summary -----------------------------------------------------------------

prelim <- function(data, t, m, mu.true) {
  
  # Combine the results from each model
  results <- mapply(
    c,
    prelim.exp(data, t, mu.true),
    prelim.gamma(data, t, m, mu.true),
    prelim.llog(data, t, m, mu.true),
    prelim.weibull(data, t, m, mu.true),
    prelim.lnorm(data, t, mu.true),
    prelim.bpt(data, t, m, mu.true),
    USE.NAMES = T,
    SIMPLIFY = F
  )
  
  # Mean inter-event times 
  mu <- results$mean
  
  # Pr(X > mu)
  pr <- results$pr
  
  #Estimated hazard rates
  hazard <- matrix(results$hazard,length(t),6)
  
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
  
  # Extract the mean inter-event times, Pr(X > mu) and hazard from both best and generating model
  vals <-
    list("weights.AIC" = aic.weight,
         "weights.BIC" = bic.weight,
         "ID" = index.aic,
         "mean.best" = mu[index.aic],
         "prob.best" = pr[index.aic],
         "hazard.best" = hazard[,index.aic],
         "mean.gen" = mu[5],
         "prob.gen" = pr[5],
         "hazard.gen" = hazard[,5]
  )
  output <- append(results, vals)
  
  
  return(output)
  
}



bstrp <- function(par, n, B, t, m, R, mean, prob, haz, mu.true) {
  vals <- mapply(
    rbind,
    bstrp.exp(par, n, B, t, R, mean, prob, haz,mu.true),
    bstrp.llog(par, n, B, t, m, R, mean, prob, haz,mu.true),
    bstrp.gamma(par, n, B, t, m, R, mean, prob, haz,mu.true),
    bstrp.weibull(par, n, B, t, m, R, mean, prob, haz,mu.true),
    bstrp.lnorm(par, n, B, t, R, mean, prob, haz,mu.true),
    bstrp.bpt(par, n, B, t, m, R, mean, prob, haz,mu.true)
  )
  
  return(
    list(
      "mu.star" = vals$mu.star,
      "pr.star" = vals$pr.star,
      "haz.star" = array(vals$haz.star, c(length(t), 6, B)),
      "mu.var.hat" = vals$mu.var.hat,
      "pr.var.hat" = vals$pr.var.hat,
      "haz.var.hat" = array(vals$haz.var.hat, c(length(t), 6, 1)),
      "mu.var.double" = vals$mu.var.double,
      "pr.var.double" = vals$pr.var.double,
      "haz.var.double" = array(vals$haz.var.double, c(length(t), 6, B)),
      "mu.Tstar" = vals$mu.Tstar,
      "pr.Tstar" = vals$pr.Tstar,
      "haz.Tstar" = array(vals$haz.Tstar, c(length(t), 6, B))
    )
  )
  
}



# Assesment ---------------------------------------------------------------

rmse <- function(hat, true) {
  
  
  #------------------------------------------ Variance  -----------------------------------------------------------------------#     
  var <- apply(hat,1,function(x) var(x))
  
  #------------------------------------------ Bias  -----------------------------------------------------------------------#     
  bias <- apply(hat,1,function(x) mean(x)) - true
  
  #------------------------------------------ R.M.S.E  -----------------------------------------------------------------------#     
  rmse <- sqrt(var+bias^2)
  
  return(list("SE" = sqrt(var), "Bias"= bias, "RMSE" = rmse))
  
}

cover <- function(lower,upper,true,size){
  
  cover <- sum(true<=upper&true>=lower)/size
  
  error.lower <- sum(true<lower)/size
  
  error.upper <- sum(true>upper)/size
  
  width.lower <- mean((lower - true)/true)
  
  width.upper <- mean((upper - true)/true)
  
  return(list("coverage" = cover, "error.lower" = error.lower, "error.upper" = error.upper,"width.lower" = width.lower,"width.upper" = width.upper))
} 

