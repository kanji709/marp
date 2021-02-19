#' A function to fit Poisson renewal model
#' @param data input inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting Poisson renewal model
#' @export

poisson_rp <- function(data, t, y) {
  ## parameters (one-parameter model)
  par1 <- 1 / mean(data) # lambda of exponential distribution
  par2 <- NA
  ## log-likelihood, AIC and BIC
  logl <- sum(stats::dexp(data, par1, log = TRUE)) # log-likelihood
  aic <- -2 * logl + 2 # AIC
  bic <- -2 * logl + log(length(data)) # BIC
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_hat <- 1 / par1
  logitp <- gtools::logit(stats::pexp(y, par1))
  loghaz <- log(rep(par1, length(t)))
  return(list("par1" = par1,"par2" = par2,"logL" = logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz))
}
