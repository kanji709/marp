#' A function to fit Weibull renewal model #' @import weibull_logl
#' @param data input inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param m the number of iterations in nlm
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting Weibull renewal model
#' @export

weibull_rp <- function(data, t, m, y) {
  ## find MLE via numerical optimization (nlm)
  i <- 1
  inits <- NULL
  loop <- NULL
  data1 <- sort(data)
  Fhat <- stats::ppoints(data1)
  k0 <-
    as.numeric(stats::lm(log(-log(1 - Fhat)) ~ log(data1))$coefficients[2])
  while (i < m) {
    tmp_init <-
      cbind(stats::runif(1, 0.8 * mean(data), 1.1 * mean(data)),
            stats::runif(1, 0.8 * k0, 1.1 * k0))
    tryCatch({
      tmp <- stats::nlm(weibull_logl, log(tmp_init), x = data)
      if (tmp$code <= 2.5) {
        eval(parse(text = paste("temp", i, '=tmp', sep = "")))
        loop <- c(loop, tmp$minimum)
        i <- i + 1
      }
    }, error = function(e) {

    })
  }
  index <- which.min(loop)
  mle <- get(paste("temp", index, sep = ""))
  ## log-likelihood, AIC and BIC
  logl <- mle$minimum
  aic <- 2 * logl + 4
  bic <- 2 * logl + 2 * log(length(data))
  ## parameters
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_hat <- par1 * gamma(1 + 1 / par2)
  logitp <- gtools::logit(stats::pweibull(y, par2, par1))
  loghaz <-
    log(stats::dweibull(t, par2, par1) / stats::pweibull(t, par2, par1, lower.tail = FALSE))
  return(list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz))
}
