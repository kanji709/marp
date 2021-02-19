#' A function to fit Gamma renewal model
#' @param data input inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param m the number of iterations in nlm
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting Gamma renewal model
#' @export

gamma_rp <- function(data, t, m, y) {
  ## find MLE via numerical optimization (nlm)
  i <- 1
  inits <- NULL
  loop <- NULL
  while (i < m) {
    tmp_init <- cbind(
      stats::runif(
        1,
        0.8 * mean(data) ^ 2 / stats::var(data) ,
        1.2 * mean(data) ^ 2 / stats::var(data)
      ),
      stats::runif(1, 0.8 * mean(data) / stats::var(data), 1.2 * mean(data) / stats::var(data))
    )
    tryCatch({
      tmp <- stats::nlm(gamma_logl, log(tmp_init), x = data)
      if (tmp$code <= 2.5) {
        eval(parse(text = paste("tmp", i, '=tmp', sep = "")))
        loop <- c(loop, tmp$minimum)
        i <- i + 1
      }
    }, error = function(e) {

    })
  }
  index <- which.min(loop)
  mle <- get(paste("tmp", index, sep = ""))
  ## log-likelihood, AIC and BIC
  logl <- mle$minimum
  aic <- 2 * logl + 4
  bic <- 2 * logl + 2 * log(length(data))
  ## parameters
  par1 <- exp(mle$estimate[1])
  par2 <- exp(mle$estimate[2])
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_hat <- par1 / par2
  logitp <- gtools::logit(stats::pgamma(y, par1, par2))
  loghaz <- log(stats::dgamma(t, par1, par2) / stats::pgamma(t, par1, par2, lower.tail = F))
  return(list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz))
}
