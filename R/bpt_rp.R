#' A function to fit BPT renewal model
#' @param data input inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param m the number of iterations in nlm
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting BPT renewal model
#' @export


bpt_rp <- function(data, t, m, y) {
  ## find MLE via numerical optimization (nlm)
  i <- 1
  inits <- NULL
  loop <- NULL
  while (i < m) {
    tmp.init <-
      cbind(stats::runif(1, 0, 2 * mean(data)), sqrt(stats::runif(
        1, 0.5 * sum(data) / sum(abs(data - mean(data))), 1.5 * sum(data) / sum(abs(data - mean(data)))
      )))
    tryCatch({
      tmp <- stats::nlm(bpt_logl, log(tmp.init), x = data)
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
  par2 <- sqrt(exp(mle$estimate[2]))
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_hat <- par1
  logitp <- gtools::logit(statmod::pinvgauss(y, par1, par1 / par2 ^ 2))
  loghaz <-
    log(statmod::dinvgauss(t, par1, par1 / par2 ^ 2) / statmod::pinvgauss(t, par1, par1 / par2 ^ 2, lower.tail = FALSE))
  return(list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz))
}
