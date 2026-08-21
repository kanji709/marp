#' A function to fit Weibull renewal model
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param m A positive integer controlling repeated random-start optimizations;
#'   the current implementation seeks `m - 1` acceptable `nlm()` fits.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#'
#' @return An object of class `marp_model_fit` retaining the following eight
#'   named list components:
#' \describe{
#' \item{par1}{Estimated scale parameter of the Weibull model}
#' \item{par2}{Estimated shape parameter of the Weibull model}
#' \item{logL}{Maximized log-likelihood}
#' \item{AIC}{Akaike information criterion (AIC)}
#' \item{BIC}{Bayesian information criterion (BIC)}
#' \item{mu_hat}{Estimated mean inter-event time}
#' \item{pr_hat}{Logit-transformed cumulative event probability at `y`}
#' \item{haz_hat}{Log-hazard values at `t`}
#' }
#'
#' @examples
#' set.seed(42)
#' data <-  rgamma(100,3,0.01)
#'
#' # set some parameters
#' m = 10  # repeated random-start optimization setting
#' t = seq(100, 200, by=10)  # time intervals
#' y = 304  # cut-off year for estimating probablity
#'
#' # fit Weibull renewal model
#' fit <- marp::weibull_rp(data, t, m, y)
#' fit
#' summary(fit)
#'
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
  out <- list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz)
  .new_marp_model_fit(
    out,
    model = "Weibull",
    parameter_names = c("scale", "shape"),
    call = match.call(),
    nobs = length(data),
    t = t,
    y = y
  )
}
