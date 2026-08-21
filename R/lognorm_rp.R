#' A function to fit Log-Normal renewal model
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#'
#' @return An object of class `marp_model_fit` retaining the following eight
#'   named list components:
#' \describe{
#' \item{par1}{Estimated mean (on the log scale) of the Log-Normal model}
#' \item{par2}{Estimated standard deviation (on the log scale)of the Log-Normal model}
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
#' t = seq(100, 200, by=10)  # time intervals
#' y = 304  # cut-off year for estimating probablity
#'
#' # fit Log-Normal renewal model
#' fit <- marp::lognorm_rp(data, t, y)
#' fit
#' summary(fit)
#'
#' @export

lognorm_rp <- function(data, t, y) {
  ## parameters
  par1 <- sum(log(data)) / length(data)
  par2 <- sqrt(sum((log(data) - par1) ^ 2) / length(data))
  ## log-likelihood, AIC and BIC
  logl <- sum(stats::dlnorm(data, par1, par2, log = TRUE))
  aic <- -2 * logl + 4
  bic <- -2 * logl + 2 * log(length(data))
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_hat <- exp(par1 + par2 ^ 2 / 2)
  logitp <- gtools::logit(stats::plnorm(y, par1, par2))
  loghaz <- log(stats::dlnorm(t, par1, par2) / stats::plnorm(t, par1, par2, lower.tail = FALSE))
  out <- list("par1" = par1,"par2" = par2,"logL" = logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz)
  .new_marp_model_fit(
    out,
    model = "Log-normal",
    parameter_names = c("meanlog", "sdlog"),
    call = match.call(),
    nobs = length(data),
    t = t,
    y = y
  )
}
