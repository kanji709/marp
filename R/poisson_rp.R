#' A function to fit Poisson renewal model
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#'
#' @return An object of class `marp_model_fit` retaining the following eight
#'   named list components:
#' \describe{
#' \item{par1}{Estimated exponential rate parameter of the Poisson renewal model}
#' \item{par2}{N/A, only keep it as a place holder for output formatting purpose}
#' \item{logL}{Maximized log-likelihood}
#' \item{AIC}{Akaike information criterion (AIC)}
#' \item{BIC}{Bayesian information criterion (BIC)}
#' \item{mu_hat}{Estimated mean inter-event time}
#' \item{pr_hat}{Logit-transformed cumulative event probability at `y`}
#' \item{haz_hat}{Log-hazard values at `t`}
#' }
#' @examples
#' set.seed(42)
#' data <-  rgamma(100,3,0.01)
#'
#' # set some parameters
#' t = seq(100, 200, by=10)  # time intervals
#' y = 304  # cut-off year for estimating probablity
#'
#' # fit Poisson renewal model
#' fit <- marp::poisson_rp(data, t, y)
#' fit
#' summary(fit)
#'
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
  out <- list("par1" = par1,"par2" = par2,"logL" = logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz)
  .new_marp_model_fit(
    out,
    model = "Poisson renewal (exponential waiting times)",
    parameter_names = c("rate", "unused"),
    call = match.call(),
    nobs = length(data),
    t = t,
    y = y
  )
}
