#' A function to fit BPT renewal model
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param m A positive integer controlling repeated random-start optimizations;
#'   the current implementation seeks `m - 1` acceptable `nlm()` fits and caps
#'   the number of attempted starts at `20 * m`.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#'
#' @return An object of class `marp_model_fit` retaining the following eight
#'   named list components:
#' \describe{
#' \item{par1}{Estimated mean parameter (mu) of the BPT model}
#' \item{par2}{Estimated aperiodicity parameter (alpha) of the BPT model}
#' \item{logL}{Maximized log-likelihood}
#' \item{AIC}{Akaike information criterion (AIC)}
#' \item{BIC}{Bayesian information criterion (BIC)}
#' \item{mu_hat}{Estimated mean inter-event time}
#' \item{pr_hat}{Logit-transformed cumulative event probability at `y`}
#' \item{haz_hat}{Log-hazard values at `t`}
#' }
#' If the requested plausible fits cannot be obtained within the existing
#' attempt limit, the function warns and returns the same components with
#' unavailable estimates set to `NA` and information criteria set to `Inf`.
#'
#' @examples
#' set.seed(42)
#' data <-  rgamma(30,3,0.01)
#'
#' # set some parameters
#' m <- 10  # repeated random-start optimization setting
#' t <- seq(100, 200, by=10)  # time intervals
#' y <- 304  # cut-off year for estimating probablity
#'
#' # fit BPT renewal model
#' fit <- marp::bpt_rp(data, t, m, y)
#' fit
#' summary(fit)
#'
#' @export


bpt_rp <- function(data, t, m, y) {
  fit_call <- match.call()
  failed_fit <- function(message) {
    warning(message, call. = FALSE)
    out <- list(
      "par1" = NA_real_,
      "par2" = NA_real_,
      "logL" = NA_real_,
      "AIC" = Inf,
      "BIC" = Inf,
      "mu_hat" = NA_real_,
      "pr_hat" = NA_real_,
      "haz_hat" = rep(NA_real_, length(t))
    )
    .new_marp_model_fit(
      out,
      model = "Brownian passage time",
      parameter_names = c("mean", "alpha"),
      call = fit_call,
      nobs = length(data),
      t = t,
      y = y,
      status = "failed"
    )
  }

  ## find MLE via numerical optimization (nlm)
  i <- 1
  inits <- NULL
  loop <- NULL
  attempts <- 0
  max_attempts <- m * 20
  while (i < m) {
    attempts <- attempts + 1
    if (attempts > max_attempts) {
      return(failed_fit(paste0(
        "bpt_rp: could not find ", m,
        " plausible BPT fits after ", max_attempts, " attempts"
      )))
    }
    tmp.init <-
      cbind(stats::runif(1, 0, 2 * mean(data)), sqrt(stats::runif(
        1, 0.5 * sum(data) / sum(abs(data - mean(data))), 1.5 * sum(data) / sum(abs(data - mean(data)))
      )))
    tryCatch({
      tmp <- stats::nlm(bpt_logl, log(tmp.init), x = data)
      if (tmp$code <= 2.5) {
        par1_candidate <- exp(tmp$estimate[1])
        if (is.finite(tmp$minimum) && is.finite(par1_candidate) && par1_candidate < 10 * max(data)) {
          eval(parse(text = paste("temp", i, '=tmp', sep = "")))
          loop <- c(loop, tmp$minimum)
          i <- i + 1
        }
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
  out <- list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz)
  .new_marp_model_fit(
    out,
    model = "Brownian passage time",
    parameter_names = c("mean", "alpha"),
    call = fit_call,
    nobs = length(data),
    t = t,
    y = y
  )
}
