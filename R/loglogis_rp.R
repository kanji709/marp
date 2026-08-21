#' A function to fit Log-Logistics renewal model
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param m A positive integer controlling repeated random-start optimizations;
#'   the current implementation seeks `m - 1` acceptable `nlm()` fits.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#'
#' @return An object of class `marp_model_fit` retaining the following eight
#'   named list components:
#'
#' \describe{
#' \item{par1}{Estimated shape parameter of the Log-Logistics model}
#' \item{par2}{Estimated scale parameter of the Log-Logistics model}
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
#' # fit Log-Logistic renewal model
#' fit <- marp::loglogis_rp(data, t, m, y)
#' fit
#' summary(fit)
#'
#' @export

loglogis_rp <- function(data, t, m, y) {
  ## find MLE via numerical optimization (nlm)
  i <- 1
  inits <- NULL
  loop <- NULL
  data_log_median <- log(stats::median(data))
  data_log_mean <- log(mean(data))
  while (i < m) {
    tmp_init <-
      cbind(stats::runif(1, 0, 2 * data_log_median), stats::runif(1, 0, 2 * data_log_mean / data_log_median))
    tryCatch({
      tmp <- stats::nlm(loglogis_logl, log(tmp_init), x = data)
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
  par1 <- exp(mle$estimate[2])
  par2 <- exp(mle$estimate[1])

  if (par1 > 1) { ## mean is undefined if par1 (shape parameter) < 1
    mu_hat <- (par2 * pi / par1) / sin(pi / par1)
  } else{ ## set to be sample mean
    mu_hat <- mean(data)
  }
  logitp <- gtools::logit(pllog(y, par1, par2))
  loghaz <- log(dllog(t, par1, par2) / pllog(t, par1, par2, lower.tail = FALSE))
  out <- list("par1" = par1,"par2" = par2,"logL" = -logl,"AIC" = aic,"BIC" = bic,"mu_hat" = mu_hat,"pr_hat" = logitp,"haz_hat" = loghaz)
  .new_marp_model_fit(
    out,
    model = "Log-logistic",
    parameter_names = c("shape", "scale"),
    call = match.call(),
    nobs = length(data),
    t = t,
    y = y
  )
}
