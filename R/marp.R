#' A function to apply model-averaged renewal process
#' @param data A numeric vector of positive inter-event times.
#' @param t A numeric vector of time points at which log-hazards are evaluated.
#' @param m A positive integer controlling repeated random-start optimizations for
#'   the Gamma, log-logistic, Weibull, and BPT candidate models. The current
#'   implementation seeks `m - 1` acceptable `nlm()` fits for each such model.
#' @param y A time point at which the logit-transformed cumulative event
#'   probability is evaluated.
#' @param which.model Integer identifying a reference or generating model:
#'   1 = Poisson, 2 = Gamma, 3 = log-logistic, 4 = Weibull, 5 = log-normal,
#'   and 6 = Brownian passage time (BPT). This is mainly useful when the
#'   generating model is known, such as in simulations.
#'
#' @return An object of class `marp_fit`: a 19-component named list of
#'   model-specific, selected-model, reference-model, and AIC-weighted estimates.
#' \describe{
#' \item{par1}{First fitted parameter for each model, in the model order listed
#'   under `which.model`. Its meaning is model-specific.}
#' \item{par2}{Second fitted parameter for each model in the same order; `NA`
#'   for the one-parameter Poisson model. Its meaning is model-specific.}
#' \item{logL}{Maximized log-likelihood for each model.}
#' \item{AIC}{Akaike information criterion (AIC)}
#' \item{BIC}{Bayesian information criterion (BIC)}
#' \item{mu_hat}{Estimated mean inter-event time for each model.}
#' \item{pr_hat}{Logit-transformed cumulative event probability at `y` for each model.}
#' \item{haz_hat}{Log-hazard values at `t` for each model.}
#' \item{weights_AIC}{Model weights calculated based on AIC}
#' \item{model_best}{Model selected based on the lowest AIC}
#' \item{mu_best}{Estimated mean obtained from the model with the lowest AIC}
#' \item{pr_best}{Estimated logit event probability from the model with the lowest AIC.}
#' \item{haz_best}{Estimated log-hazards from the model with the lowest AIC.}
#' \item{mu_gen}{Estimated mean obtained from the (true or hypothetical) generating model }
#' \item{pr_gen}{Estimated logit event probability from the reference model.}
#' \item{haz_gen}{Estimated log-hazards from the reference model.}
#' \item{mu_aic}{Estimated mean obtained from model-averaging (using AIC weights)}
#' \item{pr_aic}{AIC-weighted average of the model-specific logit event probabilities.}
#' \item{haz_aic}{AIC-weighted averages of the model-specific log-hazards at `t`.}
#' }
#'
#' @examples
#' set.seed(42)
#' data <-  rgamma(100,3,0.01)
#'
#' # set some parameters
#' m = 10  # repeated random-start optimization setting
#' t = seq(100, 200, by=10)  # time intervals
#' y = 304  # cut-off year for estimating probability
#' which.model <- 2 # specify the generating model
#'
#' # model selection and averaging
#' result <- marp::marp(data, t, m, y, which.model)
#' result
#' summary(result)
#'
#' # Bootstrap confidence intervals delegate to the existing marp_confint()
#' # engine. The original data and bootstrap sizes are supplied explicitly.
#' \dontrun{
#' confint(result, data = data, B = 99, BB = 99)
#' }
#'
#' @export

### Model-Averaged Renewal Processes -------------------------------
marp <- function(data,t,m,y,which.model=1) {
  out <- mapply(c,
                poisson_rp(data, t, y), ## 1. Poisson renewal model
                gamma_rp(data, t, m, y), ## 2. Gamma renewal model
                loglogis_rp(data, t, m, y), ## 3. Log-Logistics renewal model
                weibull_rp(data, t, m, y), ## 4. Weibull renewal model
                lognorm_rp(data, t, y), ## 5. Log-Normal renewal model
                bpt_rp(data, t, m, y), ## 6. BPT renewal model
                USE.NAMES = TRUE, SIMPLIFY = FALSE)
  ## estimated mean, (logit) probability and (log) hazard rates from six renewal models
  mu_hat <- out$mu_hat
  pr_hat <- out$pr_hat
  haz_hat <- matrix(out$haz_hat, length(t), 6)
  aic <- as.numeric(out$AIC)
  finite_aic <- is.finite(aic)
  if (!any(finite_aic)) {
    stop("No candidate renewal model returned a finite AIC.", call. = FALSE)
  }
  ## find the best model (minimum AIC value)
  which_aic <- which.min(aic)
  ## min. AIC
  min_aic <- aic[which_aic]
  ## delta AIC for each model (diff bewteen AIC of each model and the min. AIC)
  delta_aic <- aic - min_aic
  ## compute AIC weights; failed candidates with AIC = Inf receive weight 0
  aic_weight_raw <- exp(-.5 * delta_aic)
  aic_weight_raw[!finite_aic] <- 0
  aic_weight <- aic_weight_raw / sum(aic_weight_raw)
  # which_bic <- which.min(out$BIC)
  # min_bic <-  out$BIC[which_bic]
  # delta_bic <- out$BIC - min_bic
  # bic_weight <- exp(-.5 * delta_bic) / sum(exp(-.5 * delta_bic))
  ## model-averaged mean, (logit) probability and (log) hazard rates
  ## using AIC weights. Failed candidates have zero weight and NA estimates.
  weighted_sum <- function(values, weights) {
    keep <- weights > 0 & is.finite(values)
    sum(values[keep] * weights[keep])
  }
  mu_aic <- weighted_sum(mu_hat, aic_weight)
  pr_aic <- weighted_sum(pr_hat, aic_weight)
  haz_aic <- apply(haz_hat, 1, weighted_sum, weights = aic_weight)
  out1 <-list("weights_AIC" = aic_weight,
              # "weights_BIC" = bic_weight,
              "model_best" = which_aic,
              ## best model is defined as the model with the lowest AIC
              "mu_best" = mu_hat[which_aic],
              "pr_best" = pr_hat[which_aic],
              "haz_best" = haz_hat[, which_aic],
              ## default generating model is a Poisson renewal model
              "mu_gen" = mu_hat[which.model],
              "pr_gen" = pr_hat[which.model],
              "haz_gen" = haz_hat[, which.model],
              ## model-averaging using AIC weights
              "mu_aic" = mu_aic,
              "pr_aic" = pr_aic,
              "haz_aic" = as.numeric(haz_aic))
  result <- append(out, out1)
  .new_marp_fit(
    result,
    call = match.call(),
    nobs = length(data),
    t = t,
    y = y,
    m = m,
    which.model = which.model
  )
}
