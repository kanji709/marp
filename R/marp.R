#' A function to apply model-averaged renewal process
#' @param data input inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param m the number of iterations in nlm
#' @param y user-specified time point (used to compute time-to-event probability)
#' @param which.model user-specified genearting (or true underlying if known) model
#' @return returns list of estimates obtained from different renewal processes and after applying model-averaging
#'
#' @examples
#' # load example dataset (generated with rgamma(100,3,0.01))
#' data_file <- system.file("extdata", "small.txt", package = "marp", mustWork = TRUE)
#' data <- read.table(data_file)$V1
#'
#' # set some parameters
#' m = 10  # number of iterations for MLE optimization
#' t = seq(100, 200, by=10)  # time intervals
#' y = 304  # cut-off year for estimating probablity
#' model_gen = 2  # underlying true model
#'
#' # model selection and averaging
#' result <- marp::marp(data, t, m, y, which.model = model_gen)
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
                USE.NAMES = T, SIMPLIFY = F)
  ## estimated mean, (logit) probability and (log) hazard rates from six renewal models
  mu_hat <- out$mu_hat
  pr_hat <- out$pr_hat
  haz_hat <- matrix(out$haz_hat, length(t), 6)
  ## find the best model (minimum AIC value)
  which_aic <- which.min(out$AIC)
  ## min. AIC
  min_aic <- out$AIC[which_aic]
  ## delta AIC for each model (diff bewteen AIC of each model and the min. AIC)
  delta_aic <- out$AIC - min_aic
  ## compute AIC weights
  aic_weight <- exp(-.5 * delta_aic) / sum(exp(-.5 * delta_aic))
  # which_bic <- which.min(out$BIC)
  # min_bic <-  out$BIC[which_bic]
  # delta_bic <- out$BIC - min_bic
  # bic_weight <- exp(-.5 * delta_bic) / sum(exp(-.5 * delta_bic))
  ## model-averaged mean, (logit) probability and (log) hazard rates
  ## using AIC weights
  mu_aic <- mu_hat %*% aic_weight
  pr_aic <- pr_hat %*% aic_weight
  haz_aic <- haz_hat %*% aic_weight
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
  return(append(out, out1))
}
