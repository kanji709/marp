#' A function to calculate Studentized bootstrap confidence interval
#' @param n number of inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param B number of bootstrap samples
#' @param BB number of double-bootstrap samples
#' @param m the number of iterations in nlm
#' @param par_hat estimated parameters
#' @param mu_hat estimated mean inter-event times
#' @param pr_hat estimated time to event probability
#' @param haz_hat estimated hazard rates
#' @param weights model weights
#' @param alpha significance level
#' @param y user-specified time point (used to compute time-to-event probability)
#' @param best.model best model based on information criterion (i.e. AIC)
#' @param which.model user-specified generating (or true underlying if known) model
#'
#' @return returns list of Studentized bootstrap intervals (including the model-averaged approach).
#' \describe{
#' \item{mu_lower_gen}{Lower limit of the studentized bootstrap confidence interval of the estimated mean based on the generating model}
#' \item{mu_upper_gen}{Upper limit of the studentized bootstrap confidence interval of the estimated mean based on the generating model}
#' \item{mu_lower_best}{Lower limit of the studentized bootstrap confidence interval of the estimated mean based on the best model}
#' \item{mu_upper_best}{Upper limit of the studentized bootstrap confidence interval of the estimated mean based on the best model}
#' \item{pr_lower_gen}{Lower limit of the studentized bootstrap confidence interval of the estimated probabilities  based on the generating model}
#' \item{pr_upper_gen}{Upper limit of the studentized bootstrap confidence interval of the estimated probabilities  based on the generating model}
#' \item{pr_lower_best}{Lower limit of the studentized bootstrap confidence interval of the estimated probabilities  based on the best model}
#' \item{pr_upper_best}{Upper limit of the studentized bootstrap confidence interval of the estimated probabilities  based on the best model}
#' \item{haz_lower_gen}{Lower limit of the studentized bootstrap confidence interval of the estimated hazard rates  based on the generating model}
#' \item{haz_upper_gen}{Upper limit of the studentized bootstrap confidence interval of the estimated hazard rates  based on the generating model}
#' \item{haz_lower_best}{Lower limit of the studentized bootstrap confidence interval of the estimated hazard rates  based on the best model}
#' \item{haz_upper_best}{Upper limit of the studentized bootstrap confidence interval of the estimated hazard rates  based on the best model}
#' \item{mu_lower_ma}{Lower limit of model-averaged studentized bootstrap confidence interval of the estimated mean }
#' \item{mu_upper_ma}{Upper limit of model-averaged studentized bootstrap confidence interval of the estimated mean }
#' \item{pr_lower_ma}{Lower limit of model-averaged studentized bootstrap confidence interval of the estimated probabilities  }
#' \item{pr_upper_ma}{Upper limit of model-averaged studentized bootstrap confidence interval of the estimated probabilities  }
#' \item{haz_lower_ma}{Lower limit of model-averaged studentized bootstrap confidence interval of the estimated hazard rates  }
#' \item{haz_upper_ma}{Upper limit of model-averaged studentized bootstrap confidence interval of the estimated hazard rates  }
#' }
#'
#' @examples
#' \dontrun{
#' # generate random data
#' set.seed(42)
#' data <- rgamma(30, 3, 0.01)
#'
#' # set some parameters
#' n <- 30 # sample size
#' m <- 10 # number of iterations for MLE optimization
#' t <- seq(100,200,by=10) # time intervals
#' y <- 304 # cut-off year for estimating probablity
#' B <- 100 # number of bootstraps
#' BB <- 100 # number of double bootstraps
#' par_hat <- c(
#'   3.41361e-03, 2.76268e+00, 2.60370e+00, 3.30802e+02, 5.48822e+00, 2.92945e+02, NA,
#'   9.43071e-03, 2.47598e+02, 1.80102e+00, 6.50845e-01, 7.18247e-01)
#' mu_hat <- c(292.94512, 292.94513, 319.72017, 294.16945, 298.87286, 292.94512)
#' pr_hat <- c(0.60039, 0.42155, 0.53434, 0.30780, 0.56416, 0.61795)
#' haz_hat <-   matrix(c(
#'   -5.67999, -5.67999, -5.67999, -5.67999, -5.67999, -5.67999,
#'   -5.67999, -5.67999, -5.67999, -5.67999, -5.67999, -6.09420,
#'   -5.99679, -5.91174, -5.83682, -5.77031, -5.71085, -5.65738,
#'   -5.60904, -5.56512, -5.52504, -5.48833, -6.09902, -5.97017,
#'   -5.85769, -5.75939, -5.67350, -5.59856, -5.53336, -5.47683,
#'   -5.42805, -5.38621, -5.35060, -6.17146, -6.09512, -6.02542,
#'   -5.96131, -5.90194, -5.84668, -5.79498, -5.74642, -5.70064,
#'   -5.65733, -5.61624, -5.92355, -5.80239, -5.70475, -5.62524,
#'   -5.55994, -5.50595, -5.46106, -5.42359, -5.39222, -5.36591,
#'   -5.34383, -5.79111, -5.67660, -5.58924, -5.52166, -5.46879,
#'   -5.42707, -5.39394, -5.36751, -5.34637, -5.32946, -5.31596
#' ),length(t),6)
#' weights <- c(0.00000, 0.21000, 0.02000, 0.55000, 0.00000, 0.22000) # model weights
#' alpha <- 0.05 # confidence level
#' y <- 304 # cut-off year for estimating probablity
#' best.model <- 2
#' which.model <- 2 # specify the generating model#'
#'
#' # construct Studentized bootstrap confidence interval
#' marp::student_confint(
#'   n,B,t,m,BB,par_hat,mu_hat,pr_hat,haz_hat,weights,alpha,y,best.model,which.model
#' )
#' }
#'
#' @export

### Model-Averaged Renewal Processes -------------------------------
student_confint <- function(n,B,t,m,BB,par_hat,mu_hat,pr_hat,haz_hat,weights,alpha,y,best.model,which.model=1) {
  ## parametric double-bootstraps & fit all six renewal models
  double <- marp_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y)
  mu_var_hat <- double$mu_var_hat
  pr_var_hat <- double$pr_var_hat
  haz_var_hat <- double$haz_var_hat
  ## from all six renewal models
  ## extract t-statistics of estimated mean, (logit) probability and (log) hazard rates
  mu_Tstar <- double$mu_Tstar
  pr_Tstar <- double$pr_Tstar
  haz_Tstar <- double$haz_Tstar

  active_models <- which(weights > 0 & is.finite(weights))
  if (length(active_models) == 0L) {
    stop("No positive-weight candidate models are available for Studentized CIs.", call. = FALSE)
  }

  # --- Remove invalid bootstrap replicates for active models only. Failed
  # zero-weight candidates are ignored in model-averaged interval construction.

  # For the estimated mean (mu) bootstrap t-statistics:
  # Columns = bootstrap replicates, rows = models.
  valid_mu <- apply(mu_Tstar[active_models, , drop = FALSE], 2, function(x) all(is.finite(x)))
  removed_mu <- sum(!valid_mu)
  total_B <- ncol(mu_Tstar)
  percent_mu <- removed_mu / total_B * 100
  cat(sprintf("Removed %d mu replicates (%.1f%%)\n", removed_mu, percent_mu))
  mu_Tstar <- mu_Tstar[, valid_mu, drop = FALSE]

  # For the estimated time-to-event probability (pr):
  valid_pr <- apply(pr_Tstar[active_models, , drop = FALSE], 2, function(x) all(is.finite(x)))
  removed_pr <- sum(!valid_pr)
  percent_pr <- removed_pr / total_B * 100
  cat(sprintf("Removed %d pr replicates (%.1f%%)\n", removed_pr, percent_pr))
  pr_Tstar <- pr_Tstar[, valid_pr, drop = FALSE]

  # For the hazard rate (haz):
  # 3rd dimension = bootstrap replicates.
  valid_haz <- apply(haz_Tstar[, active_models, , drop = FALSE], 3, function(x) all(is.finite(x)))
  removed_haz <- sum(!valid_haz)
  total_B_haz <- dim(haz_Tstar)[3]
  percent_haz <- removed_haz / total_B_haz * 100
  cat(sprintf("Removed %d hazard replicates (%.1f%%)\n", removed_haz, percent_haz))
  haz_Tstar <- haz_Tstar[, , valid_haz, drop = FALSE]

  ## using the generating model
  ## find studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  mu_lower_gen <- mu_hat[which.model] - stats::quantile(mu_Tstar[which.model,], probs = 1 - alpha / 2) * sqrt(mu_var_hat[which.model])
  mu_upper_gen <- mu_hat[which.model] - stats::quantile(mu_Tstar[which.model,], probs = alpha / 2) * sqrt(mu_var_hat[which.model])
  pr_lower_gen <- pr_hat[which.model] - stats::quantile(pr_Tstar[which.model,], probs = 1 - alpha / 2) * sqrt(pr_var_hat[which.model])
  pr_upper_gen <- pr_hat[which.model] - stats::quantile(pr_Tstar[which.model,], probs = alpha / 2) * sqrt(pr_var_hat[which.model])
  haz_lower_gen <- haz_hat[, which.model] - apply(haz_Tstar[, which.model, ], 1, function(x) stats::quantile(x, prob = 1 - alpha / 2)) * sqrt(haz_var_hat[, which.model, 1])
  haz_upper_gen <- haz_hat[, which.model] - apply(haz_Tstar[, which.model, ], 1, function(x) stats::quantile(x, prob = alpha / 2)) * sqrt(haz_var_hat[, which.model, 1])
  ## using the best model
  ## find studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  mu_lower_best <- mu_hat[best.model] - stats::quantile(mu_Tstar[best.model,], probs = 1 - alpha / 2) * sqrt(mu_var_hat[best.model])
  mu_upper_best <- mu_hat[best.model] - stats::quantile(mu_Tstar[best.model,], probs = alpha / 2) * sqrt(mu_var_hat[best.model])
  pr_lower_best <- pr_hat[best.model] - stats::quantile(pr_Tstar[best.model,], probs = 1 - alpha / 2) * sqrt(pr_var_hat[best.model])
  pr_upper_best <- pr_hat[best.model] - stats::quantile(pr_Tstar[best.model,], probs = alpha / 2) * sqrt(pr_var_hat[best.model])
  haz_lower_best <- haz_hat[, best.model] - apply(haz_Tstar[, best.model, ], 1, function(x) stats::quantile(x, prob = 1 - alpha / 2)) * sqrt(haz_var_hat[, best.model, 1])
  haz_upper_best <- haz_hat[, best.model] - apply(haz_Tstar[, best.model, ], 1, function(x) stats::quantile(x, prob = alpha / 2)) * sqrt(haz_var_hat[, best.model, 1])
  ## using all six renewal models
  ## find model-averaged studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  safe_uniroot <- function(expr, label) {
    tryCatch(
      expr,
      error = function(e) {
        warning(
          "student_confint: returning NA for ", label, ": ",
          conditionMessage(e),
          call. = FALSE
        )
        NA_real_
      }
    )
  }

  mu_active <- active_models[is.finite(mu_hat[active_models]) & is.finite(mu_var_hat[active_models])]
  pr_active <- active_models[is.finite(pr_hat[active_models]) & is.finite(pr_var_hat[active_models])]

  if (length(mu_active) > 0L && ncol(mu_Tstar) > 0L) {
    mu_lower_ma <- safe_uniroot(
      stats::uniroot(
        function(low) lowerT(
          low,
          mu_hat[mu_active],
          mu_var_hat[mu_active],
          mu_Tstar[mu_active, , drop = FALSE],
          weights[mu_active],
          B,
          alpha
        ),
        lower = min(sapply(mu_active, function(i) mu_hat[i] - max(mu_Tstar[i,]) * sqrt(mu_var_hat[i]))),
        upper = max(sapply(mu_active, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.9) * sqrt(mu_var_hat[i])))
      )$root,
      "mu_lower_ma"
    )
    mu_upper_ma <- safe_uniroot(
      stats::uniroot(
        function(up) upperT(
          up,
          mu_hat[mu_active],
          mu_var_hat[mu_active],
          mu_Tstar[mu_active, , drop = FALSE],
          weights[mu_active],
          B,
          alpha
        ),
        lower = min(sapply(mu_active, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.1) * sqrt(mu_var_hat[i]))),
        upper = max(sapply(mu_active, function(i) mu_hat[i] - min(mu_Tstar[i,]) * sqrt(mu_var_hat[i])))
      )$root,
      "mu_upper_ma"
    )
  } else {
    mu_lower_ma <- NA_real_
    mu_upper_ma <- NA_real_
  }

  if (length(pr_active) > 0L && ncol(pr_Tstar) > 0L) {
    pr_lower_ma <- safe_uniroot(
      stats::uniroot(
        function(low) lowerT(
          low,
          pr_hat[pr_active],
          pr_var_hat[pr_active],
          pr_Tstar[pr_active, , drop = FALSE],
          weights[pr_active],
          B,
          alpha
        ),
        lower = min(sapply(pr_active, function(i) pr_hat[i] - max(pr_Tstar[i,]) * sqrt(pr_var_hat[i]))),
        upper = max(sapply(pr_active, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.9) * sqrt(pr_var_hat[i])))
      )$root,
      "pr_lower_ma"
    )
    pr_upper_ma <- safe_uniroot(
      stats::uniroot(
        function(up) upperT(
          up,
          pr_hat[pr_active],
          pr_var_hat[pr_active],
          pr_Tstar[pr_active, , drop = FALSE],
          weights[pr_active],
          B,
          alpha
        ),
        lower = min(sapply(pr_active, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.1) * sqrt(pr_var_hat[i]))),
        upper = max(sapply(pr_active, function(i) pr_hat[i] - min(pr_Tstar[i,]) * sqrt(pr_var_hat[i])))
      )$root,
      "pr_upper_ma"
    )
  } else {
    pr_lower_ma <- NA_real_
    pr_upper_ma <- NA_real_
  }

  haz_lower_ma <- sapply(seq_along(t), function(i) {
    haz_active <- active_models[is.finite(haz_hat[i, active_models]) & is.finite(haz_var_hat[i, active_models, 1])]
    if (length(haz_active) == 0L || dim(haz_Tstar)[3] == 0L) {
      return(NA_real_)
    }
    safe_uniroot(
      stats::uniroot(
        function(low) lowerT(
          low,
          haz_hat[i, haz_active],
          haz_var_hat[i, haz_active, 1],
          matrix(haz_Tstar[i, haz_active, ], nrow = length(haz_active)),
          weights[haz_active],
          B,
          alpha
        ),
        lower = min(sapply(haz_active, function(j) haz_hat[i, j] - max(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1]))),
        upper = max(sapply(haz_active, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.9) * sqrt(haz_var_hat[i, j, 1])))
      )$root,
      paste0("haz_lower_ma[", i, "]")
    )
  })
  haz_upper_ma <- sapply(seq_along(t), function(i) {
    haz_active <- active_models[is.finite(haz_hat[i, active_models]) & is.finite(haz_var_hat[i, active_models, 1])]
    if (length(haz_active) == 0L || dim(haz_Tstar)[3] == 0L) {
      return(NA_real_)
    }
    safe_uniroot(
      stats::uniroot(
        function(up) upperT(
          up,
          haz_hat[i, haz_active],
          haz_var_hat[i, haz_active, 1],
          matrix(haz_Tstar[i, haz_active, ], nrow = length(haz_active)),
          weights[haz_active],
          B,
          alpha
        ),
        lower = min(sapply(haz_active, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.1) * sqrt(haz_var_hat[i, j, 1]))),
        upper = max(sapply(haz_active, function(j) haz_hat[i, j] - min(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1])))
      )$root,
      paste0("haz_upper_ma[", i, "]")
    )
  })
  return(list("mu_lower_gen" = unname(mu_lower_gen),
              "mu_upper_gen" = unname(mu_upper_gen),
              "pr_lower_gen" = unname(pr_lower_gen),
              "pr_upper_gen" = unname(pr_upper_gen),
              "haz_lower_gen" = haz_lower_gen,
              "haz_upper_gen" = haz_upper_gen,
              "mu_lower_best" = unname(mu_lower_best),
              "mu_upper_best" = unname(mu_upper_best),
              "pr_lower_best" = unname(pr_lower_best),
              "pr_upper_best" = unname(pr_upper_best),
              "haz_lower_best" = haz_lower_best,
              "haz_upper_best" = haz_upper_best,
              "mu_lower_ma" = mu_lower_ma,
              "mu_upper_ma" = mu_upper_ma,
              "pr_lower_ma" = pr_lower_ma,
              "pr_upper_ma" = pr_upper_ma,
              "haz_lower_ma" = haz_lower_ma,
              "haz_upper_ma" = haz_upper_ma))
}
