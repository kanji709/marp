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
#' @param which.model user-specified genearting (or true underlying if known) model
#' @return returns list of Studentized bootstrap intervals (including the model-averaged approach).
#' @export

### Model-Averaged Renewal Processes -------------------------------
student_confint <- function(n,B,t,m,BB,par_hat,mu_hat,pr_hat,haz_hat,weights,alpha,y,best.model,which.model=1) {
  ## parametric double-bootstraps & fit all six renewal models
  double <- marp_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y)
  mu_var_hat <- double$mu_var_hat
  pr_var_hat <- double$pr_var_hat
  haz_var_hat <- double$haz_var_hat
  ## from all six renewal modelss
  ## extract t-statistics of estimated mean, (logit) probability and (log) hazard rates
  mu_Tstar <- double$mu_Tstar
  pr_Tstar <- double$pr_Tstar
  haz_Tstar <- double$haz_Tstar
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
  mu_lower_ma <- stats::uniroot(function(low) lowerT(low, mu_hat, mu_var_hat, mu_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) mu_hat[i] - max(mu_Tstar[i,]) * sqrt(mu_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.9) * sqrt(mu_var_hat[i]))))$root
  mu_upper_ma <- stats::uniroot(function(up) upperT(up, mu_hat, mu_var_hat, mu_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.1) * sqrt(mu_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) mu_hat[i] - min(mu_Tstar[i,]) * sqrt(mu_var_hat[i]))))$root
  pr_lower_ma <- stats::uniroot(function(low) lowerT(low, pr_hat, pr_var_hat, pr_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) pr_hat[i] - max(pr_Tstar[i,]) * sqrt(pr_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.9) * sqrt(pr_var_hat[i]))))$root
  pr_upper_ma <- stats::uniroot(function(up) upperT(up, pr_hat, pr_var_hat, pr_Tstar, weights, B, alpha),
                         lower =  min(sapply(1:6, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.1) * sqrt(pr_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) pr_hat[i] - min(pr_Tstar[i,]) * sqrt(pr_var_hat[i]))))$root
  haz_lower_ma <- sapply(1:length(t), function(i)
    stats::uniroot(function(low) lowerT(low, haz_hat[i, ], haz_var_hat[i, , 1], haz_Tstar[i, ,], weights, B, alpha),
            lower = min(sapply(1:6, function(j) haz_hat[i, j] - max(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1]))),
            upper = max(sapply(1:6, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.9) * sqrt(haz_var_hat[i, j, 1]))))$root)
  haz_upper_ma <- sapply(1:length(t), function(i)
    stats::uniroot(function(up) upperT(up, haz_hat[i, ], haz_var_hat[i, , 1], haz_Tstar[i, ,], weights, B, alpha),
            lower = min(sapply(1:6, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.1) * sqrt(haz_var_hat[i, j, 1]))),
            upper = max(sapply(1:6, function(j) haz_hat[i, j] - min(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1]))))$root)
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
