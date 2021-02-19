#' A function to generate (double) bootstrap samples and fit Log-Normal renewal model
#' @param n number of inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param B number of bootstrap samples
#' @param BB number of double-bootstrap samples
#' @param par_hat estimated parameters
#' @param mu_hat estimated mean inter-event times
#' @param pr_hat estimated time to event probability
#' @param haz_hat estimated hazard rates
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting Log-Normal renewal model on (double) bootstarp samples
#' @export


lognorm_bstrp <- function(n, t, B, BB, par_hat, mu_hat, pr_hat, haz_hat, y) {
  ## bootstraps
  bstrp <- replicate(B, stats::rlnorm(n, par_hat[5], par_hat[11]))
  ## fit a Gamma Log-Normal model
  star <- apply(bstrp, 2, function(x) lognorm_rp(x, t, y))
  ## parameters
  par1_star <- sapply(star, '[[', 1)
  par2_star <- sapply(star, '[[', 2)
  ## estimated mean, (logit) probability and (log) hazard rates
  mu_star <- sapply(star, '[[', 6)
  pr_star <- sapply(star, '[[', 7)
  haz_star <- sapply(star, '[[', 8)
  ## variance of estimated mean, (logit) probability and (log) hazard rates
  mu_var_hat <- stats::var(mu_star)
  pr_var_hat <- stats::var(pr_star)
  haz_var_hat <- apply(haz_star, 1, stats::var)
  ## double-bootstraps & fit a Log-Normal renewal model
  double <- sapply(1:B, function(i) replicate(BB, stats::rlnorm(n, par1_star[i], par2_star[i])), simplify = "array")
  star_double <- apply(double, c(2, 3), function(x) lognorm_rp(x, t, y))
  ## estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_double <- apply(star_double, c(1, 2), function(x) sapply(x, '[[', 6))
  pr_double <- apply(star_double, c(1, 2), function(x) sapply(x, '[[', 7))
  haz_double <- apply(star_double, c(1, 2), function(x) sapply(x, '[[', 8))
  ## variance of estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_var_double <- apply(mu_double, 2, stats::var)
  pr_var_double <- apply(pr_double, 2, stats::var)
  haz_var_double <- apply(haz_double, c(1, 3), stats::var)
  ## t-statistics of estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_Tstar <- (mu_star - mu_hat[5]) / sqrt(mu_var_double)
  pr_Tstar <- (pr_star - pr_hat[5]) / sqrt(pr_var_double)
  haz_Tstar <- (haz_star - haz_hat[, 5]) / sqrt(haz_var_double)
  return(list('mu_star' = mu_star,
              'pr_star' = pr_star,
              'haz_star' = haz_star,
              'mu_var_hat' =  mu_var_hat,
              'pr_var_hat' =  pr_var_hat,
              'haz_var_hat' = matrix(haz_var_hat, length(t), 1),
              'mu_var_double' = mu_var_double,
              'pr_var_double' = pr_var_double,
              'haz_var_double' = haz_var_double,
              'mu_Tstar' = mu_Tstar,
              'pr_Tstar' = pr_Tstar,
              'haz_Tstar' = haz_Tstar))
}
