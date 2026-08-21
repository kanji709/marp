#' A function to generate (double) bootstrap samples and fit Gamma renewal model
#' @param n Number of inter-event times generated in each bootstrap sample.
#' @param t Time points at which log-hazards are evaluated.
#' @param B Number of bootstrap samples.
#' @param BB Number of double-bootstrap samples per bootstrap sample.
#' @param m Positive integer controlling repeated random-start optimizations;
#'   see [marp()].
#' @param par_hat Length-12 vector containing `par1` for models 1--6 followed
#'   by `par2` for models 1--6; Gamma uses shape and rate.
#' @param mu_hat Length-6 vector of model-specific mean estimates.
#' @param pr_hat Length-6 vector of model-specific logit event probabilities.
#' @param haz_hat Matrix of model-specific log-hazards with `length(t)` rows
#'   and six model columns.
#' @param y Time point at which logit event probabilities are evaluated.
#' @return A list of bootstrap estimates and variance/T-statistic quantities.
#'   Components beginning with `pr_` use the logit-probability scale and those
#'   beginning with `haz_` use the log-hazard scale.
#' \describe{
#' \item{mu_star}{Estimated mean from bootstrapped samples }
#' \item{pr_star}{Logit event probabilities from bootstrap samples}
#' \item{haz_star}{Log-hazards from bootstrap samples}
#' \item{mu_var_hat}{Variance of estimated mean}
#' \item{pr_var_hat}{Variance of logit event probabilities}
#' \item{haz_var_hat}{Variance of log-hazards}
#' \item{mu_var_double}{Variance of estimated mean of bootstrapped samples (via double-bootstrapping)}
#' \item{pr_var_double}{Variance of estimated probability of bootstrapped samples (via double-bootstrapping)}
#' \item{haz_var_double}{Variance of estimated hazard rates of bootstrapped samples (via double-bootstrapping)}
#' \item{mu_Tstar}{Pivot quantity of the estimated mean}
#' \item{pr_Tstar}{Pivot quantity of the estimated probability }
#' \item{haz_Tstar}{Pivot quantity of the estimated hazard rates }
#' }
#'
#' @examples
#' \dontrun{
#' # set some parameters
#' n <- 30 # sample size
#' t <- seq(100, 200, by = 10) # time intervals
#' B <- 100 # number of bootstraps
#' BB <- 100 # number of double-bootstraps
#' m <- 10 # repeated random-start optimization setting
#' par_hat <- c(
#'   3.41361e-03, 2.76268e+00, 2.60370e+00, 3.30802e+02, 5.48822e+00, 2.92945e+02, NA,
#'   9.43071e-03, 2.47598e+02, 1.80102e+00, 6.50845e-01, 7.18247e-01
#' )
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
#' y <- 304 # cut-off year for estimating probablity
#'
#' # generate bootstrapped samples then fit renewal model
#' res <- marp::gamma_bstrp(n, t, B, BB, m, par_hat, mu_hat, pr_hat, haz_hat, y)
#' }
#'
#' @export


gamma_bstrp <- function(n, t, B, BB, m, par_hat, mu_hat, pr_hat, haz_hat, y) {
  ## bootstraps
  bstrp <- replicate(B, stats::rgamma(n, par_hat[2], par_hat[8]))
  ## fit a Gamma renewal model
  star <- apply(bstrp, 2, function(x) gamma_rp(x, t, m, y))
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
  ## double-bootstraps & fit a Gamma renewal model
  double <-
    sapply(1:B, function(i)
      replicate(BB, stats::rgamma(n, par1_star[i], par2_star[i])), simplify = "array")
  star_double <-
    apply(double, c(2, 3), function(x)
      gamma_rp(x, t, m, y))
  ## estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_double <-
    apply(star_double, c(1, 2), function(x)
      sapply(x, '[[', 6))
  pr_double <-
    apply(star_double, c(1, 2), function(x)
      sapply(x, '[[', 7))
  haz_double <-
    apply(star_double, c(1, 2), function(x)
      sapply(x, '[[', 8))
  ## variance of estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_var_double <- apply(mu_double, 2, stats::var)
  pr_var_double <- apply(pr_double, 2, stats::var)
  haz_var_double <- apply(haz_double, c(1, 3), stats::var)
  ## t-statistics of estimated mean, (logit) probability and (log) hazard rates from double bootstraps
  mu_Tstar <- (mu_star - mu_hat[2]) / sqrt(mu_var_double)
  pr_Tstar <- (pr_star - pr_hat[2]) / sqrt(pr_var_double)
  haz_Tstar <- (haz_star - haz_hat[, 2]) / sqrt(haz_var_double)
  return(
    list(
      'mu_star' = mu_star,
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
      'haz_Tstar' = haz_Tstar
    )
  )
}
