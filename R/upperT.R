#' Evaluate the upper-limit equation for a model-averaged T statistic
#' @param up Candidate upper confidence limit.
#' @param hat Vector of model-specific estimates.
#' @param sigmasq Vector of model-specific variance estimates.
#' @param Tstar Matrix of bootstrap T statistics, with models in rows.
#' @param weights Vector of model weights.
#' @param B Number of bootstrap samples represented in `Tstar`.
#' @param alpha Significance level.
#' @return A numeric scalar giving the upper-limit root equation value.
#' @examples
#' # set some parameters
#' up <- 100 # upper bound
#' hat <- rep(150, 6) # estimates obtained from each model
#' sigmasq <- 10 # variance
#' Tstar <- matrix(rep(100,600),6,100) # T statistics estimated from bootstrap samples
#' weights <- rep(1/6, 6) # model weights
#' B <- 100 # number of bootstrapped samples
#' alpha <- 0.05 # confidence level
#'
#' # calculate the upper limit of T statistics
#' res <-  marp::upperT(up, hat, sigmasq, Tstar, weights, B, alpha)
#'
#' res
#'
#' @export

upperT <-  function(up, hat, sigmasq, Tstar, weights, B, alpha) {
  upperT <- (hat - up) / sqrt(sigmasq)
  temp <- sapply(seq_along(weights), function(i) weights[i] * sum(Tstar[i, ] <= upperT[i]) / B)
  return(sum(temp) - alpha / 2)
}
