#' Evaluate the lower-limit equation for a model-averaged T statistic
#' @param low Candidate lower confidence limit.
#' @param hat Vector of model-specific estimates.
#' @param sigmasq Vector of model-specific variance estimates.
#' @param Tstar Matrix of bootstrap T statistics, with models in rows.
#' @param weights Vector of model weights.
#' @param B Number of bootstrap samples represented in `Tstar`.
#' @param alpha Significance level.
#' @return A numeric scalar giving the lower-limit root equation value.
#' @examples
#' # set some parameters
#' low <- 100 # lower bound
#' hat <- rep(150, 6) # estimates obtained from each model
#' sigmasq <- 10 # variance
#' Tstar <- matrix(rep(100,600),6,100) # T statistics estimated from bootstrap samples
#' weights <- rep(1/6, 6) # model weights
#' B <- 100 # number of bootstrapped samples
#' alpha <- 0.05 # confidence level
#'
#' # calculate the upper limit of T statistics
#' res <- marp::lowerT(low, hat, sigmasq, Tstar, weights, B, alpha)
#'
#' res
#'
#' @export

lowerT <-  function(low, hat, sigmasq, Tstar, weights, B, alpha) {
  upperT <- (hat - low) / sqrt(sigmasq)
  temp <- sapply(seq_along(weights), function(i) weights[i] * sum(Tstar[i, ] >= upperT[i]) / B)
  return(sum(temp) - alpha / 2)
}
