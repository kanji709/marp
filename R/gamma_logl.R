#' A function to calculate the log-likelihood of Gamma model
#' @param param parameters of Gamma model
#' @param x input data for Gamma model
#' @return returns the value of negative log-likelihood of the Gamma model
#' @export

gamma_logl <- function(param, x) {
  alpha <- exp(param[1]) # shape
  beta <- exp(param[2]) # rate
  logl <- sum(stats::dgamma(x, alpha, beta, log = TRUE)) # log-likelihood
  return(-logl)
}
