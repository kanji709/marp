#' A function to calculate the log-likelihood of Weibull model
#' @param param parameters of Weibull model
#' @param x input data for Weibull model
#' @return returns the value of negative log-likelihood of the Weibull model
#' @export

weibull_logl <- function(param, x) {
  lambda <- exp(param[1])
  k <- exp(param[2])
  logl <- sum(stats::dweibull(x, k, lambda, log = T))
  return(-logl)
}
