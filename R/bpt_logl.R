#' A function to calculate the log-likelihood of BPT model
#' @param param parameters of BPT model
#' @param x input data for BPT model
#' @return returns the value of negative log-likelihood of the BPT model
#' @export

bpt_logl <- function(param, x) {
  mu <- exp(param[1])
  alpha <- exp(param[2])
  n <- length(x)
  logl <-log(1 - 0) + (n / 2) * (log(mu) - log(2 * pi * alpha)) - sum((x - mu) ^ 2 / (2 * x * mu * alpha)) - 3 * sum(log(x)) / 2
  return(-logl)
}
