#' A function to calculate the log-likelihood of Log-Logistics model
#' @param param parameters of Log-Logistics model
#' @param x input data for Log-Logistics model
#' @return returns the value of negative log-likelihood of the Log-Logistics model
#' @export


loglogis_logl <- function(param, x) {
  alpha <- exp(param[1])
  beta <- exp(param[2])
  logl <- sum(dllog(x, beta, alpha, log = T))
  return(-logl)
}
