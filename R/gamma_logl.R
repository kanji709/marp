#' A function to calculate the log-likelihood of Gamma model
#' @param param Length-2 numeric vector containing log(shape) and log(rate).
#' @param x Numeric vector of positive observations.
#' @return returns the value of negative log-likelihood of the Gamma model
#'
#' @examples
#' set.seed(42)
#' data <-  rgamma(30,3,0.01)
#'
#' # set some parameters
#' par_hat <- c(2.7626793657057762, 0.0094307059277139432) # estimated parameters
#' param <- log(par_hat) # input parameters for logl function
#'
#' # calculate log-likelihood
#' result <- marp::gamma_logl(param, data)
#'
#' result
#'
#' @export

gamma_logl <- function(param, x) {
  alpha <- exp(param[1]) # shape
  beta <- exp(param[2]) # rate
  logl <- sum(stats::dgamma(x, alpha, beta, log = TRUE)) # log-likelihood
  return(-logl)
}
