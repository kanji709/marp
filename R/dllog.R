#' Density function of Log-Logistics model
#' @param x Numeric vector of positive quantiles.
#' @param shape Positive shape parameter of the log-logistic distribution.
#' @param scale Positive scale parameter of the log-logistic distribution.
#' @param log Logical; if `TRUE`, return log-densities.
#' @return A numeric vector of densities, or log-densities when `log = TRUE`.
#' @examples
#' x <- as.numeric(c(350., 450., 227., 352., 654.))
#' # set paramters
#' shape <- 5
#' scale <- 3
#' log <- FALSE
#' result_1 <- marp::dllog(x, shape, scale, log)
#'
#' # alternatively, set log == TRUE
#' log <- TRUE
#' result_2 <- marp::dllog(x, shape, scale, log)
#'
#' @export

dllog <- function (x,shape = 1,scale = 1,log = FALSE) {
  fx <- (shape / scale) * (x / scale) ^ {shape - 1} / (1 + (x / scale) ^ shape) ^ 2
  if (log)
    return(log(fx))
  else
    return(fx)
}
