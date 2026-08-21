#' Probability function of Log-Logistics model
#' @param q Numeric vector of positive quantiles.
#' @param shape Positive shape parameter of the log-logistic distribution.
#' @param scale Positive scale parameter of the log-logistic distribution.
#' @param lower.tail Logical; if `TRUE`, return lower-tail probabilities,
#'   otherwise return upper-tail probabilities.
#' @param log.p Logical; if `TRUE`, return probabilities on the log scale.
#' @return A numeric vector of probabilities, or log-probabilities when
#'   `log.p = TRUE`.
#' @examples
#' q <- c(1, 2, 3, 4)
#' # set paramters
#' shape <- 5
#' scale <- 3
#' result_1 <- marp::pllog(q, shape, scale)
#'
#' # alternatively, return log-probabilities
#' result_2 <- marp::pllog(q, shape, scale, log.p = TRUE)
#'
#' @export

pllog <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  Fx <- 1 / (1 + (q / scale) ^ { -shape})
  if (!lower.tail)
    Fx <- 1 - Fx
  if (log.p)
    Fx <- log(Fx)
  return(Fx)
}
