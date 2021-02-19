#' Propbability function of Log-Logistics model
#' @param q input quantile for Log-Logistics model
#' @param shape shape parameter of Log-Logistics model
#' @param scale scale parameter of Log-Logistics model
#' @param lower.tail logic function to determine whether lower tail probability to be returned
#' @param log.p logic function to determine whether log of logistics to be returned
#' @return returns the probability of the Log-Logistics model
#' @export

pllog <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  Fx <- 1 / (1 + (q / scale) ^ { -shape})
  if (!lower.tail)
    Fx <- 1 - Fx
  if (log.p)
    Fx <- log(Fx)
  return(Fx)
}
