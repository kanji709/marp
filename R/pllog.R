#' Propbability function of Log-Logistics model
#' @param q input quantile for Log-Logistics model
#' @param shape shape parameter of Log-Logistics model
#' @param scale scale parameter of Log-Logistics model
#' @param lower.tail logic function to determine whether lower tail probability to be returned
#' @param log.p logic function to determine whether log of logistics to be returned
#' @return returns the probability of the Log-Logistics model
#' @examples
#' q <- c(1, 2, 3, 4)
#' # set paramters
#' shape <- 5
#' scale <- 3
#' log <- FALSE
#' result_1 <- marp::pllog(q, shape, scale, log)
#'
#' # alternatively, set log == TRUE
#' log <- TRUE
#' result_2 <-  marp::pllog(q, shape, scale, log)
#'
#' # print result
#' cat("result_1 = ", c(0.995901639344262 , 0.883636363636364 , 0.500000000000000 , 0.191791633780584), "\n")
#' cat("result_2 = ", c(0.0040983606557377 , 0.1163636363636363 , 0.5000000000000000 , 0.8082083662194159), "\n")
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
