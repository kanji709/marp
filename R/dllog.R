#' Density function of Log-Logistics model
#' @param x input data for Log-Logistics model
#' @param shape shape parameter of Log-Logistics model
#' @param scale scale parameter of Log-Logistics model
#' @param log logic function to determine whether log of logistics to be returned
#' @return returns the denisty of the Log-Logistics model
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
#' # print result
#' cat("result_1 = ", c(6.609490942787837e-13, 1.463191586609840e-13, 8.880167025529705e-12, 6.387343862593732e-13, 1.552779497074985e-14), "\n")
#' cat("result_2 = ", c(-28.04509957121864, -29.55298614083788, -25.44720074992009, -28.07928769790387, -31.79614475297261), "\n")
#'
#' @export

dllog <- function (x,shape = 1,scale = 1,log = FALSE) {
  fx <- (shape / scale) * (x / scale) ^ {shape - 1} / (1 + (x / scale) ^ shape) ^ 2
  if (log)
    return(log(fx))
  else
    return(fx)
}
