#' Density function of Log-Logistics model
#' @param x input data for Log-Logistics model
#' @param shape shape parameter of Log-Logistics model
#' @param scale scale parameter of Log-Logistics model
#' @param log logic function to determine whether log of logistics to be returned
#' @return returns the denisty of the Log-Logistics model
#' @export

dllog <- function (x,shape = 1,scale = 1,log = FALSE) {
  fx <- (shape / scale) * (x / scale) ^ {shape - 1} / (1 + (x / scale) ^ shape) ^ 2
  if (log)
    return(log(fx))
  else
    return(fx)
}
