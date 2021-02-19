#' An utility function to calculate upper limit of T statistics
#' @param up lower limit
#' @param hat estimates
#' @param sigmasq variance
#' @param Tstar T statistics estimated from bootstrap samples
#' @param weights model weights
#' @param B number of bootstraps
#' @param alpha significance level
#' @return returns list of percentile bootstrap intervals (including the best and generating model approach).
#' @export

upperT <-  function(up, hat, sigmasq, Tstar, weights, B, alpha) {
  lowerT <- (hat - up) / sqrt(sigmasq)
  temp <- sapply(1:6, function(i) weights[i] * sum(Tstar[i, ] <= lowerT[i]) / B)
  return(sum(temp) - alpha / 2)
}
