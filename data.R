#--------------Simulate data for p times ------------------------------------#
gen.data <- function(n, mu, sig, p) {
  
  da <- NULL
  
  for (i in 1:p) {
    
    da <- cbind(da, rlnorm(n, mu, sig))
    
  }
  
  
  return(da)
  
}
