dllog <- function (x, shape = 1, scale = 1, log = FALSE) {
  
  fx <- (shape/scale)*(x/scale)^{shape - 1}/(1 + (x/scale)^shape)^2
  
  if (log) 
    
    return(log(fx))
  
  else return(fx)

}

pllog <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE){
  
  Fx <- 1/(1+(q/scale)^{-shape})
  
  if (!lower.tail)
    
    Fx <- 1 - Fx
  
  if (log.p) 
    
    Fx <- log(Fx)
  
  return(Fx)

}
