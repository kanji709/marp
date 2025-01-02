# 
# 
# x <- 0.03
# 
# lowert = function(x) {
#   ((haz.hat - x) / sqrt(var))
#   
# }
# 
# temp <- NULL
# x <- 0.03
# upperT = function(x) {
#   lowert <- ((haz.hat - x) / sqrt(var))
#   temp <- NULL
#   
#   for (i in 1:6) {
#     temp[i] <- sum(Tstar[, i] <= lowert[i]) / B
#     #sapply(1:6,function(i) weights[i]*sum(Tstar[,i] <= lowert(x)[i])/B)
#   }
#   
#   
#   return(aic.weights %*% temp - 0.025)
# }

upperT = function(x){
  
  temp <- NULL
  
  lowert<-((haz.hat - x)/sqrt(var))
  
  for(i in 1:6){
    temp[i] <- aic.weights[i]*sum(Tstar[,i] <= lowert[i])/B
    #sapply(1:6,function(i) weights[i]*sum(Tstar[,i] <= lowert(x)[i])/B)
    
  }
  
  return(-sum(temp)-0.05/2))
}


uniroot(upperT,lower= 0,upper=100)







