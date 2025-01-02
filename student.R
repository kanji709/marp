student <- function(par,n,B,t,m,R,haz,weights){
  
  double <- bstrp(par, n, B, t, m, R, haz)
  
  Tstar <- matrix(double$Tstar,B*length(t),6)
  
  var <- matrix(double$var.hat,length(t),6)
      
}
