upperT <-  function(x,hat,var,Tstar,weights,B,alpha){
  
  lowert <- (hat - x)/sqrt(var)
  
  temp <- sapply(1:6,function(i) weights[i]*sum(Tstar[i,] <= lowert[i])/B)
  
  return(sum(temp)-alpha/2)
  
}

lowerT <-  function(x,hat,var,Tstar,weights,B,alpha){
  
  uppert <- (hat - x)/sqrt(var)
  
  temp <- sapply(1:6,function(i) weights[i]*sum(Tstar[i,] >= uppert[i])/B)
  
  return(sum(temp)-alpha/2)
  
}


student <- function(par, n, B, t, m, R, mean, prob, haz, weights, alpha, mu.true) {
  
  double <- bstrp(par, n, B, t, m, R, mean, prob, haz, mu.true)
  
  mu.var.hat <- double$mu.var.hat
  
  pr.var.hat <- double$pr.var.hat
  
  haz.var.hat <- double$haz.var.hat
  
  mu.Tstar <- double$mu.Tstar
  
  pr.Tstar <- double$pr.Tstar
  
  haz.Tstar <- double$haz.Tstar
  
  mu.lower <- uniroot(
    function(x)
      lowerT(x, mean, mu.var.hat, mu.Tstar, weights, B, alpha),
    lower = -100,
    upper = 100
  )$root
  
  mu.upper <- uniroot(
    function(x)
      upperT(x, mean, mu.var.hat, mu.Tstar, weights, B, alpha),
    lower = 0,
    upper = 100  
  )$root
  
  pr.lower <- uniroot(
    function(x)
      lowerT(x, prob, pr.var.hat, pr.Tstar, weights, B, alpha),
    lower = -10,
    upper = 10
  )$root
  
  pr.upper <- uniroot(
    function(x)
      upperT(x, prob, pr.var.hat, pr.Tstar, weights, B, alpha),
    lower = -10,
    upper = 10
  )$root
  
  
  haz.lower <- sapply(1:length(t), function(i)
    uniroot(
      function(x)
        lowerT(x, haz[i, ], haz.var.hat[i,,1], haz.Tstar[i,,], weights, B, alpha),
      lower = -10,
      upper = 10
    )$root)
  
  haz.upper <- sapply(1:length(t), function(i)
    uniroot(
      function(x)
        upperT(x, haz[i, ], haz.var.hat[i,,1], haz.Tstar[i,,], weights, B, alpha),
      lower = -10,
      upper = 10
    )$root)
  
  return(
    list(
      "mu.lower" = mu.lower,
      "mu.upper" = mu.upper,
      "pr.lower" = pr.lower,
      "pr.upper" = pr.upper,
      "haz.lower" = haz.lower,
      "haz.upper" = haz.upper
    )
  )
  
}


