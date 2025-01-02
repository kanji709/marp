#--------------------------------------- Parametric and Double Bootstrapping -------------------------------------------------------------------#
#----------------------------------------- Exp Distribution--------------------------------------------#
bstrp.exp <- function(par, n, B, t, R, haz) {
  ### Exp Distribution
  bstrp <- replicate(B, rexp(n, par[1]))
  
  intro <- sapply(1:B, function(i)
    prelim.exp(bstrp[, i], t))
  
  par.star <- unlist(intro[1, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rexp(n, par.star[i])), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.exp(x, t)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[1]) / sqrt(var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}
###  Log-Logistic Distribution
bstrp.llog <- function(par, n, B, t, m, R, haz) {
  bstrp <- replicate(B, rllog(n, par[3], par[9]))
  
  intro <- sapply(1:B, function(i)
    prelim.llog(bstrp[, i], t, m))
  
  par1.star <- unlist(intro[1, ])
  
  par2.star <- unlist(intro[2, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rllog(n, par1.star[i], par2.star[i])), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.llog(x, t, m)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[3]) / sqrt(var.double)
  
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}


### Gamma Distribution
bstrp.gamma <- function(par, n, B, t, m, R, haz) {
  bstrp <- replicate(B, rgamma(n, par[2], par[8]))
  
  intro <- sapply(1:B, function(i)
    prelim.gamma(bstrp[, i], t, m))
  
  par1.star <- unlist(intro[1, ])
  par2.star <- unlist(intro[2, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rgamma(n, par1.star[i], par2.star[i])), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.gamma(x, t, m)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[2]) / sqrt(var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}


### Weibull Distribution
bstrp.weibull <- function(par, n, B, t, m, R, haz) {
  bstrp <- replicate(B, rweibull(n, par[10], par[4]))
  
  intro <- sapply(1:B, function(i)
    prelim.weibull(bstrp[, i], t, m))
  
  par1.star <- unlist(intro[2, ])
  
  par2.star <- unlist(intro[1, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rweibull(n, par1.star[i], par2.star[i])), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.weibull(x, t, m)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[4]) / sqrt(var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}


### Log-Normal Distribution
bstrp.lnorm <- function(par, n, B, t, R, haz) {
  bstrp <- replicate(B, rlnorm(n, par[5], par[11]))
  
  intro <- sapply(1:B, function(i)
    prelim.lnorm(bstrp[, i], t))
  
  par1.star <- unlist(intro[1, ])
  
  par2.star <- unlist(intro[2, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rlnorm(n, par1.star[i], par2.star[i])), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.lnorm(x, t)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[5]) / sqrt(var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}

### BPT Distribution
bstrp.bpt <- function(par, n, B, t, m, R, haz) {
  bstrp <- replicate(B, rinvgauss(n, par[6], par[6] / par[12] ^ 2))
  
  intro <- sapply(1:B, function(i)
    prelim.bpt(bstrp[, i], t, m))
  
  par1.star <- unlist(intro[1, ])
  
  par2.star <- unlist(intro[2, ])
  
  haz.star <- matrix(unlist(intro[6, ]), length(t), B)
  
  var <- apply(haz.star, 1, var)
  
  double <-
    sapply(1:B, function(i)
      replicate(R, rinvgauss(n, par1.star[i], par1.star / par2.star ^ 2)), simplify = F)
  
  intro.double <-
    sapply(1:B, function(i)
      apply(double[[i]], 2, function(x)
        prelim.bpt(x, t, m)), simplify = F)
  
  haz.double <-
    sapply(1:B, function(j)
      sapply(1:R, function(i)
        intro.double[[j]][[i]]$hazard), simplify = F)
  
  var.double <-
    sapply(1:B, function(i)
      apply(matrix(unlist(haz.double[[i]]), length(t), R), 1, function(x)
        var(x)))
  
  Tstar <- (haz.star - haz[6]) / sqrt(var.double)
  
  return(
    list(
      'haz.star' = haz.star,
      'var.hat' = var,
      'var.double' = var.double,
      'Tstar' = Tstar
    )
  )
  
}

bstrp <-
  function(par, n, B, t, m, R, haz) {
    mapply(
      c,
      bstrp.exp(par, n, B, t, R, haz),
      bstrp.llog(par, n, B, t, m, R, haz),
      bstrp.gamma(par, n, B, t, m, R, haz),
      bstrp.weibull(par, n, B, t, m, R, haz),
      bstrp.lnorm(par, n, B, t, R, haz),
      bstrp.bpt(par, n, B, t, m, R, haz)
    )
  }
