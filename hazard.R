##################### Hazard functions ########################


#Exp Distribution
hazard_exp <- function(x, par) {
  return(rep(par, length(x)))
}


##Log-logistic Distribution
hazard_llog <-
  function(x, par) {
    ((par[2] / par[1]) * (x / par[1]) ^ (par[2] - 1)) / (1 + (x / par[1]) ^
                                                           par[2])
  }
#  dllog(x,par[2],par[1])/(1-pllog(x,par[2],par[1]))


##Gamma Distribution
hazard_gamma <-
  function(x, par) {
    dgamma(x, par[1], par[2]) / (1 - pgamma(x, par[1], par[2]))
  }
#(par[2]*x)^(par[1]-1)*exp(-par[2]*x)/(pgamma(par[2]*x,par[1],lower.tail = F)*gamma(par[1]))}


##Weibull Distribution
hazard_weibull <-
  function(x, par) {
    (par[2] / par[1]) * (x / par[1]) ^ (par[2] - 1)
  }


##Log-normal Distribution
hazard_logn <-
  function(x, par) {
    exp(-(log(x) - par[1]) ^ 2 / (2 * par[2] ^ 2)) / (x * par[2] * sqrt(2 *
                                                                          pi) * (1 - pnorm((log(
                                                                            x
                                                                          ) - par[1]) / par[2])))
  }


##BPT Distribution
hazard_bpt <-
  function(x, par) {
    sqrt(par[1] / (2 * pi * par[2] ^ 2 * x ^ 3)) * exp(-(x - par[1]) ^ 2 / (2 *
                                                                              par[1] * par[2] ^ 2 * x)) / (1 - (pnorm((x - par[1]) / (
                                                                                par[2] * sqrt(par[1] * x)
                                                                              )) + exp(2 / (par[2] ^ 2)) * pnorm(-(x + par[1]) / (
                                                                                par[2] * sqrt(par[1] * x)
                                                                              ))))
  }


##Main Hazard function compute the hazard for all 6 models
hazard <- function(par,x) {
  
  
  #Compute hazard for each model with functions above
  haz <- mapply(
    c,
    hazard_exp(x, par[1]),
    hazard_gamma(x, c(par[2],par[8])),
    hazard_llog(x, c(par[3],par[9])),
    hazard_weibull(x, c(par[4],par[10])),
    hazard_logn(x, c(par[5],par[11])),
    hazard_bpt(x, c(par[6],par[12])),
    USE.NAMES = F,
    SIMPLIFY = T
  )
  
  #Combine the results in matrix
  out <- t(matrix(unlist(haz),6,length(x)))
  
  #Function return hazard rates
  return(out)
  
  
}


