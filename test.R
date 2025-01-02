require(FAdist)
require(statmod)

n <- 50 
mu <- 3
sig <- 1

m = 10
B = 10
R = 5

P = 1
t = 1
set <- rlnorm(n,mu,sig)

basic <- prelim(set,t,m)

par <- c(basic$par1,basic$par2)

haz.hat <- c(basic$hazard)

double <- bstrp(par,n,B,t,m,R,haz = haz.hat)

var <- c(double$var.hat)

Tstar <- matrix(double$Tstar,B,6)

aic.weights <- basic$weights.AIC
