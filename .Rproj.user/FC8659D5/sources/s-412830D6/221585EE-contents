#' A function to fit model-averaged renewal process
#' @param n number of inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param B number of bootstrap samples
#' @param BB number of double-bootstrap samples
#' @param m the number of iterations in nlm
#' @param par_hat estimated parameters
#' @param mu_hat estimated mean inter-event times
#' @param pr_hat estimated time to event probability
#' @param haz_hat estimated hazard rates
#' @param y user-specified time point (used to compute time-to-event probability)
#' @return returns list of estimates after fitting different renewal models on (double) bootstarp samples
#' @export

### Model-Averaged Renewal Processes -------------------------------
marp_bstrp <- function(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y){
  out <- mapply(rbind,
                poisson_bstrp(n,t,B,BB,par_hat,mu_hat,pr_hat,haz_hat,y), ## 1. Poisson renewal model
                gamma_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y), ## 2. Gamma renewal model
                loglogis_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y), ## 3. Log-Logistics renewal model
                weibull_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y), ## 4. Weibull renewal model
                lognorm_bstrp(n,t,B,BB,par_hat,mu_hat,pr_hat,haz_hat,y), ## 5. Log-Normal renewal model
                bpt_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y) ## 6. BPT renewal model
  )
  return(list("mu_star" = out$mu_star,
              "pr_star" = out$pr_star,
              "haz_star" = array(out$haz_star, c(length(t), 6, B)),
              "mu_var_hat" = out$mu_var_hat,
              "pr_var_hat" = out$pr_var_hat,
              "haz_var_hat" = array(out$haz_var_hat, c(length(t), 6, 1)),
              "mu_var_double" = out$mu_var_double,
              "pr_var_double" = out$pr_var_double,
              "haz_var_double" = array(out$haz_var_double, c(length(t), 6, B)),
              "mu_Tstar" = out$mu_Tstar,
              "pr_Tstar" = out$pr_Tstar,
              "haz_Tstar" = array(out$haz_Tstar, c(length(t), 6, B))))
}
