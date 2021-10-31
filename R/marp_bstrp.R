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
#'
#' @examples
#' \dontrun{
#' # set some parameters
#' n <- 30 # sample size
#' t <- seq(100, 200, by = 10) # time intervals
#' B <- 100 # number of bootstraps
#' BB <- 100 # number of double-bootstraps
#' m <- 10 # number of iterations for MLE optimization
#' par_hat <- c(3.4136086430979953e-03, 2.7626793657057762e+00, 2.6037039674870583e+00, 3.3080162440951688e+02, 5.4882183788378658e+00, 2.9294512422957860e+02, NA, 9.4307059277139432e-03, 2.4759796859031687e+02, 1.8010183507666513e+00, 6.5084541680686814e-01, 7.1824719073918109e-01)
#' mu_hat <- c(292.94512187913182, 292.94512912200048, 319.72017228620746, 294.16945213908519, 298.87285747700128, 292.94512422957860)
#' pr_hat <- c(0.60038574701819891, 0.42154974433034809, 0.53433568234281148, 0.30779792692414687, 0.56416103510057725, 0.61794524610544410)
#' haz_hat <-   matrix(c(
#'   -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -5.6799852941338829,
#'   -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -5.6799852941338829, -6.0942031084732298,
#'   -5.9967873794574516, -5.9117418563554684, -5.8368230853439300, -5.7703089176306639, -5.7108525626839901, -5.6573839062669986,
#'   -5.6090408956082456, -5.5651206740587922, -5.5250440506799734, -5.4883291920475745, -6.0990192429336094, -5.9701664705134210,
#'   -5.8576899644670348, -5.7593884711134971, -5.6734972529860741, -5.5985621349393231, -5.5333565788683616, -5.4768259914915305,
#'   -5.4280496904694857, -5.3862145095364315, -5.3505961502861927, -6.1714638710963881, -6.0951186680582552, -6.0254209583640863,
#'   -5.9613052806725335, -5.9019434350392981, -5.8466788789061646, -5.7949823391436279, -5.7464209045603756, -5.7006359661738628,
#'   -5.6573271297614109, -5.6162402596857071, -5.9235521978533958, -5.8023896004395645, -5.7047473880293342, -5.6252373537796752,
#'   -5.5599409055534252, -5.5059486025117375, -5.4610610586440487, -5.4235891601883868, -5.3922173604047572, -5.3659081375131672,
#'   -5.3438339586221275, -5.7911126719889303, -5.6765973314326752, -5.5892417143301261, -5.5216608261560411, -5.4687921205249133,
#'   -5.4270729562323066, -5.3939387902533049, -5.3675067327627373, -5.3463701567645607, -5.3294619641245422, -5.3159614865560094
#' ),length(t),6)
#' y <- 304 # cut-off year for estimating probablity
#'
#' # fix the random seed
#' set.seed(42)
#'
#' # generate bootstrapped samples then fit renewal model
#' res <- marp::marp_bstrp(n, t, B, BB, m, par_hat, mu_hat, pr_hat, haz_hat, y)
#' }
#'
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
