#' A function to calculate Studentized bootstrap confidence interval
#' @param n number of inter-event times
#' @param t user-specified time intervals (used to compute hazard rate)
#' @param B number of bootstrap samples
#' @param BB number of double-bootstrap samples
#' @param m the number of iterations in nlm
#' @param par_hat estimated parameters
#' @param mu_hat estimated mean inter-event times
#' @param pr_hat estimated time to event probability
#' @param haz_hat estimated hazard rates
#' @param weights model weights
#' @param alpha significance level
#' @param y user-specified time point (used to compute time-to-event probability)
#' @param best.model best model based on information criterion (i.e. AIC)
#' @param which.model user-specified genearting (or true underlying if known) model
#' @return returns list of Studentized bootstrap intervals (including the model-averaged approach).
#' @examples
#' \dontrun{
#' # generate random data
#' set.seed(42)
#' data <- rgamma(30, 3, 0.01)
#'
#' # set some parameters
#' m <- 10 # number of iterations for MLE optimization
#' t <- seq(100,200,by=10) # time intervals
#' y <- 304 # cut-off year for estimating probablity
#' B <- 100 # number of bootstraps
#' BB <- 100 # number of double bootstraps
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
#' weights <- c(0.00000000000000000, 0.20999999999999999, 0.02000000000000000, 0.55000000000000004, 0.00000000000000000, 0.22000000000000000) # model weights
#' alpha <- 0.05 # confidence level
#' y <- 304 # cut-off year for estimating probablity
#' best.model <- 2
#' which.model <- 2 # specify the generating model#'
#'
#' # construct Studentized bootstrap confidence interval
#' marp::student_confint(n,B,t,m,BB,par_hat,mu_hat,pr_hat,haz_hat,weights,alpha,y,best.model,which.model)
#' }
#'
#' @export

### Model-Averaged Renewal Processes -------------------------------
student_confint <- function(n,B,t,m,BB,par_hat,mu_hat,pr_hat,haz_hat,weights,alpha,y,best.model,which.model=1) {
  ## parametric double-bootstraps & fit all six renewal models
  double <- marp_bstrp(n,t,B,BB,m,par_hat,mu_hat,pr_hat,haz_hat,y)
  mu_var_hat <- double$mu_var_hat
  pr_var_hat <- double$pr_var_hat
  haz_var_hat <- double$haz_var_hat
  ## from all six renewal modelss
  ## extract t-statistics of estimated mean, (logit) probability and (log) hazard rates
  mu_Tstar <- double$mu_Tstar
  pr_Tstar <- double$pr_Tstar
  haz_Tstar <- double$haz_Tstar
  ## using the generating model
  ## find studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  mu_lower_gen <- mu_hat[which.model] - stats::quantile(mu_Tstar[which.model,], probs = 1 - alpha / 2) * sqrt(mu_var_hat[which.model])
  mu_upper_gen <- mu_hat[which.model] - stats::quantile(mu_Tstar[which.model,], probs = alpha / 2) * sqrt(mu_var_hat[which.model])
  pr_lower_gen <- pr_hat[which.model] - stats::quantile(pr_Tstar[which.model,], probs = 1 - alpha / 2) * sqrt(pr_var_hat[which.model])
  pr_upper_gen <- pr_hat[which.model] - stats::quantile(pr_Tstar[which.model,], probs = alpha / 2) * sqrt(pr_var_hat[which.model])
  haz_lower_gen <- haz_hat[, which.model] - apply(haz_Tstar[, which.model, ], 1, function(x) stats::quantile(x, prob = 1 - alpha / 2)) * sqrt(haz_var_hat[, which.model, 1])
  haz_upper_gen <- haz_hat[, which.model] - apply(haz_Tstar[, which.model, ], 1, function(x) stats::quantile(x, prob = alpha / 2)) * sqrt(haz_var_hat[, which.model, 1])
  ## using the best model
  ## find studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  mu_lower_best <- mu_hat[best.model] - stats::quantile(mu_Tstar[best.model,], probs = 1 - alpha / 2) * sqrt(mu_var_hat[best.model])
  mu_upper_best <- mu_hat[best.model] - stats::quantile(mu_Tstar[best.model,], probs = alpha / 2) * sqrt(mu_var_hat[best.model])
  pr_lower_best <- pr_hat[best.model] - stats::quantile(pr_Tstar[best.model,], probs = 1 - alpha / 2) * sqrt(pr_var_hat[best.model])
  pr_upper_best <- pr_hat[best.model] - stats::quantile(pr_Tstar[best.model,], probs = alpha / 2) * sqrt(pr_var_hat[best.model])
  haz_lower_best <- haz_hat[, best.model] - apply(haz_Tstar[, best.model, ], 1, function(x) stats::quantile(x, prob = 1 - alpha / 2)) * sqrt(haz_var_hat[, best.model, 1])
  haz_upper_best <- haz_hat[, best.model] - apply(haz_Tstar[, best.model, ], 1, function(x) stats::quantile(x, prob = alpha / 2)) * sqrt(haz_var_hat[, best.model, 1])
  ## using all six renewal models
  ## find model-averaged studentized bootstrap CIs of estimated mean, (logit) probability and (log) hazard rates
  mu_lower_ma <- stats::uniroot(function(low) lowerT(low, mu_hat, mu_var_hat, mu_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) mu_hat[i] - max(mu_Tstar[i,]) * sqrt(mu_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.9) * sqrt(mu_var_hat[i]))))$root
  mu_upper_ma <- stats::uniroot(function(up) upperT(up, mu_hat, mu_var_hat, mu_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) mu_hat[i] - stats::quantile(mu_Tstar[i,], prob = 0.1) * sqrt(mu_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) mu_hat[i] - min(mu_Tstar[i,]) * sqrt(mu_var_hat[i]))))$root
  pr_lower_ma <- stats::uniroot(function(low) lowerT(low, pr_hat, pr_var_hat, pr_Tstar, weights, B, alpha),
                         lower = min(sapply(1:6, function(i) pr_hat[i] - max(pr_Tstar[i,]) * sqrt(pr_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.9) * sqrt(pr_var_hat[i]))))$root
  pr_upper_ma <- stats::uniroot(function(up) upperT(up, pr_hat, pr_var_hat, pr_Tstar, weights, B, alpha),
                         lower =  min(sapply(1:6, function(i) pr_hat[i] - stats::quantile(pr_Tstar[i,], prob = 0.1) * sqrt(pr_var_hat[i]))),
                         upper = max(sapply(1:6, function(i) pr_hat[i] - min(pr_Tstar[i,]) * sqrt(pr_var_hat[i]))))$root
  haz_lower_ma <- sapply(1:length(t), function(i)
    stats::uniroot(function(low) lowerT(low, haz_hat[i, ], haz_var_hat[i, , 1], haz_Tstar[i, ,], weights, B, alpha),
            lower = min(sapply(1:6, function(j) haz_hat[i, j] - max(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1]))),
            upper = max(sapply(1:6, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.9) * sqrt(haz_var_hat[i, j, 1]))))$root)
  haz_upper_ma <- sapply(1:length(t), function(i)
    stats::uniroot(function(up) upperT(up, haz_hat[i, ], haz_var_hat[i, , 1], haz_Tstar[i, ,], weights, B, alpha),
            lower = min(sapply(1:6, function(j) haz_hat[i, j] - stats::quantile(haz_Tstar[i, j,], prob = 0.1) * sqrt(haz_var_hat[i, j, 1]))),
            upper = max(sapply(1:6, function(j) haz_hat[i, j] - min(haz_Tstar[i, j,]) * sqrt(haz_var_hat[i, j, 1]))))$root)
  return(list("mu_lower_gen" = unname(mu_lower_gen),
              "mu_upper_gen" = unname(mu_upper_gen),
              "pr_lower_gen" = unname(pr_lower_gen),
              "pr_upper_gen" = unname(pr_upper_gen),
              "haz_lower_gen" = haz_lower_gen,
              "haz_upper_gen" = haz_upper_gen,
              "mu_lower_best" = unname(mu_lower_best),
              "mu_upper_best" = unname(mu_upper_best),
              "pr_lower_best" = unname(pr_lower_best),
              "pr_upper_best" = unname(pr_upper_best),
              "haz_lower_best" = haz_lower_best,
              "haz_upper_best" = haz_upper_best,
              "mu_lower_ma" = mu_lower_ma,
              "mu_upper_ma" = mu_upper_ma,
              "pr_lower_ma" = pr_lower_ma,
              "pr_upper_ma" = pr_upper_ma,
              "haz_lower_ma" = haz_lower_ma,
              "haz_upper_ma" = haz_upper_ma))
}
