
library(Rcpp)
library(microbenchmark)

cxx_verbose <- FALSE

# make it reproducible
set.seed(12345)

# original function
dllog <- function (x,shape = 1,scale = 1,log = FALSE) {
  fx <- (shape / scale) * (x / scale) ^ {shape - 1} / (1 + (x / scale) ^ shape) ^ 2
  if (log)
    return(log(fx))
  else
    return(fx)
}

# Rcpp function
cppFunction(depends=c("Rcpp"), showOutput=cxx_verbose, '
    Rcpp::NumericVector rcpp_dllog(const Rcpp::NumericVector &x, const double shape, const double scale, const bool dolog) {
        long n = x.size();
        Rcpp::NumericVector fx(n);
        for (long i = 0; i < n; i++) {
            fx[i] = (shape / scale) * pow(x[i] / scale, shape - 1) / pow(1.0 + pow(x[i] / scale, shape), 2);
            if (dolog) {
                fx[i] = log(fx[i]);
            }
        }
        return fx;
    }'
)

# Rcpp Eigen function
cppFunction(depends=c("RcppEigen"), showOutput=cxx_verbose, '
    Eigen::ArrayXd eigen_dllog(Eigen::Map<Eigen::ArrayXd> &x, const double shape, const double scale, const bool dolog) {
        Eigen::ArrayXd fx(x.rows(), x.cols());
        fx = (shape / scale) * (x / scale).pow(shape - 1) / (1.0 + (x / scale).pow(shape)).pow(2);
        if (dolog) {
            fx = fx.log();
        }
        return fx;
    }'
)



# load the data and initial estimate
data <- read.table('../marp/data/small.txt')$V1 # rgamma(100,3,0.01)
tmp_init <- cbind(stats::runif(1, 0, 2 * log(stats::median(data))), stats::runif(1, 0, 2 * log(mean(data)) / log(stats::median(data))))

# run the original function
r1 <- dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), log = T)
print(sum(r1))

# run the Rcpp C++ function
r2 <- rcpp_dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), T)
print(sum(r2))
stopifnot(all.equal(r1, r2))

# run the Eigen C++ function
r3 <- eigen_dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), T)
print(sum(r3))
stopifnot(all.equal(r1, r3))

# microbenchmarks
microbenchmark(
    R=dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), log = T),
    Rcpp=rcpp_dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), T),
    RcppEigen=eigen_dllog(data, exp(tmp_init[2]), exp(tmp_init[1]), T),
    times = 10000
)
