# Internal constructors and metadata helpers for marp S3 objects.

.marp_model_names <- c(
  "Poisson renewal (exponential waiting times)",
  "Gamma",
  "Log-logistic",
  "Weibull",
  "Log-normal",
  "Brownian passage time"
)

.new_marp_model_fit <- function(x, model, parameter_names, call, nobs, t, y,
                                status = "ok") {
  stopifnot(is.list(x))
  attr(x, "model") <- model
  attr(x, "parameter_names") <- parameter_names
  attr(x, "call") <- call
  attr(x, "nobs") <- as.integer(nobs)
  attr(x, "t") <- t
  attr(x, "y") <- y
  attr(x, "fit_status") <- status
  class(x) <- c("marp_model_fit", "list")
  x
}

.new_marp_fit <- function(x, call, nobs, t, y, m, which.model) {
  stopifnot(is.list(x))
  attr(x, "call") <- call
  attr(x, "nobs") <- as.integer(nobs)
  attr(x, "t") <- t
  attr(x, "y") <- y
  attr(x, "m") <- m
  attr(x, "which.model") <- which.model
  attr(x, "model_names") <- .marp_model_names
  attr(x, "fit_status") <- ifelse(is.finite(x$AIC), "ok", "failed")
  class(x) <- c("marp_fit", "list")
  x
}

.new_marp_confint <- function(x, call, nobs, level) {
  stopifnot(is.list(x))
  attr(x, "call") <- call
  attr(x, "nobs") <- as.integer(nobs)
  attr(x, "level") <- level
  class(x) <- c("marp_confint", "list")
  x
}
