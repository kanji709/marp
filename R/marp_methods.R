#' Print a model-specific renewal-process fit
#'
#' @param x An object returned by one of the model-specific fitting functions.
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @importFrom stats coef confint logLik nobs
#' @export
print.marp_model_fit <- function(x, ...) {
  model <- attr(x, "model") %||% "Renewal-process model"
  status <- attr(x, "fit_status") %||% "unknown"
  parameter_names <- attr(x, "parameter_names") %||% c("par1", "par2")
  parameters <- stats::setNames(c(x$par1, x$par2), parameter_names)
  parameters <- parameters[!is.na(parameters)]

  cat(model, " fit\n", sep = "")
  if (!is.null(attr(x, "nobs"))) {
    cat("Observations:", attr(x, "nobs"), "\n")
  }
  cat("Status:", status, "\n")
  cat("\nParameter estimates:\n")
  print(parameters)
  cat("\nLog-likelihood:", format(x$logL), "\n")
  cat("AIC:", format(x$AIC), "  BIC:", format(x$BIC), "\n")
  cat("Estimated mean:", format(x$mu_hat), "\n")
  cat("Logit event probability at y =", format(attr(x, "y")), ":",
      format(x$pr_hat), "\n")

  finite_hazard <- sum(is.finite(x$haz_hat))
  cat("Log-hazard evaluations: ", length(x$haz_hat), " stored (",
      finite_hazard, " finite)\n", sep = "")
  invisible(x)
}

#' Summarise a model-specific renewal-process fit
#'
#' @param object An object returned by one of the model-specific fitting
#'   functions.
#' @param ... Additional arguments, currently unused.
#' @return A structured object of class `summary_marp_model_fit`.
#' @export
summary.marp_model_fit <- function(object, ...) {
  parameter_names <- attr(object, "parameter_names") %||% c("par1", "par2")
  parameters <- stats::setNames(c(object$par1, object$par2), parameter_names)
  parameters <- parameters[!is.na(parameters)]
  t <- attr(object, "t")

  out <- list(
    call = attr(object, "call"),
    model = attr(object, "model") %||% "Renewal-process model",
    nobs = attr(object, "nobs"),
    status = attr(object, "fit_status") %||% "unknown",
    parameters = parameters,
    fit_statistics = c(
      logLik = object$logL,
      AIC = object$AIC,
      BIC = object$BIC
    ),
    estimates = list(
      mean = object$mu_hat,
      logit_probability = object$pr_hat,
      y = attr(object, "y"),
      log_hazard = data.frame(
        t = if (is.null(t)) seq_along(object$haz_hat) else t,
        estimate = as.numeric(object$haz_hat)
      )
    )
  )
  class(out) <- c("summary_marp_model_fit", "list")
  out
}

#' Print a summary of a model-specific renewal-process fit
#'
#' @param x A `summary_marp_model_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @export
print.summary_marp_model_fit <- function(x, ...) {
  cat(x$model, " fit summary\n", sep = "")
  if (!is.null(x$nobs)) {
    cat("Observations:", x$nobs, "\n")
  }
  cat("Status:", x$status, "\n")
  cat("\nParameter estimates:\n")
  print(x$parameters)
  cat("\nFit statistics:\n")
  print(x$fit_statistics)
  cat("\nDerived estimates:\n")
  cat("Mean:", format(x$estimates$mean), "\n")
  cat("Logit event probability at y =", format(x$estimates$y), ":",
      format(x$estimates$logit_probability), "\n")
  cat("Log-hazard evaluations:", nrow(x$estimates$log_hazard), "stored\n")
  invisible(x)
}

#' Extract parameters from a model-specific fit
#'
#' @param object A `marp_model_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return A named numeric vector of fitted distribution parameters.
#' @export
coef.marp_model_fit <- function(object, ...) {
  parameter_names <- attr(object, "parameter_names") %||% c("par1", "par2")
  out <- stats::setNames(c(object$par1, object$par2), parameter_names)
  out[!is.na(out)]
}

#' Extract the log-likelihood from a model-specific fit
#'
#' @param object A `marp_model_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return An object of class `logLik`.
#' @export
logLik.marp_model_fit <- function(object, ...) {
  out <- object$logL
  attr(out, "df") <- sum(!is.na(c(object$par1, object$par2)))
  attr(out, "nobs") <- attr(object, "nobs")
  class(out) <- "logLik"
  out
}

#' Number of observations in a model-specific fit
#'
#' @param object A `marp_model_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return The number of observations used for fitting.
#' @export
nobs.marp_model_fit <- function(object, ...) {
  attr(object, "nobs")
}

#' Print a model-averaged renewal-process fit
#'
#' @param x An object returned by [marp()].
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @export
print.marp_fit <- function(x, ...) {
  model_names <- attr(x, "model_names") %||% .marp_model_names
  status <- attr(x, "fit_status") %||% ifelse(is.finite(x$AIC), "ok", "failed")
  comparison <- data.frame(
    Model = model_names,
    AIC = as.numeric(x$AIC),
    `AIC weight` = as.numeric(x$weights_AIC),
    Status = status,
    check.names = FALSE
  )

  cat("Model-averaged renewal-process fit\n")
  if (!is.null(attr(x, "nobs"))) {
    cat("Observations:", attr(x, "nobs"), "\n")
  }
  cat("\nModel comparison:\n")
  print(comparison, row.names = FALSE)
  cat("\nBest model:", model_names[x$model_best], "\n")
  cat("Best-model mean:", format(x$mu_best), "\n")
  cat("Model-averaged mean:", format(x$mu_aic), "\n")
  cat("Model-averaged logit event probability:", format(x$pr_aic), "\n")
  cat("Model-averaged log-hazard evaluations:", length(x$haz_aic), "stored\n")
  invisible(x)
}

#' Summarise a model-averaged renewal-process fit
#'
#' @param object An object returned by [marp()].
#' @param ... Additional arguments, currently unused.
#' @return A structured object of class `summary_marp_fit`.
#' @export
summary.marp_fit <- function(object, ...) {
  model_names <- attr(object, "model_names") %||% .marp_model_names
  status <- attr(object, "fit_status") %||%
    ifelse(is.finite(object$AIC), "ok", "failed")
  which_model <- attr(object, "which.model")
  t <- attr(object, "t")

  out <- list(
    call = attr(object, "call"),
    nobs = attr(object, "nobs"),
    model_comparison = data.frame(
      model = model_names,
      par1 = as.numeric(object$par1),
      par2 = as.numeric(object$par2),
      logLik = as.numeric(object$logL),
      AIC = as.numeric(object$AIC),
      BIC = as.numeric(object$BIC),
      AIC_weight = as.numeric(object$weights_AIC),
      status = status,
      check.names = FALSE
    ),
    selected = list(
      model = model_names[object$model_best],
      model_index = object$model_best,
      mean = object$mu_best,
      logit_probability = object$pr_best,
      log_hazard = data.frame(
        t = if (is.null(t)) seq_along(object$haz_best) else t,
        estimate = as.numeric(object$haz_best)
      )
    ),
    model_averaged = list(
      mean = object$mu_aic,
      logit_probability = object$pr_aic,
      log_hazard = data.frame(
        t = if (is.null(t)) seq_along(object$haz_aic) else t,
        estimate = as.numeric(object$haz_aic)
      )
    ),
    reference = list(
      model = if (is.null(which_model)) NA_character_ else model_names[which_model],
      model_index = which_model,
      mean = object$mu_gen,
      logit_probability = object$pr_gen,
      log_hazard = data.frame(
        t = if (is.null(t)) seq_along(object$haz_gen) else t,
        estimate = as.numeric(object$haz_gen)
      )
    ),
    y = attr(object, "y")
  )
  class(out) <- c("summary_marp_fit", "list")
  out
}

#' Print a summary of a model-averaged renewal-process fit
#'
#' @param x A `summary_marp_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @export
print.summary_marp_fit <- function(x, ...) {
  cat("Model-averaged renewal-process fit summary\n")
  if (!is.null(x$nobs)) {
    cat("Observations:", x$nobs, "\n")
  }
  cat("\nModel comparison:\n")
  print(x$model_comparison, row.names = FALSE)
  cat("\nSelected model:", x$selected$model, "\n")
  cat("Selected-model mean:", format(x$selected$mean), "\n")
  cat("Selected-model logit event probability:",
      format(x$selected$logit_probability), "\n")
  cat("\nModel-averaged estimates:\n")
  cat("Mean:", format(x$model_averaged$mean), "\n")
  cat("Logit event probability at y =", format(x$y), ":",
      format(x$model_averaged$logit_probability), "\n")
  cat("Log-hazard evaluations:",
      nrow(x$model_averaged$log_hazard), "stored\n")
  if (!is.null(x$reference$model_index)) {
    cat("\nReference model:", x$reference$model, "\n")
  }
  invisible(x)
}

#' Number of observations in a model-averaged fit
#'
#' @param object A `marp_fit` object.
#' @param ... Additional arguments, currently unused.
#' @return The number of observations used for fitting.
#' @export
nobs.marp_fit <- function(object, ...) {
  attr(object, "nobs")
}

#' Confidence intervals for a model-averaged renewal-process fit
#'
#' This delegates to [marp_confint()] without changing its bootstrap
#' calculations. The original data are deliberately not stored in `object`, so
#' they must be supplied explicitly.
#'
#' @param object A `marp_fit` object.
#' @param parm Currently unsupported. It must be `NULL`.
#' @param level Confidence level. The default is 0.95.
#' @param data The original inter-event times used to fit `object`.
#' @param B Number of bootstrap samples.
#' @param BB Number of double-bootstrap samples.
#' @param ... Additional arguments passed to [marp_confint()].
#' @return The existing nested confidence-interval result returned by
#'   [marp_confint()], with class `marp_confint`.
#' @export
confint.marp_fit <- function(object, parm = NULL, level = 0.95, data, B, BB, ...) {
  if (!is.null(parm)) {
    stop("`parm` selection is not supported by confint.marp_fit().", call. = FALSE)
  }
  if (missing(data)) {
    stop("`data` must be supplied because marp_fit objects do not store the original data.",
         call. = FALSE)
  }
  if (missing(B) || missing(BB)) {
    stop("Both `B` and `BB` must be supplied for the existing bootstrap procedure.",
         call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  required <- c("t", "m", "y", "which.model")
  metadata <- stats::setNames(lapply(required, function(name) attr(object, name)), required)
  if (any(vapply(metadata, is.null, logical(1)))) {
    stop("The fitted object does not contain the metadata required by marp_confint().",
         call. = FALSE)
  }

  marp_confint(
    data = data,
    m = metadata$m,
    t = metadata$t,
    B = B,
    BB = BB,
    alpha = 1 - level,
    y = metadata$y,
    which.model = metadata$which.model,
    ...
  )
}

#' Print a model-averaged confidence-interval result
#'
#' @param x An object returned by [marp_confint()] or `confint.marp_fit()`.
#' @param ... Additional arguments, currently unused.
#' @return The input object, invisibly.
#' @export
print.marp_confint <- function(x, ...) {
  cat("Model-averaged renewal-process confidence intervals\n")
  if (!is.null(attr(x, "nobs"))) {
    cat("Observations:", attr(x, "nobs"), "\n")
  }
  if (!is.null(attr(x, "level"))) {
    cat("Confidence level:", format(attr(x, "level")), "\n")
  }
  cat("\nModel-averaged point estimates:\n")
  cat("Mean:", format(x$out$mu_aic), "\n")
  cat("Logit event probability:", format(x$out$pr_aic), "\n")
  cat("Log-hazard evaluations:", length(x$out$haz_aic), "stored\n")
  cat("Studentized model-averaged mean interval:",
      format(x$student_CI$mu_lower_ma), "to",
      format(x$student_CI$mu_upper_ma), "\n")
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
