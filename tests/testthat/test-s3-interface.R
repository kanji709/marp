model_fit_names <- c(
  "par1", "par2", "logL", "AIC", "BIC", "mu_hat", "pr_hat", "haz_hat"
)

marp_fit_names <- c(
  "par1", "par2", "logL", "AIC", "BIC", "mu_hat", "pr_hat", "haz_hat",
  "weights_AIC", "model_best", "mu_best", "pr_best", "haz_best", "mu_gen",
  "pr_gen", "haz_gen", "mu_aic", "pr_aic", "haz_aic"
)

test_that("model-specific fits retain their list contract and receive a common class", {
  set.seed(42)
  data <- rgamma(30, 3, 0.01)
  t <- seq(100, 120, by = 10)
  y <- 304

  fits <- suppressWarnings(list(
    poisson_rp(data, t, y),
    gamma_rp(data, t, 3, y),
    loglogis_rp(data, t, 3, y),
    weibull_rp(data, t, 3, y),
    lognorm_rp(data, t, y),
    bpt_rp(data, t, 3, y)
  ))

  for (fit in fits) {
    expect_s3_class(fit, "marp_model_fit")
    expect_named(fit, model_fit_names, ignore.order = FALSE)
    expect_length(fit, 8L)
    expect_identical(fit$par1, unclass(fit)$par1)
    expect_identical(fit$haz_hat, unclass(fit)$haz_hat)
    expect_equal(nobs(fit), length(data))

    fit_printed <- capture.output(printed <- withVisible(print(fit)))
    expect_gt(length(fit_printed), 0L)
    expect_false(printed$visible)
    expect_identical(printed$value, fit)

    fit_summary <- summary(fit)
    expect_s3_class(fit_summary, "summary_marp_model_fit")
    summary_printed <- capture.output(summary_visible <- withVisible(print(fit_summary)))
    expect_gt(length(summary_printed), 0L)
    expect_false(summary_visible$visible)
    expect_identical(summary_visible$value, fit_summary)
  }

  expect_equal(fits[[1]]$par1, 0.0034136086430979953)
  expect_equal(fits[[1]]$AIC, 402.79911764803296)
  expect_equal(fits[[1]]$mu_hat, 292.94512187913182)
  expect_equal(coef(fits[[2]]), c(shape = fits[[2]]$par1, rate = fits[[2]]$par2))
  expect_s3_class(logLik(fits[[2]]), "logLik")
})

test_that("failed model-specific fits print and summarise safely", {
  failed <- marp:::.new_marp_model_fit(
    list(
      par1 = NA_real_, par2 = NA_real_, logL = NA_real_, AIC = Inf, BIC = Inf,
      mu_hat = NA_real_, pr_hat = NA_real_, haz_hat = rep(NA_real_, 3)
    ),
    model = "Brownian passage time",
    parameter_names = c("mean", "alpha"),
    call = quote(bpt_rp(data, t, m, y)),
    nobs = 30,
    t = c(100, 110, 120),
    y = 304,
    status = "failed"
  )

  expect_s3_class(failed, "marp_model_fit")
  expect_named(failed, model_fit_names, ignore.order = FALSE)
  expect_match(paste(capture.output(print(failed)), collapse = "\n"), "failed")
  expect_s3_class(summary(failed), "summary_marp_model_fit")
  expect_output(print(summary(failed)), "failed")
})

test_that("marp fits retain all 19 components and use the aggregate methods", {
  set.seed(42)
  data <- rgamma(30, 3, 0.01)
  t <- seq(100, 120, by = 10)
  fit <- suppressWarnings(marp(data, t, 3, 304, which.model = 2))

  expect_s3_class(fit, "marp_fit")
  expect_named(fit, marp_fit_names, ignore.order = FALSE)
  expect_length(fit, 19L)
  expect_identical(fit$weights_AIC, unclass(fit)$weights_AIC)
  expect_identical(fit$model_best, unclass(fit)$model_best)
  expect_equal(nobs(fit), length(data))

  printed_text <- capture.output(printed <- withVisible(print(fit)))
  expect_false(printed$visible)
  expect_identical(printed$value, fit)
  expect_match(paste(printed_text, collapse = "\n"), "Gamma")
  expect_match(paste(printed_text, collapse = "\n"), "AIC weight")

  fit_summary <- summary(fit)
  expect_s3_class(fit_summary, "summary_marp_fit")
  expect_named(
    fit_summary$model_comparison,
    c("model", "par1", "par2", "logLik", "AIC", "BIC", "AIC_weight", "status")
  )
  expect_equal(fit_summary$model_comparison$AIC_weight, fit$weights_AIC)
  expect_output(print(fit_summary), "Model comparison")

  failed_candidate <- fit
  failed_candidate$AIC[6] <- Inf
  failed_candidate$BIC[6] <- Inf
  failed_candidate$weights_AIC[6] <- 0
  attr(failed_candidate, "fit_status")[6] <- "failed"
  expect_output(print(failed_candidate), "failed")
  expect_equal(summary(failed_candidate)$model_comparison$status[6], "failed")
})

test_that("confint.marp_fit validates inputs and delegates without changing the engine", {
  set.seed(42)
  data <- rgamma(30, 3, 0.01)
  fit <- suppressWarnings(marp(data, c(100, 110), 3, 304, which.model = 2))

  expect_error(confint(fit, B = 2, BB = 2), "data.*must be supplied")
  expect_error(confint(fit, data = data), "Both `B` and `BB`")
  expect_error(confint(fit, data = data, B = 2, BB = 2, parm = "mean"),
               "parm.*not supported")

  called <- new.env(parent = emptyenv())
  local_mocked_bindings(
    marp_confint = function(data, m, t, B, BB, alpha, y, which.model, ...) {
      called$arguments <- list(
        data = data, m = m, t = t, B = B, BB = BB, alpha = alpha,
        y = y, which.model = which.model
      )
      marp:::.new_marp_confint(
        list(
          out = list(mu_aic = 1, pr_aic = 2, haz_aic = 3),
          percent_CI = list(),
          student_CI = list(mu_lower_ma = 0, mu_upper_ma = 2)
        ),
        call = quote(marp_confint()),
        nobs = length(data),
        level = 1 - alpha
      )
    },
    .package = "marp"
  )

  result <- confint(fit, data = data, B = 9, BB = 7, level = 0.9)
  expect_s3_class(result, "marp_confint")
  expect_equal(called$arguments$B, 9)
  expect_equal(called$arguments$BB, 7)
  expect_equal(called$arguments$alpha, 0.1)
  expect_identical(called$arguments$data, data)
  expect_output(print(result), "confidence intervals")
})
