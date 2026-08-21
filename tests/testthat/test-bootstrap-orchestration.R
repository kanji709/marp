test_that("marp_bstrp combines existing model bootstrap results", {
  B <- 2
  t <- c(100, 110)

  one_model <- function(...) {
    list(
      mu_star = c(1, 2),
      pr_star = c(3, 4),
      haz_star = matrix(5:8, nrow = length(t)),
      mu_var_hat = 9,
      pr_var_hat = 10,
      haz_var_hat = matrix(c(11, 12), ncol = 1),
      mu_var_double = c(13, 14),
      pr_var_double = c(15, 16),
      haz_var_double = matrix(17:20, nrow = length(t)),
      mu_Tstar = c(21, 22),
      pr_Tstar = c(23, 24),
      haz_Tstar = matrix(25:28, nrow = length(t))
    )
  }

  local_mocked_bindings(
    poisson_bstrp = one_model,
    gamma_bstrp = one_model,
    loglogis_bstrp = one_model,
    weibull_bstrp = one_model,
    lognorm_bstrp = one_model,
    bpt_bstrp = one_model,
    .package = "marp"
  )

  result <- marp_bstrp(
    n = 5, t = t, B = B, BB = 2, m = 3,
    par_hat = rep(1, 12), mu_hat = rep(1, 6),
    pr_hat = rep(1, 6), haz_hat = matrix(1, length(t), 6), y = 120
  )

  expect_named(result, c(
    "mu_star", "pr_star", "haz_star", "mu_var_hat", "pr_var_hat",
    "haz_var_hat", "mu_var_double", "pr_var_double", "haz_var_double",
    "mu_Tstar", "pr_Tstar", "haz_Tstar"
  ))
  expect_equal(dim(result$haz_star), c(length(t), 6, B))
  expect_equal(dim(result$haz_Tstar), c(length(t), 6, B))
  expect_equal(dim(result$mu_Tstar), c(6, B))
})

test_that("marp_confint preserves its nested result contract", {
  set.seed(42)
  dat <- rgamma(20, 3, 0.01)
  t <- c(100, 110)
  base_fit <- suppressWarnings(marp(dat, t, 3, 304, which.model = 2))

  percentile <- list(
    weights_bstp = rep(1 / 6, 6),
    mu_gen = 1, mu_gen_lower = 0, mu_gen_upper = 2,
    mu_best = 1, mu_best_lower = 0, mu_best_upper = 2,
    pr_gen = 1, pr_gen_lower = 0, pr_gen_upper = 2,
    pr_best = 1, pr_best_lower = 0, pr_best_upper = 2,
    haz_gen = rep(1, length(t)), haz_gen_lower = rep(0, length(t)),
    haz_gen_upper = rep(2, length(t)), haz_best = rep(1, length(t)),
    haz_best_lower = rep(0, length(t)), haz_best_upper = rep(2, length(t))
  )
  studentized <- list(mu_lower_ma = 0, mu_upper_ma = 2)

  local_mocked_bindings(
    marp = function(...) base_fit,
    percent_confint = function(...) percentile,
    student_confint = function(...) studentized,
    .package = "marp"
  )

  result <- marp_confint(dat, 3, t, 5, 5, 0.05, 304, 2)

  expect_s3_class(result, "marp_confint")
  expect_named(result, c("out", "percent_CI", "student_CI"))
  expect_identical(result$percent_CI, percentile)
  expect_identical(result$student_CI, studentized)
  expect_named(
    result$out,
    c(
      "par1", "par2", "logL", "AIC", "BIC", "mu_hat", "pr_hat",
      "haz_hat", "weights_AIC", "model_best", "mu_best", "pr_best",
      "haz_best", "mu_gen", "pr_gen", "haz_gen", "mu_aic", "pr_aic",
      "haz_aic", "mu_bstrp", "pr_bstrp", "haz_bstrp"
    ),
    ignore.order = FALSE
  )
})

test_that("student_confint rejects an empty active-model set", {
  t <- c(100, 110)
  mocked_bootstrap <- list(
    mu_var_hat = rep(1, 6),
    pr_var_hat = rep(1, 6),
    haz_var_hat = array(1, c(length(t), 6, 1)),
    mu_Tstar = matrix(0, 6, 2),
    pr_Tstar = matrix(0, 6, 2),
    haz_Tstar = array(0, c(length(t), 6, 2))
  )

  local_mocked_bindings(
    marp_bstrp = function(...) mocked_bootstrap,
    .package = "marp"
  )

  expect_error(
    student_confint(
      n = 5, B = 2, t = t, m = 3, BB = 2,
      par_hat = rep(1, 12), mu_hat = rep(1, 6), pr_hat = rep(1, 6),
      haz_hat = matrix(1, length(t), 6), weights = rep(0, 6),
      alpha = 0.05, y = 120, best.model = 1, which.model = 1
    ),
    "No positive-weight candidate models"
  )
})
