## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
# devtools::install_github("kanji709/marp")

## ----data---------------------------------------------------------------------
if (requireNamespace("marp", quietly = TRUE)) {
  library(marp)
} else {
  devtools::load_all("..")
}

set.seed(42)
dat <- rgamma(50, shape = 3, rate = 0.01)
summary(dat)

## ----inputs-------------------------------------------------------------------
t <- seq(100, 200, by = 20)
y <- 304
m <- 3

## ----fit----------------------------------------------------------------------
set.seed(42)
fit <- marp(dat, t, m, y, which.model = 2)
fit

## ----summarize----------------------------------------------------------------
fit_summary <- summary(fit)
fit_summary

## ----extract------------------------------------------------------------------
fit$weights_AIC
fit$model_best
fit$mu_best
fit$mu_aic
fit$pr_aic
fit$haz_aic

## ----reference----------------------------------------------------------------
fit$mu_gen
fit$pr_gen
fit$haz_gen

## ----confidence-intervals, eval=FALSE-----------------------------------------
# set.seed(42)
# ci <- confint(
#   fit,
#   data = dat,
#   B = 99,
#   BB = 99,
#   level = 0.95
# )
# ci

## ----direct-confidence-intervals, eval=FALSE----------------------------------
# ci_direct <- marp_confint(
#   data = dat,
#   m = m,
#   t = t,
#   B = 99,
#   BB = 99,
#   alpha = 0.05,
#   y = y,
#   which.model = 2
# )

