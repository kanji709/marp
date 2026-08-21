# marp

<!-- badges: start -->
![R CMD check](https://github.com/kanji709/marp/actions/workflows/check-standard.yml/badge.svg)
<!-- badges: end -->

`marp` fits six parametric renewal-process models and provides AIC-based model
selection and model-averaged estimates. It also implements percentile and
studentized bootstrap confidence intervals.

![Overview of the marp workflow](https://github.com/kanji709/marp/blob/master/inst/extdata/chart.png?raw=true)

## Installation

Install the released version from CRAN:

```r
install.packages("marp")
```

The development version is available from GitHub:

```r
# install.packages("remotes")
remotes::install_github("kanji709/marp")
```

## Basic workflow

```r
library(marp)

set.seed(42)
dat <- rgamma(100, shape = 3, rate = 0.01)

# m controls repeated random-start optimizations for candidate models that
# use nlm(); t contains the times for log-hazard evaluation.
m <- 10
t <- seq(100, 200, by = 10)
y <- 304

# Model codes are 1 Poisson, 2 Gamma, 3 log-logistic, 4 Weibull,
# 5 log-normal, and 6 Brownian passage time (BPT).
fit <- marp(dat, t, m, y, which.model = 2)
fit
summary(fit)
```

The fitted object retains the original named list components for backward
compatibility. For example, AIC weights and model-averaged estimates remain
available through `$`:

```r
fit$weights_AIC
fit$mu_aic
fit$pr_aic    # logit-transformed event probability at y
fit$haz_aic   # model-averaged log-hazards at t
```

Bootstrap confidence intervals use the original data and can be
computationally expensive, especially when both `B` and `BB` are large. The
standard S3 interface delegates to the existing bootstrap implementation:

```r
## Not run:
# ci <- confint(fit, data = dat, B = 99, BB = 99, level = 0.95)
# ci
```

See `vignette("marp-workflow")` for a fuller workflow and interpretation of
the returned quantities.
