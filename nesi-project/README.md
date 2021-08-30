# marp optimisation

## Model fitting (small dataset)

Running just the model fitting with the 6 functions within a single loop with
1000 iterations (00_model_fitting_all_profvis.{R,html}). The table below shows
the proportion of time spent in each function. Times were taken from the
profvis output, profvis is a profiler. It can be used as follows:

```r
# this line is important to get source code info in the profiling output when
# running with Rscript
options(keep.source=TRUE)

library(profvis)

p <- profvis({
    # put the code here to be profiled
})

# the following line requires pandoc, make sure you load it (e.g. "module load pandoc")
htmlwidgets::saveWidget(p, "profile.html")

```

See output in [00_model_fitting_all_profvis.html](00_model_fitting_all_profvis.html)

- loglogis_rp takes most time
- then gamma_rp, weibull_rp, bpt_rp
- lognorm_rp and poisson_rp take very little time

### loglogis_rp

- 75% time spent in `stats::nlm`
  - 77% of the time spent in `stats::nlm` is spent in `dllog` (58% of `loglogis_rp`)
- 11% time spent in `median`

#### median

Noticed the call to `median` and `mean` of `data` occurs within a `while` loop but (I think)
`data` does not change within the loop. Therefore, we can precompute the mean and
median values and reuse them.

See output in [02_median_profvis.html](02_median_profvis.html)

- approx 15 % reduction in run time

#### dllog

Investigating use of Rcpp. Ongoing.


