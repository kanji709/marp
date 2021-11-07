## Test environments

* local R installation, R 4.0.0
* via the rcmdcheck package on Linux, macOS and Windows with the current,
  development, and previous versions of R (on GitHub Actions)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* The GitHub repository pointed to by URL and BugReports is currently private
  but will be made public when the package is published

## Longer examples and tests

* some of the examples take a few minutes each to run, these have been wrapped
  with `\donttest{}`
* likewise, longer running tests are skipped with `skip_on_cran()` and
  `skip_on_ci()`

## References

* A paper that describes the methods used in this package is under preparation 
  and will include it as a reference once it has been published
