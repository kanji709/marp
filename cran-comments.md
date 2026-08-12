## Resubmission

This package was previously on CRAN (0.1.0) and was archived on 2025-03-02
because a paper describing the methods used in this package, which we had
been given a deadline to add as a reference, was not published in time.

A paper describing these methods is being submitted to the R Journal, but
the R Journal requires the software to be on CRAN before it will consider
the submission. We will add the reference once the paper is published.

## Test environments

* local macOS installation, R 4.6.1
* via the rcmdcheck package on Linux (release, devel, and previous release),
  macOS (release), and Windows (release) (on GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission, package was previously archived (see "Resubmission" above).

## Longer examples and tests

* some of the examples take a few minutes each to run, these have been wrapped
  with `\donttest{}`
* likewise, longer running tests are skipped with `skip_on_cran()` and
  `skip_on_ci()`
