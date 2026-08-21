## Resubmission

This package was previously on CRAN (0.1.0) and was archived on 2025-03-02
because a paper describing the methods used in this package, which we had
been given a deadline to add as a reference, was not published in time.

A paper describing these methods will be submitted to the R Journal, but
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

* examples that invoke the nested bootstrap can take substantially longer than
  five seconds, so these are wrapped with `\dontrun{}`; the ordinary fitting
  and S3-interface examples run during checks
* likewise, longer running tests are skipped with `skip_on_cran()` and
  `skip_on_ci()`; fast tests cover bootstrap orchestration and returned
  structures on CRAN and CI
