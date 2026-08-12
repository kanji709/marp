# marp 0.1.1

* Resubmission to CRAN following archival of 0.1.0.
* Robustified BPT model fitting and confidence interval handling, including
  handling of `NA` results from `uniroot()`, exclusion of invalid bootstrap
  replications (with the percentage removed reported), rejection of extreme
  BPT fitted means, and a cap on the number of fit retries after blow-up.
* Documentation and example fixes (typos, return value descriptions, example
  runtime tweaks).

# marp 0.1.0

* Initial release of the package.
* Added a `NEWS.md` file to track changes to the package.
