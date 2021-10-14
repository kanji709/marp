test_that("pllog", {
  q <- c(1,2,3,4)
  shape <- 5
  scale <- 3
  log <- FALSE
  result <- pllog(q, shape, scale, log)
  expected_result <- as.numeric(c(0.9959016, 0.8836364, 0.5000000, 0.1917916))
  expect_true(all.equal(result, expected_result, tolerance = 1e-7))
  
  q <- c(1,2,3,4)
  shape <- 5
  scale <- 3
  log <- TRUE
  result <- pllog(q, shape, scale, log)
  expected_result <- as.numeric(c(0.004098361, 0.116363636, 0.500000000, 0.808208366))
  expect_true(all.equal(result, expected_result, tolerance = 1e-7))
})
