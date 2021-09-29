test_that("poisson_rp", {
    # load the small test dataset
    data_file <- system.file("extdata", "small.txt", package = "marp", mustWork = TRUE)
    data <- read.table(data_file)$V1

    # set some parameters
    m = 10 # number of iterations for MLE optimization
    t = seq(100,200,by=10) # time intervals
    y = 304 # cut-off year for estimating probablity

    # fix the random seed
    set.seed(42)

    # fit renewal model
    res <- marp::poisson_rp(data, t, y)

    # check result
    expect_equal(res$par1, 0.002818136384740506)
    expect_equal(res$par2, NA)
    expect_equal(res$logL, -206.1503840690702)
    expect_equal(res$AIC, 414.3007681381404)
    expect_equal(res$BIC, 415.7019655198025)
    expect_equal(res$mu_hat, 354.8444303174064)
    expect_equal(res$pr_hat, 0.3041016472842401)
    expect_true(all.equal(res$haz_hat, rep(c(-5.871679468969006), times=11)))
})
