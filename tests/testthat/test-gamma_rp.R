test_that("gamma_rp", {
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
    res <- marp::gamma_rp(data, t, m, y)

    # check result
    expect_equal(res$par1, 3.584577307540262)
    expect_equal(res$par2, 0.01010182534216774)
    expect_equal(res$logL, -196.5790836929589)
    expect_equal(res$AIC, 397.1581673859178)
    expect_equal(res$BIC, 399.9605621492421)
    expect_equal(res$mu_hat, 354.8445143450731)
    expect_equal(res$pr_hat, -0.1696303553719428)

    haz_hat_expected <- c(-6.837250511298184, -6.680292814300809, -6.542946469927704, -6.421757960418091, -6.314055002359465, -6.217731760485419, -6.131102966898272, -6.052801955691350, -5.981707621905834, -5.916890979365526, -5.857575330456532)
    expect_true(all.equal(res$haz_hat, haz_hat_expected))
})
