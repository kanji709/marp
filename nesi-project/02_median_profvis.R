
options(keep.source=TRUE)

library(profvis)
library(marp)

dat <- read.table('../marp/data/small.txt')$V1 # rgamma(100,3,0.01)

m = 10 # number of iterations for MLE optimization 
t = seq(100,500,by=10) # time intervals 
# B = 999 # number of bootstraps
# BB = 999 # number of double-bootstrapps 
# alpha = 0.05 # confidence level 
y = 304 # cut-off year for estimating probablity 
# model_gen = 2 # underlying true model

# fitting differnt renewal models
p <- profvis({
    for (i in 1:1000) {
        res1 <- marp::poisson_rp(dat,t,y)
        res2 <- marp::gamma_rp(dat,t,m,y)
        res3 <- marp::loglogis_rp(dat,t,m,y)
        res4 <- marp::weibull_rp(dat,t,m,y)
        res5 <- marp::lognorm_rp(dat,t,y)
        res6 <- marp::bpt_rp(dat,t,m,y)
    }
})
htmlwidgets::saveWidget(p, "02_median_profvis.html")
