start.time <- Sys.time()
out <- sim(n = 50, mu = 3, sig = 1, P = 10, m = 10, t = seq(10,100), B = 99, R = 99, alpha = 0.05)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



library("xtable")

load(file = 'sims.RData')


attach()


par(font=2,font.axis=2,font.lab=2)


plot(t,rmse.haz.gen$Bias,type='l',lwd=2, ylim = c(-0.005,0.005),ylab = "Bias", xlab = "Times (t)")

lines(t,rmse.haz.best$Bias,lwd=2, col='red')

lines(t,rmse.haz.aic$Bias,lwd=2, col='blue')

lines(t,rmse.haz.prop$Bias,lwd=2, col='gold')

lines(t,rmse.haz.med.gen$Bias,lwd=2, col='green')

lines(t,rmse.haz.med.best$Bias,lwd=2, col='brown')


plot(t,rmse.haz.gen$Var,type='l',lwd=2, ylim = c(0,0.00012), ylab = "Var", xlab = "Times (t)")

lines(t,rmse.haz.best$Var,lwd=2, col='red')

lines(t,rmse.haz.aic$Var,lwd=2, col='blue')

lines(t,rmse.haz.prop$Var,lwd=2, col='gold')

lines(t,rmse.haz.med.gen$Var,lwd=2, col='green')

lines(t,rmse.haz.med.best$Var,lwd=2, col='brown')



plot(t,rmse.haz.gen$RMSE,type='l',lwd=2, ylim = c(0,0.011),ylab = "RMSE", xlab = "Times (t)")

lines(t,rmse.haz.best$RMSE,lwd=2, col='red')

lines(t,rmse.haz.aic$RMSE,lwd=2, col='blue')

lines(t,rmse.haz.prop$RMSE,lwd=2, col='gold')

lines(t,rmse.haz.med.gen$RMSE,lwd=2, col='green')

lines(t,rmse.haz.med.best$RMSE,lwd=2, col='brown')

dev.off()

rmse.name <- c("Generating","Best","AIC","Prop","Med.Gen","Med.Best")

rmse.mu <- data.frame(cbind(rmse.mu.gen,rmse.mu.best,rmse.mu.aic,rmse.mu.prop,rmse.mu.med.gen,rmse.mu.med.best))

colnames(rmse.mu) <- rmse.name

rmse.pr <- data.frame(cbind(rmse.pr.gen,rmse.pr.best,rmse.pr.aic,rmse.pr.prop,rmse.pr.med.gen,rmse.pr.med.best))

colnames(rmse.pr) <- rmse.name


xtable(rmse.mu)

xtable(rmse.pr)


cover.name <- c("Percentile Gen.","Percentile  Best", "Studentized Gen.", "Studentized Best", "Studentized M-A")

cover.mu <- data.frame(cbind(coverage.mu.percent.gen,coverage.mu.percent.best,coverage.mu.student.gen,coverage.mu.student.best,coverage.mu.student.ma))

colnames(cover.mu) <- cover.name


cover.pr <- data.frame(cbind(coverage.pr.percent.gen,coverage.pr.percent.best,coverage.pr.student.gen,coverage.pr.student.best,coverage.pr.student.ma))

colnames(cover.pr) <- cover.name

xtable(cover.mu)

xtable(cover.pr)


cover.haz <- sapply(1:5, function(i) cbind(unlist(coverage.haz.percent.gen[i,]),unlist(coverage.haz.percent.best[i,]),unlist(coverage.haz.student.gen[i,]),unlist(coverage.haz.student.best[i,]),unlist(coverage.haz.student.ma[i,])),simplify = 'array')

par(font=2,font.axis=2,font.lab=2,cex.axis= 0.7)

boxplot.matrix(cover.haz[,,1],xaxt='n',main = "Coverage Rate")

axis(1,at = 1:5, labels= cover.name)


par(mfrow = c(2,1),font=2,font.axis=2,font.lab=2,cex.axis= 0.7)

boxplot.matrix(cover.haz[,,2],xaxt='n',main = "Lower Error Rate")

axis(1,at = 1:5, labels= cover.name)

boxplot.matrix(cover.haz[,,3],xaxt='n',main = "Upper Error Rate")

axis(1,at = 1:5, labels= cover.name)


par(mfrow = c(2,1),font=2,font.axis=2,font.lab=2,cex.axis= 0.7)

boxplot.matrix(cover.haz[,,4],xaxt='n',main = "Lower Half Width")

axis(1,at = 1:5, labels= cover.name)

boxplot.matrix(cover.haz[,,5],xaxt='n',main = "Upper Half width")

axis(1,at = 1:5, labels= cover.name)

dev.off()
