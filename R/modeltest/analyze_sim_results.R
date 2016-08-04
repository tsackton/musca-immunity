
setwd("~/Projects/immunity/musca/musca-immunity/R/modeltest/")
results1<-readRDS("perms.results1")
results2<-readRDS("perms.results2")
all.fp<-do.call("rbind", results1)
all.fp$model <- sub("\\.\\d+", "", rownames(all.fp), perl=T)

table(all.fp$pval < 0.05, all.fp$model)
plot(real.coef ~ est.coef, data=all.fp[all.fp$model=="poisson",])
plot(real.coef ~ est.coef, data=all.fp[all.fp$model=="quasipoisson",])
plot(real.coef ~ est.coef, data=all.fp[all.fp$model=="negbin",])
plot(real.coef ~ est.coef, data=all.fp[all.fp$model=="mm1",])
plot(real.coef ~ est.coef, data=all.fp[all.fp$model=="mm2",])

all.tp<-do.call("rbind", results2)
all.tp$model <- sub("\\.\\d+", "", rownames(all.tp), perl=T)
table(all.tp$pval < 0.05, all.tp$model)

hist(all.tp$real.coef)
hist(all.fp$real.coef)

table(abs(all.fp$real.coef) > 0.5, all.fp$pval < 0.05, all.fp$model)

all.fp$error = sqrt((all.fp$est.coef - all.fp$real.coef)^2)

#error
plot(density(all.fp$error[all.fp$model=="mm1"]), col="red", ylim=c(0,4.5), main="RMSE for different models", lwd=2)
lines(density(all.fp$error[all.fp$model=="negbin"]), col="blue", lwd=2)
lines(density(all.fp$error[all.fp$model=="poisson"]), col="black", lwd=2)

legend("topright", legend=c("Mixed model", "Negative binomial", "(Quasi)poisson"), col=c("red", "blue", "black"), lwd=1)

#significance
barplot(table(all.fp$pval < 0.05, all.fp$model), las=1, ylab="# of permutations", xlab="Model", names=c("Mixed model 1", "Mixed model 2", "Neg. binom.", "Poisson", "Quasipoisson"))
legend("top", legend=c("P < 0.05", "P >= 0.05"), col=c("gray", "black"), pch=15, horiz=T, xpd=NA, inset=-0.15, bty="n")
abline(h=950, lwd=2, lty="dashed", col="red")

#only those with small real coef
barplot(table(all.fp$pval[abs(all.fp$real.coef) < 0.05] < 0.05, all.fp$model[abs(all.fp$real.coef) < 0.05])/length(all.fp$pval[abs(all.fp$real.coef) < 0.05 & all.fp$model=="poisson"]), las=1, ylab="# of permutations", xlab="Model", names=c("Mixed model 1", "Mixed model 2", "Neg. binom.", "Poisson", "Quasipoisson"), ylim=c(0,1))
abline(h=0.95, lwd=2, col="red", lty="dashed")
