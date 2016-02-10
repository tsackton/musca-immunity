#FIGURES

#figure 1 - MA plot
library(DESeq2)
plotMA(object = des.rsem.results, alpha=0.05, main="", ylim=c(-2.5,2.5), cex=.6, las=1, bty="l")
legend("topleft", legend=c("Differentially regulated\n(5% FDR)"), bty="n", pch=16, col="red", cex=0.8)

#figure 2 - differential expression and log2 fold change by hmm class

#part A
barplot(prop.table(table(mdom.dif$difexp == 1, mdom.dif$hmm),2)[2,], col="red", las=2, ylab="Proportion Induced by Infection")

#part B
stripchart(mdom.dif$log2FoldChange[mdom.dif$padj > 0.05]  ~ mdom.dif$hmm[mdom.dif$padj > 0.05] , vertical=T, las=2, pch=16, col="black", method="jitter", jitter=0.1, ylim=c(-2,5), ylab="log2 fold change")
stripchart(mdom.dif$log2FoldChange[mdom.dif$padj <= 0.05]  ~ mdom.dif$hmm[mdom.dif$padj <= 0.05] , vertical=T, las=2, pch=16, col="red", method="jitter", jitter=0.1, add=T)

#figure 3 --  logistic regression results for uncorrected or corrected gene age
#partA
plot(with(mus.strat, prop.table(table(difexp==1, min.age),2))[2,] ~ sort(unique(mus.strat$min.age)), xlab="Gene Age (MYA)", ylab="Proportion Induced", las=1, pch=16, col="red", bty="l", cex=1.2)
points(predict(strata.logit, newdata=data.frame(min.age=seq(0,1400,10)), type="response") ~ seq(0,1400,10), type="l", lwd=3, col="black", lty="dashed")
legend(x=900, y=0.21, lwd=3, lty="dashed", legend=c("Logistic regression fit"), bty="n")

#partB
plot(with(mus.strat, prop.table(table(difexp==1, cut(min.age.norm, breaks=22, labels=F)),2))[2,] ~ c(unlist(by(as.vector(mus.strat[,c("min.age.norm")]), cut(mus.strat$min.age.norm, breaks=22, labels=F), median, na.rm=T))+308), xlab="Normalized Gene Age", ylab="Proportion Induced", las=1, pch=16, col="red", bty="l", cex=1.2)
points(predict(strata.logit.norm, newdata=data.frame(min.age.norm=seq(0,2000,10)), type="response") ~ seq(0,2000,10), type="l", lwd=3, col="black", lty="dashed")
legend(x=1300, y=0.29, lwd=3, lty="dashed", legend=c("Logistic regression fit"), bty="n")

#figure 4 -- proportion induced by gene category
require(plotrix)
plotCI(y=strat.cat.res$prop, x=c(1,2,3,4), ui=strat.cat.res$prop.u, li=strat.cat.res$prop.l, ylim=c(0,0.22), xlim=c(0.7,4.5), xlab="Category", ylab="Proportion Induced", las=1, bty="n", xaxt="n", pch=16, cex=1.2, col="darkgreen")
plotCI(y=strat.cat.res$prop.norm, x=c(1.2, 2.2, 3.2, 4.2), ui=strat.cat.res$prop.norm.u, li=strat.cat.res$prop.norm.l, add=TRUE, pch=16, cex=1.2, col="blue")
axis(1, at=c(1.1,2.1,3.1,4.1), labels=c("Young", "Intermediate", "Old", "Ancient"))
abline(h=mean(mus.strat$difexp==1, na.rm=T), lwd=2, lty="dotted")
legend(x=0.8, y=0.04, bty="n", legend=c("Raw Age", "Normalized Age"), col=c("darkgreen", "blue"), pch=16)
legend(x=0.8, y=0.02, bty="n", legend=c("P < 0.05"), col="black", pch=8)
mtext(2, at=0.08, text=("Genome\nAverage"), las=1)
mtext(3, at=c(1,1.2,2,2.2,4,4.2), text="*", cex=2, line=-4)
