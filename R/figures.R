#FIGURES

#figure 1 - MA plot
library(DESeq2)
cairo_pdf(file="Figure1.pdf", width=3.25, height=2.5)
par(mar=c(4.1,4.1,1,1), cex=0.5)
plotMA(object = des.rsem.results, alpha=0.05, main="", ylim=c(-2.5,2.5), cex=.75, las=1, bty="l", cex.lab=1.1)
legend("topleft", legend=c("Differentially regulated\n(5% FDR)"), bty="n", pch=16, col="red", cex=1)
dev.off()

##################
#figure 2 - differential expression and log2 fold change by hmm class

pdf(file="Figure2.pdf", width=6, height=5)
par(mar=c(3.1,4.1,1.8,1), cex=0.4, mfcol=c(2,1))

#order HMM class
decreasing.prop<-names(sort(prop.table(table(mdom.dif$difexp == 1, mdom.dif$hmm),2)[2,], decreasing=T))
eff.classes = c("DEF", "ATT", "CEC", "DIPT", "TSF", "LYS", "PPO", "TPX", "HPX", "GPX")
rec.classes = c("BGBP", "TEP", "PGRP", "GALE", "FREP", "CTL", "NIM", "MD2L", "SRCA", "SRCB", "SRCC")
sig.classes = c("CLIPE", "SPRN", "CLIPC", "CLIPA", "CLIPB", "TLL", "NFKB", "CLIPD")

mdom.dif$hmm2 = factor(mdom.dif$hmm, levels=c(decreasing.prop[decreasing.prop %in% rec.classes], decreasing.prop[decreasing.prop %in% sig.classes], decreasing.prop[decreasing.prop %in% eff.classes]))
difexp_dmel$hmm2 = factor(difexp_dmel$hmm, levels=c(decreasing.prop[decreasing.prop %in% rec.classes], decreasing.prop[decreasing.prop %in% sig.classes], decreasing.prop[decreasing.prop %in% eff.classes]))

colors=c(rep("darkgreen", length(rec.classes)), rep("purple", length(sig.classes)), rep("red", length(eff.classes)))

#part A
barplot(prop.table(table(mdom.dif$difexp == 1, mdom.dif$hmm2),2)[2,], col=colors, las=2, ylab="Proportion Induced", cex.lab=0.75, cex.axis=0.75, cex.names=0.75, cex.main=0.85, ylim=c(0,1), main=expression(italic("M. domestica")))
mtext(expression(bold("A")), side=2, line=3, at=1.06, las=2)
legend(x=29, y=1.15, bty="n", pch=15, legend=c("Recognition", "Signaling", "Effector"), col=c("darkgreen", "purple", "red"), cex=0.8, xpd=NA)
#part B
barplot(prop.table(table(difexp_dmel$difexp == 1, difexp_dmel$hmm2),2)[2,], col=colors, las=2, ylab="Proportion Induced", cex.lab=0.75, cex.axis=0.75, cex.names=0.75, ylim=c(0,1), cex.main=0.85, main=expression(italic("D. melanogaster")))
mtext(expression(bold("B")), side=2, line=3, at=1.06, las=2)

dev.off()


#############
#figure 3 -- 1:1 orthology comparison
mdom.comp$color="black"
mdom.comp$color[mdom.comp$difexp.dmel==1 & mdom.comp$difexp.mdom==1]="purple"
mdom.comp$color[mdom.comp$difexp.dmel==1 & mdom.comp$difexp.mdom<1]="blue"
mdom.comp$color[mdom.comp$difexp.dmel<1 & mdom.comp$difexp.mdom==1]="red"

pdf(file="Figure3.pdf", width=6, height=3, onefile=FALSE)
# setup layout
library(gridBase)
library(VennDiagram)
gl <- grid.layout(nrow=1, ncol=2)

# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1) 

# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)

# start new base graphics in first viewport
par(fig=gridFIG())
par(mar=c(3.8,3.8,1,1), cex=0.6)
#A
with(mdom.comp[mdom.comp$type.x=="1:1" & mdom.comp$color=="black",], plot(log2FoldChange.mdom, log2FoldChange.dmel, col=color, cex=0.9, pch=16, xlim=c(-2.5,4.5), ylim=c(-2.5,4.5), xlab="", ylab="", bty="n", cex.axis=0.9, cex.lab=0.90, las=1))
with(mdom.comp[mdom.comp$type.x=="1:1" & mdom.comp$color!="black",], points(log2FoldChange.mdom, log2FoldChange.dmel, col=color, cex=0.9, pch=16))
legend(y=4, x=-2.5, legend=c("Induced in both", "Induced in house fly", "Induced in Drosophila"), col=c("purple", "red", "blue"), pch=16, bty="n", cex=0.9, xpd=NA)
mtext(side=1, text="House fly fold change", line=2.2, cex=0.5)
mtext(side=2, text="Fruit fly fold change", line=2.2, cex=0.5)
popViewport()
pushViewport(vp.2)
#B
with(mdom.comp[mdom.comp$type.x=="1:1",], draw.pairwise.venn(length(gene[difexp.mdom == 1 & !is.na(difexp.mdom)]),length(gene[difexp.dmel==1 & !is.na(difexp.dmel)]), length(intersect(gene[difexp.mdom == 1 & !is.na(difexp.mdom)],gene[difexp.dmel==1 & !is.na(difexp.dmel)])), fill=c("red", "blue"), lty=2, alpha=c(0.7,0.7),height=3, width=3, category=c("House fly \ninduced", "Fruit fly \ninduced"), margin=0.05, cat.cex=c(0.6,0.6), ext.text=F, cat.dist=0.05, ind=T, cex=0.7, cat.pos=c(300,150)))
popViewport(1)
dev.off()

##########
#figure 4 --  logistic regression results for uncorrected or corrected gene age

cairo_pdf("Figure4.pdf", width=6.75, height=4)
par(mar=c(4.2,4.1,1.8,2), cex=0.4, mfrow=c(1,3))

#A
plot(with(mus.strat, prop.table(table(difexp==1, min.age),2))[2,] ~ sort(unique(mus.strat$min.age)), xlab="Gene Age (MYA)", ylab="Proportion Induced", las=1, pch=16, col="red", bty="l", cex=1.2)
points(predict(strata.logit, newdata=data.frame(min.age=seq(0,1400,10)), type="response") ~ seq(0,1400,10), type="l", lwd=2, col="black", lty="dashed")
legend(x=600, y=0.185, lwd=1, lty="dashed", legend=c("Logistic\n regression fit"), bty="n", cex=0.9)
mtext(2, at=0.19, text=expression(bold("A")), line=1, las=1)

#B
plot(with(mus.strat, prop.table(table(difexp==1, cut(min.age.norm, breaks=22, labels=F)),2))[2,] ~ c(unlist(by(as.vector(mus.strat[,c("min.age.norm")]), cut(mus.strat$min.age.norm, breaks=22, labels=F), median, na.rm=T))+308), xlab="Normalized Gene Age", ylab="Proportion Induced", las=1, pch=16, col="red", bty="l", cex=1.2)
points(predict(strata.logit.norm, newdata=data.frame(min.age.norm=seq(0,2000,10)), type="response") ~ seq(0,2000,10), type="l", lwd=2, col="black", lty="dashed")
legend(x=900, y=0.21, lwd=1, lty="dashed", legend=c("Logistic\n regression fit"), bty="n", cex=0.9)
mtext(2, at=0.225, text=expression(bold("B")), line=1, las=1)

#C -- proportion induced by gene category
require(plotrix)
plotCI(y=strat.cat.res$prop, x=c(1,2,3,4), ui=strat.cat.res$prop.u, li=strat.cat.res$prop.l, ylim=c(0,0.22), xlim=c(0.7,4.5), xlab="Category", ylab="Proportion Induced", las=1, bty="l", xaxt="n", pch=16, cex=1.2, col="darkgreen", cex.lab=0.9)
plotCI(y=strat.cat.res$prop.norm, x=c(1.2, 2.2, 3.2, 4.2), ui=strat.cat.res$prop.norm.u, li=strat.cat.res$prop.norm.l, add=TRUE, pch=16, cex=1.2, col="blue")
axis(1, at=c(1.1,2.1,3.1,4.1), labels=c("Young", "Inter.", "Old", "Ancient"))
abline(h=mean(mus.strat$difexp==1, na.rm=T), lwd=2, lty="dotted")
legend(x=0.8, y=0.04, bty="n", legend=c("Raw Age", "Normalized Age"), col=c("darkgreen", "blue"), pch=16)
legend(x=0.8, y=0.02, bty="n", legend=c("P < 0.01"), col="black", pch=8)
mtext(4, at=0.085, text=("Genome\nAverage"), las=1, line=-1, cex=0.55)
mtext(3, at=c(1,1.2,4,4.2), text="*", cex=1.5, line=-4)
mtext(2, at=0.225, text=expression(bold("C")), line=1, las=1)

dev.off()

#figure 5 (phylogeny)
library(ape)
pdf(file="Figure5.pdf", width=2, height=3)
par(mar=c(2,0.5,0.1,0.1))
dipt.plot<-dipt.ultra
species.key<-data.frame(shortname=c("musDom", "gloMor", "droWil", "droPse", "droYak", "droMel", "droAna", "droVir", "droMoj", "culQui", "aedAeg", "anoDar", "anoSte", "anoGam"), species=c("M. domestica", "G. morsitans", "D. willistoni", "D. pseudoobscura", "D. yakuba", "D. melanogaster", "D. ananassae", "D. virilis", "D. mojavensis", "C. quinquefasciatus", "A. aegypti", "A. darlingi", "A. stephensi", "A. gambiae"), stringsAsFactors=F)
rownames(species.key)=species.key$shortname
dipt.plot$tip.label = species.key[dipt.plot$tip.label,"species"]
plot(dipt.plot, cex=0.5, no.margin=F)
axisPhylo(cex=0.5, cex.axis=0.5, cex.lab=0.5, mgp=c(1,0.3,0), line=0.5)
dev.off()

#########
#figure 6

pdf(file="Figure6.pdf", height=4, width=8)
par(mar=c(4.1,4.1,1.8,1), cex=0.4, mfcol=c(1,2))

#A
plot(exp(sim.rates$log) ~ sim.rates$rate, col="red", type="l", lwd=1, log="xy", ylab="Estimated Rate (Events / MY)", xlab="Simulated Rate (Events / MY)", bty="l", ylim=c(0.0005, 0.4), xlim=c(0.0005, 0.4), las=2, cex.axis=0.7, cex.lab=0.7)
points(exp(sim.rates$tot.rate)/2~ sim.rates$rate, col="black", pch=0)
points(exp(sim.rates$dup.rate) ~ sim.rates$rate, col="blue", pch=0)
points(exp(sim.rates$loss.rate) ~ sim.rates$rate, col="darkgreen", pch=0)
legend(x=5e-04, y=0.13, legend=c("simulation expectation"), lwd=1, col="red", bty="n", cex=0.7)
legend(x=5.3e-04, y=0.60, legend=c("Duplication", "Loss", "Total"), pch=0, col=c("blue", "darkgreen", "black"), title="Poisson", bty="n", cex=0.7)
mtext(side=2, at=0.8, text=expression(bold("A")), las=1, line=3)

plot(exp(sim.rates$log) ~ sim.rates$rate, col="red", type="l", lwd=1, log="xy", ylab="Estimated Rate (Events / MY)", xlab="Simulated Rate (Events / MY)", bty="l", ylim=c(0.0005, 0.4), xlim=c(0.0005, 0.4), las=2, cex.axis=0.7, cex.lab=0.7)
points(sim.rates$cafe.tot ~ sim.rates$rate, col="black", pch = 2)
points(sim.rates$cafe.dup ~ sim.rates$rate, col="blue", pch = 2)
points(sim.rates$cafe.loss ~ sim.rates$rate, col="darkgreen", pch = 2)
legend(x=5e-04, y=0.13, legend=c("simulation expectation"), lwd=1, col="red", bty="n", cex=0.7)
legend(x=5.3e-04, y=0.60, title="CAFE", legend=c("Duplication", "Loss", "Total"), pch=2, col=c("blue", "darkgreen", "black"), bty="n", cex=0.7)
mtext(side=2, at=0.8, text=expression(bold("B")), las=1, line=3)

dev.off()
######

#figure 7
pdf(file="Figure7.pdf", width=4, height=3)
par(cex=0.7, mar=c(3,4,3,1))
library(plotrix)
fig7data<-summary(glht(musca.immclass.test.conv.ok$optimx.nlminb, linfct=byclass_cont))$test
fig7fam<-summary(glht(musca.immclass.byfam.2, linfct=byfam.contrasts))$test
plotCI(x=c(1,1.4,1.8,2.2), y=fig7data$coefficients[c(1,4,3,2)], uiw=fig7data$sigma[c(1,4,3,2)], col=c("darkgreen", "purple", "blue", "red"), pch=16, xlab="", xaxt="n", bty="n", ylim=c(-2,3), xlim=c(0.8,8.4), ylab="Duplication Rate (Musca vs Other Dipterans)", las=2, cex=1.2)
axis(1, at=c(1.5,3.5,5.5,7.5), labels=c("All Dipterans", "vs Mosquito", "vs Glossina", "vs Drosophila"), lwd=0, lwd.ticks=0, mgp=c(0,0.5,0))
mtext(1, line=1.9, text="Comparison", cex=0.75)
abline(h=0, lwd=2, lty="dotted")
plotCI(x=c(3,3.4,3.8,4.2), y=fig7fam$coef[c(3,6,9,12)], uiw=fig7fam$sigma[c(3,6,9,12)], col=c("darkgreen", "purple", "blue", "red"), pch=16, cex=1.2, add=T)
plotCI(x=c(5,5.4,5.8,6.2), y=fig7fam$coef[c(2,5,8,11)], uiw=fig7fam$sigma[c(2,5,8,11)], col=c("darkgreen", "purple", "blue", "red"), pch=16, cex=1.2, add=T)
plotCI(x=c(7,7.4,7.8,8.2), y=fig7fam$coef[c(1,4,7,10)], uiw=fig7fam$sigma[c(1,4,7,10)], col=c("darkgreen", "purple", "blue", "red"), pch=16, cex=1.2, add=T)
abline(v=2.5)
abline(v=4.5)
abline(v=6.5)
legend(x=0.8, y=4.2, legend=c("Recognition", "Signaling", "Modulation", "Effector"), pch=16, col=c("darkgreen", "purple", "blue", "red"), bty="n", bg="white", ncol=4, box.lwd=0, xpd=NA, cex=0.8)

positions=c(1,1.4,1.8,2.2,3,5,7,3.4,5.4,7.4,3.8,5.8,7.8,4.2,6.2,8.2)
p01=c(fig7data$pvalues[c(1,4,3,2)] < 0.01, fig7fam$pvalues < 0.01)
p05=c(fig7data$pvalues[c(1,4,3,2)] < 0.05, fig7fam$pvalues < 0.05)
p05[p01]=F

mtext(at=positions[p01], line=-1.3, text="**", cex=0.8)
mtext(at=positions[p05], line=-1.3, text="*", cex=0.8)

dev.off()
