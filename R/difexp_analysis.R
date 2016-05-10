#analysis of differential expression results

#basic summary
summary(des.rsem.results, alpha=0.05)

#annotated vs differetial expression
#first, make a subset of mdom.all for analysis that does not have duplication information

mdom.dif<-unique(subset(mdom.all, mdom.all$type.y == "total", select=c("gene", "strata", "hmm", "eval", "sigp", "type.x", "dmel", "dmel.imm", "dmel.sym")))
mdom.dif<-merge(mdom.dif, difexp, all=T, by.x="gene", by.y="row.names")
mdom.dif<-mdom.dif[!is.na(mdom.dif$gene),]
mdom.dif$difexp<-sign(mdom.dif$log2FoldChange)*(mdom.dif$padj<0.05)
#mdom dif includes non-coding RNAs which we have expression for but not homology annotations
#to get just protein-coding genes, ignore those where is.na(strata)==T

#combined both measures
mdom.dif$immune=F
mdom.dif$immune[!is.na(mdom.dif$dmel.imm) | !is.na(mdom.dif$hmm)]=T
#verify assignments
table(mdom.dif$dmel.imm, mdom.dif$hmm, mdom.dif$immune, useNA="ifany")

#new tests for revised manuscript:
with(mdom.dif[!is.na(mdom.dif$strata),], table(difexp, !is.na(dmel.imm))[3,2]/sum(!is.na(dmel.imm)))
with(mdom.dif[!is.na(mdom.dif$strata),], table(difexp, !is.na(hmm))[3,2]/sum(!is.na(hmm)))
with(mdom.dif[!is.na(mdom.dif$strata),], table(difexp, immune)[3,2]/sum(immune==T))

fisher.test(with(mdom.dif[!is.na(mdom.dif$strata),], table(difexp==1, !is.na(dmel.imm))))
fisher.test(with(mdom.dif[!is.na(mdom.dif$strata),], table(difexp==1, !is.na(hmm))))

##DMEL ANALYSIS
#dmel.imm has annotated dmel immunge genes
#dmel.hmm.imm has hmm results for dmel <- these have protein key not FBgn key

difexp_dmel = merge(difexp_dmel, dmel.imm, by.x="row.names", by.y="fbgn", all=T)
difexp_dmel$difexp = as.numeric(difexp_dmel$padj <= 0.05) * sign(difexp_dmel$log2FoldChange)
table(difexp_dmel$difexp, useNA="ifany")

table(difexp_dmel$difexp, !is.na(difexp_dmel$class))
sort(subset(difexp_dmel, difexp == 1 & !is.na(sym), select=c("sym"))$sym)

#dmel HMM analysis
#clean up dmel.hmm.imm
difexp_dmel = merge(difexp_dmel, dmel.hmm.imm, by.x="Row.names", by.y="gene", all.x=T)

#merge homologs
#first split multi-way dmel homologs into two lines
library(dplyr)
library(tidyr)

mdom.comp<-mdom.dif
mdom.comp<-mdom.comp %>% mutate(fbgn = strsplit(as.character(dmel), ",")) %>% unnest(fbgn)
mdom.comp<-merge(mdom.comp, difexp_dmel, by.x="fbgn", by.y="Row.names", suffixes=c(".mdom", ".dmel"))

#correlation of base expression
with(mdom.comp[mdom.comp$type.x=="1:1",], cor.test(baseMean.dmel, baseMean.mdom, method="kendall"))
with(mdom.comp[mdom.comp$type.x=="1:1",], cor.test(log2FoldChange.dmel, log2FoldChange.mdom, method="kendall"))
with(mdom.comp[mdom.comp$type.x=="1:1",], table(difexp.mdom, difexp.dmel))

#1:1 orthologs
mdom.comp$color="black"
mdom.comp$color[mdom.comp$difexp.dmel==1 & mdom.comp$difexp.mdom==1]="purple"
mdom.comp$color[mdom.comp$difexp.dmel==1 & mdom.comp$difexp.mdom<1]="blue"
mdom.comp$color[mdom.comp$difexp.dmel<1 & mdom.comp$difexp.mdom==1]="red"
with(mdom.comp[mdom.comp$type.x=="1:1" & mdom.comp$color=="black",], plot(log2FoldChange.mdom, log2FoldChange.dmel, col=color, cex=0.75, pch=16, xlim=c(-2.5,4.5), ylim=c(-2.5,4.5), xlab="House fly fold change", ylab="Drosophila fold change", bty="n"))
with(mdom.comp[mdom.comp$type.x=="1:1" & mdom.comp$color!="black",], points(log2FoldChange.mdom, log2FoldChange.dmel, col=color, cex=0.75, pch=16))
legend("topleft", legend=c("Induced in both", "Induced in house fly", "Induced in Drosophila"), col=c("purple", "red", "blue"), pch=16, bty="n")
with(mdom.comp[mdom.comp$type.x=="1:1",], venn(list(mdom.induced=gene[difexp.mdom == 1 & !is.na(difexp.mdom)],dmel.induced=gene[difexp.dmel==1 & !is.na(difexp.dmel)], dmel.repressed=gene[difexp.dmel==-1 & !is.na(difexp.dmel)],mdom.repressed=gene[difexp.mdom==-1 & !is.na(difexp.mdom)])))


#statistical tests
fisher.test(with(mdom.comp[mdom.comp$type.x=="1:1",], table(difexp.mdom > 0, difexp.dmel > 0)))
fisher.test(with(mdom.comp[mdom.comp$type.x=="1:1",], table(difexp.mdom < 0, difexp.dmel < 0)))

#test of mean fold change in dmel for genes not induced in dmel but induced in mdom
with(mdom.comp[mdom.comp$type.x=="1:1",], t.test(log2FoldChange.dmel[difexp.mdom == 1 & difexp.dmel == 0], log2FoldChange.dmel[difexp.mdom == 0 & difexp.dmel == 0]))

with(mdom.comp[mdom.comp$type.x=="1:1",], plot(density(log2FoldChange.dmel[difexp.mdom == 1 & difexp.dmel == 0], na.rm=T), col="red", lwd=2, main="", las=1, xlab="Drosophila log2 fold change", bty="n"))
with(mdom.comp[mdom.comp$type.x=="1:1",], lines(density(log2FoldChange.dmel[difexp.mdom == 0 & difexp.dmel == 0], na.rm=T), col="black", lwd=2))
legend("topleft", legend=c("Induced in house fly, not regulated in Drosophila", "Not regulated in either species"), col=c("red", "black"), lwd=2, bty="n", inset=c(0,-.1), xpd=NA)

#hmm analysis
