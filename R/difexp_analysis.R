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

#these generate tables and FET for immune anntations vs. differential expression
#dmel annotation
with(mdom.dif[!is.na(mdom.dif$strata),], table(!is.na(dmel.imm), difexp))
fisher.test(with(mdom.dif[!is.na(mdom.dif$strata),], table(!is.na(dmel.imm), difexp))[,2:3])

#hmm
with(mdom.dif[!is.na(mdom.dif$strata),], table(!is.na(hmm), difexp))
fisher.test(with(mdom.dif[!is.na(mdom.dif$strata),], table(!is.na(hmm), difexp))[,2:3])

#combined both measures
mdom.dif$immune=F
mdom.dif$immune[!is.na(mdom.dif$dmel.imm) | !is.na(mdom.dif$hmm)]=T
#verify assignments
table(mdom.dif$dmel.imm, mdom.dif$hmm, mdom.dif$immune, useNA="ifany")

#combined table
with(mdom.dif[!is.na(mdom.dif$strata),], table(immune, difexp))
fisher.test(with(mdom.dif[!is.na(mdom.dif$strata),], table(immune, difexp))[,2:3])
