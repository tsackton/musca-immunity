#analyze phylostrata

#make dataset for analysis
mus.strat<-merge(mdom.all, difexp, by.x="gene", by.y="row.names", all=T)
mus.strat<-unique(mus.strat[!is.na(mus.strat$gene),c("gene", "isoform", "ogsid", "class", "strata", "hmm", "sigp", "type.x", "dmel", "dmel.imm", "dmel.sym", "baseMean", "log2FoldChange", "pvalue", "padj")])
mus.strat$difexp=as.numeric(mus.strat$padj<0.05)*as.numeric(sign(mus.strat$log2FoldChange))
mus.strat<-merge(mus.strat, strata.age)
mus.strat$immune<-as.factor(!is.na(mus.strat$dmel.imm) | !is.na(mus.strat$hmm))
mus.strat<-merge(mus.strat, protlen)

#proportion of young genes induced compared to old genes
#first, try a regular logistic regression

strata.logit<-glm(as.numeric(difexp==1) ~ min.age, family="binomial", data=mus.strat)

#also try with age ranks instead of ages, which treats age as a step function
strata.ordered<-glm(as.numeric(difexp==1) ~ rank(min.age), family="binomial", data=mus.strat)

#finally, without order
strata.unordered<-glm(as.numeric(difexp==1) ~ relevel(strata,ref="Opisthokont"), family="binomial", data=mus.strat)

#now generate corrected values of min age based on expression and length
#first, rescale
mus.strat$scaled.length=scale(log10(mus.strat$length), scale=F)
mus.strat$scaled.exp=mus.strat$baseMean
mus.strat$scaled.exp[mus.strat$baseMean==0]=NA
mus.strat$scaled.exp=scale(log10(mus.strat$scaled.exp),scale=F)
mus.strat$scaled.exp[is.na(mus.strat$scaled.exp)]=mean(mus.strat$scaled.exp,na.rm=T)
#compute coefs
exp.mod<-coef(lm(min.age ~ scaled.exp, data=mus.strat))[2]
len.mod<-coef(lm(min.age ~ scaled.length, data=mus.strat))[2]
#predicted values
mus.strat$min.age.norm<-mus.strat$min.age+(mus.strat$scaled.exp*exp.mod)+(mus.strat$scaled.length*len.mod)
#look at plot to verify reasonableness
plot(mus.strat$min.age.norm ~ mus.strat$min.age, xlab="Raw age", ylab="Normalized age", col="gray30", cex=0.5)

#now use the normalized age in analysis too

