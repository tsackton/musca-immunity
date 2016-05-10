#analyze phylostrata

#make dataset for analysis
mus.strat<-merge(mdom.all, difexp, by.x="gene", by.y="row.names", all=T)
mus.strat<-unique(mus.strat[!is.na(mus.strat$gene),c("gene", "protein", "ogsid", "class", "strata", "hmm", "sigp", "type.x", "dmel", "dmel.imm", "dmel.sym", "baseMean", "log2FoldChange", "pvalue", "padj")])
mus.strat$difexp=as.numeric(mus.strat$padj<0.05)*as.numeric(sign(mus.strat$log2FoldChange))
mus.strat<-merge(mus.strat, strata.age)
mus.strat$immune<-as.factor(!is.na(mus.strat$dmel.imm) | !is.na(mus.strat$hmm))
mus.strat<-merge(mus.strat, protlen, by.x="protein", by.y="isoform")

#proportion of young genes induced compared to old genes
#first, try a regular logistic regression

strata.logit<-glm(as.numeric(difexp==1) ~ min.age, family="binomial", data=mus.strat)

#also try with age ranks instead of ages, which treats age as a step function
strata.ordered<-glm(as.numeric(difexp==1) ~ rank(min.age), family="binomial", data=mus.strat)


#now generate corrected values of min age based on expression and length
#first, rescale
mus.strat$scaled.length=scale(log10(mus.strat$length), scale=T)
mus.strat$scaled.exp=mus.strat$baseMean
mus.strat$scaled.exp[mus.strat$baseMean==0]=NA
mus.strat$scaled.exp=scale(log10(mus.strat$scaled.exp),scale=T)
mus.strat$scaled.exp[is.na(mus.strat$scaled.exp)]=mean(mus.strat$scaled.exp,na.rm=T)
#compute coefs
exp.mod<-coef(lm(min.age ~ scaled.exp, data=mus.strat))[2]
len.mod<-coef(lm(min.age ~ scaled.length, data=mus.strat))[2]
#predicted values
mus.strat$min.age.norm<-mus.strat$min.age-(mus.strat$scaled.exp*exp.mod)
#look at plot to verify reasonableness
plot(mus.strat$min.age.norm ~ mus.strat$min.age, xlab="Raw age", ylab="Normalized age", col="gray30", cex=0.5)

#now use the normalized age in analysis too
strata.logit.norm<-glm(as.numeric(difexp==1) ~ as.vector(min.age.norm), family="binomial", data=mus.strat)

#also look at broader cats: ancient (>~900), old (~600-800), intermediate (~100-400), young(<~100) 
mus.strat$age.cat=cut(mus.strat$min.age, breaks=c(-400,(85+275)/2, (355+642)/2, (790+910)/2, 2000), label=c("young", "intermediate", "old", "ancient"), ordered_result=T)
mus.strat$age.cat.norm=cut(mus.strat$min.age.norm, breaks=c(-400,(85+275)/2, (355+642)/2, (790+910)/2, 2000), label=c("young", "intermediate", "old", "ancient"), ordered_result=T)

#in both cases chisq is highly significant
chisq.test(mus.strat$difexp==1, mus.strat$age.cat)
chisq.test(mus.strat$difexp==1, mus.strat$age.cat.norm)

#let's also do a logistic regression with our age categories
strata.logit.cat<-glm(difexp==1 ~ relevel(as.factor(as.character(age.cat)), ref="ancient"), family="binomial", data=mus.strat)
strata.logit.cat.norm<-glm(difexp==1 ~ relevel(as.factor(as.character(age.cat.norm)), ref="ancient"), family="binomial", data=mus.strat)

#do a post-hoc test of each cat vs. all others pooled
#this code creates a dataframe (strat.cat.res) that has the prop induced in each age category, 
#the chisq p-value for that age cat vs all others, and the CI of the proportion calculated by 
#prop.test()

strat.cat.res<-data.frame(cat=factor(c("young", "intermediate", "old", "ancient"), levels=c("young", "intermediate", "old", "ancient"), ordered=T), prop.l=numeric(4), prop=numeric(4), prop.u=numeric(4), prop.norm.l=numeric(4), prop.norm=numeric(4), prop.norm.u=numeric(4), p.value=numeric(4), p.value.norm=numeric(4))
for (i in c(1,2,3,4)) {
  strat.cat.res$p.value[i]<-chisq.test(mus.strat$difexp==1, mus.strat$age.cat==strat.cat.res$cat[i])$p.value
  strat.cat.res$p.value.norm[i]<-chisq.test(mus.strat$difexp==1, mus.strat$age.cat.norm==strat.cat.res$cat[i])$p.value
  strat.cat.res$prop[i]=with(droplevels(mus.strat[mus.strat$age.cat==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat))))$estimate
  strat.cat.res$prop.l[i]=with(droplevels(mus.strat[mus.strat$age.cat==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat)), conf.level=0.99))$conf.int[1]
  strat.cat.res$prop.u[i]=with(droplevels(mus.strat[mus.strat$age.cat==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat)), conf.level=0.99))$conf.int[2] 
  strat.cat.res$prop.norm[i]=with(droplevels(mus.strat[mus.strat$age.cat.norm==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat.norm))))$estimate
  strat.cat.res$prop.norm.l[i]=with(droplevels(mus.strat[mus.strat$age.cat.norm==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat.norm)), conf.level=0.99))$conf.int[1]
  strat.cat.res$prop.norm.u[i]=with(droplevels(mus.strat[mus.strat$age.cat.norm==strat.cat.res$cat[i],]), prop.test(t(table(difexp!=1, age.cat.norm)), conf.level=0.99))$conf.int[2] 
  
}

#BELOW HERE NOT RUN FOR PAPER##

#finally, let's look at the effect of gene age on annotated immune genes separately from genes without an immune annotation
#we'll repeat the logistic regression on data subsets first
strata.logit.norm.imm<-glm(as.numeric(difexp==1) ~ min.age.norm, family="binomial", data=mus.strat[mus.strat$immune==T,])
strata.logit.norm.ni<-glm(as.numeric(difexp==1) ~ min.age.norm, family="binomial", data=mus.strat[mus.strat$immune==F,])
strata.logit.imm<-glm(as.numeric(difexp==1) ~ min.age, family="binomial", data=mus.strat[mus.strat$immune==T,])
strata.logit.ni<-glm(as.numeric(difexp==1) ~ min.age, family="binomial", data=mus.strat[mus.strat$immune==F,])

#can also look at the ordered version for the non-normalized age
strata.ordered.imm<-glm(as.numeric(difexp==1) ~ rank(min.age), family="binomial", data=mus.strat[mus.strat$immune==T,])
strata.ordered.ni<-glm(as.numeric(difexp==1) ~ rank(min.age), family="binomial", data=mus.strat[mus.strat$immune==F,])

#finally, correlation between fold change and gene age
by(mus.strat, mus.strat$difexp, function(x) cor.test(x$min.age, x$log2FoldChange, method="k"))
by(mus.strat, mus.strat$difexp, function(x) cor.test(x$min.age.norm, x$log2FoldChange, method="k"))
