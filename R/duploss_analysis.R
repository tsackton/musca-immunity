setwd("~/Projects/immunity/musca/musca-immunity/R/")

#code to analyis gene family evolution using Poisson linear models, and also export data subsets for CAFE analysis

##organize data
#all.pois has the main dataset, need to add information about differential expression and immune function

#first get immune annotations
immune.ogs<-unique(mdom.all[,c("ogsid", "hmm", "dmel.imm")])
#need to collapse duplicated ogs
immune.ogs<-immune.ogs[order(immune.ogs$ogsid, immune.ogs$hmm, immune.ogs$dmel.imm),]
#NAs sort last
immune.ogs<-subset(immune.ogs, !duplicated(ogsid))
#fix a few specific cases
immune.ogs$dmel.imm[immune.ogs$ogsid=="14740.1"]="recognition"
immune.ogs$hmm[immune.ogs$ogsid=="8194"]="DIPT"

immune.ogs$hmm<-as.character(immune.ogs$hmm)
immune.ogs$hmm[is.na(immune.ogs$hmm)]="none"
immune.ogs$dmel.imm[is.na(immune.ogs$dmel.imm)]="none"
immune.ogs$dmel.imm[immune.ogs$dmel.imm=="modulation;signaling"] = "signaling"
immune.ogs$dmel.imm[immune.ogs$dmel.imm=="effector;recognition"] = "recognition"

#add to all.pois
all.pois<-merge(all.pois, immune.ogs, by.x="ogs", by.y="ogsid", all.x=T, all.y=F)
all.pois$dmel.imm = factor(all.pois$dmel.imm, levels=c("none", "recognition", "signaling", "modulation", "effector"))
all.pois$musca = factor(all.pois$musca, levels=c("other_dipt", "musca"))
all.pois$hmm = factor(all.pois$hmm)
all.pois$hmm = relevel(all.pois$hmm, ref="none")
all.pois = subset(all.pois, events > 0)
all.pois$immune.narrow = as.factor(all.pois$dmel.imm!="none")
all.pois$immune.broad = all.pois$immune.narrow
all.pois$immune.broad[all.pois$hmm != "none"] = TRUE


#now run models on all.pois
#first, look at main results -- recognition, modulation, signaling, and effector differences on musca branch and overall
require(multcomp)
require(lme4)

#result 1 -- musca vs other dipts for immune vs non-immune
musca.imm.test.dup<-glmer(count ~ musca*immune.narrow+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"))
musca.imm.test.loss<-glmer(count ~ musca*immune.narrow+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="loss" & nodeclass != "root"))
musca.imm.test.tot<-glmer(count ~  musca*immune.narrow+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))

#note convergence issues in total, use allFit in afex package to examine
library(afex)
library(optimx)
library(nloptr)
musca.imm.test.tot.conv<-allFit(musca.imm.test.tot)
is.OK<-sapply(musca.imm.test.tot.conv, is, "merMod")
musca.imm.test.tot.conv.ok<-musca.imm.test.tot.conv[is.OK]
lapply(musca.imm.test.tot.conv.ok, function(x) x@optinfo$conv$lme4$messages)

#convergence looks good using bobyqa, look at that result
summary(musca.imm.test.tot.conv.ok$bobyqa)

summary(musca.imm.test.dup)
summary(musca.imm.test.loss)

#broadly defined immune genes
musca.imm.test.dup.br<-glmer(count ~ musca*immune.broad+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"))
musca.imm.test.loss.br<-glmer(count ~ musca*immune.broad+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="loss" & nodeclass != "root"))
musca.imm.test.tot.br<-glmer(count ~  musca*immune.broad+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))

#total model also has convergence problems with broad immmune definition, use same approach as before

musca.imm.test.tot.br.conv<-allFit(musca.imm.test.tot.br)
is.OK<-sapply(musca.imm.test.tot.br.conv, is, "merMod")
musca.imm.test.tot.br.conv.ok<-musca.imm.test.tot.br.conv[is.OK]
lapply(musca.imm.test.tot.br.conv.ok, function(x) x@optinfo$conv$lme4$messages)

#convergence looks good using bobyqa, look at that result
summary(musca.imm.test.tot.br.conv.ok$bobyqa)

#is this driven by something funny going on in Drosophila (streamlined genomes)? Look at different families separately
#do this with family-based linear regression and specific contrasts
all.pois$family = factor(all.pois$family, levels=c("Muscidae", "Diptera", "Drosophilidae", "Glossinidae", "Culicidae"))

fam.imm.test.dup<-glmer(count ~ family*immune.narrow + offset(log(br)) + (1|ogs), family="poisson", data=droplevels(subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root" & family != "Diptera")), glmerControl(optimizer=c("bobyqa")))
fam.imm.test.loss<-glmer(count ~ family*immune.narrow + offset(log(br)) + (1|ogs), family="poisson", data=droplevels(subset(all.pois, subset=="conserved" & type=="loss" & nodeclass != "root" & family != "Diptera")), glmerControl(optimizer=c("bobyqa")))
fam.imm.test.tot<-glmer(count ~ family*immune.narrow + offset(log(br)) + (1|ogs), family="poisson", data=droplevels(subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root" & family != "Diptera")), glmerControl(optimizer=c("bobyqa")))

#now test by immune class
musca.immclass.test<-glmer(count ~ musca*dmel.imm+offset(log(br)) + (1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"), glmerControl(optimizer=c("bobyqa")))

#try other optimizers
musca.immclass.test.conv<-allFit(musca.immclass.test)
is.OK<-sapply(musca.immclass.test.conv, is, "merMod")
musca.immclass.test.conv.ok<-musca.immclass.test.conv[is.OK]
lapply(musca.immclass.test.conv.ok, function(x) x@optinfo$conv$lme4$messages)

#musca.immclass.test.conv.ok$optimx.nlminb

muscaVdipt_rec_dup=matrix(c(0,0,1,0,0,0,1,0,0,0)-c(0,0,1,0,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_rec_dup")))
muscaVdipt_sig_dup=matrix(c(0,0,0,1,0,0,0,1,0,0)-c(0,0,0,1,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_sig_dup")))
muscaVdipt_mod_dup=matrix(c(0,0,0,0,1,0,0,0,1,0)-c(0,0,0,0,1,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_mod_dup")))
muscaVdipt_eff_dup=matrix(c(0,0,0,0,0,1,0,0,0,1)-c(0,0,0,0,0,1,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_eff_dup")))
byclass_cont=rbind(muscaVdipt_rec_dup,muscaVdipt_eff_dup, muscaVdipt_mod_dup, muscaVdipt_sig_dup)
colnames(byclass_cont)=names(fixef(musca.immclass.test))
summary(glht(musca.immclass.test.conv.ok$optimx.nlminb, linfct=byclass_cont))
summary(musca.immclass.test.conv.ok$optimx.nlminb)

#rec, sig, mod, eff by family - duplication
musca.immclass.byfam<-glmer(count ~ family*dmel.imm+offset(log(br)) + (1|ogs), family="poisson", data=droplevels(subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root" & family != "Diptera")))

#convergence checks
musca.immclass.byfam.conv<-allFit(musca.immclass.byfam)
is.OK<-sapply(musca.immclass.byfam.conv, is, "merMod")
musca.immclass.byfam.conv.ok<-musca.immclass.byfam.conv[is.OK]
lapply(musca.immclass.byfam.conv.ok, function(x) x@optinfo$conv$lme4$messages)

#still no convergence so try rerunning with optimx.nlimnb which is closest
ss <- getME(musca.immclass.byfam.conv$optimx.nlminb,c("theta","fixef"))
musca.immclass.byfam.2 <- update(musca.immclass.byfam.conv$optimx.nlminb,start=ss,control=glmerControl(optCtrl=list(method="nlminb"),optimizer="optimx"))
  
byfam.contrasts<-matrix(data=0, nrow=12, ncol=length(fixef(musca.immclass.byfam.2)))
rownames(byfam.contrasts)=c("Dros_rec", "Glos_rec", "Cul_rec", "Dros_sig", "Glos_sig", "Cul_sig", "Dros_mod", "Glos_mod", "Cul_mod", "Dros_eff", "Glos_eff", "Cul_eff")
for (i in 1:12) {
  byfam.contrasts[i, i+8] = -1
}
colnames(byfam.contrasts)=names(fixef(musca.immclass.byfam.2))
summary(glht(musca.immclass.byfam.2, linfct=byfam.contrasts))

#analysis by HMM-defined gene family
musca.immhmm.test<-glmer(count ~ musca*hmm+offset(log(br))+(1|ogs), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"))
fams=c("BGBP", "CEC", "CLIPA", "CLIPB", "CLIPC", "CLIPD", "CTL", "FREP", "GALE", "HPX", "IGSF", "LYS", "MD2L", "NFKB", "NIM", "PGRP", "PPO", "SPRN", "SRCA", "SRCB", "SRCC", "TEP", "TLL", "TPX", "TSF")
hmm_conts<-matrix(0, ncol=length(names(coef(musca.immhmm.test))), nrow=length(fams), dimnames=list(fams,names(coef(musca.immhmm.test))))
hmm_conts[,3]=0
for (i in 1:25) {
  hmm_conts[i,i+54]=1  
}
summary(glht(musca.immhmm.test, linfct=hmm_conts))

#confirm with total, no type
musca.immhmmtot.test<-glm(count ~ musca*hmm, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))
hmm_conts2<-matrix(0, ncol=length(names(coef(musca.immhmmtot.test))), nrow=length(fams), dimnames=list(fams,names(coef(musca.immhmmtot.test))))
hmm_conts2[,2]=0
for (i in 1:25) {
  hmm_conts2[i,i+27]=1  
}
summary(glht(musca.immhmmtot.test, linfct=hmm_conts2))

#by lineage for a specific gene family
for (test in c("TEP", "LYS", "CEC")) {
  selected_fam=test
  musca.immhmm.family<-glm(count ~ type*family*(hmm==selected_fam), offset=log(br), family="poisson", data=droplevels(subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root" & family != "Diptera" & family != "Glossinidae")))
  #make tests
  fams=c("Dros", "Cul")
  hmm_conts<-matrix(0, ncol=length(names(coef(musca.immhmm.family))), nrow=length(fams), dimnames=list(fams,names(coef(musca.immhmm.family))))
  
  for (i in 1:nrow(hmm_conts)) {
    hmm_conts[i,i+8]=-1  
  }
  print(summary(glht(musca.immhmm.family, linfct=hmm_conts)))
}


#look at per-ogs rate
per.ogs.data<-droplevels(subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))
per.ogs.results<-data.frame(ogs=character(length(unique(per.ogs.data$ogs))), rate=numeric(length(unique(per.ogs.data$ogs))), pval=numeric(length(unique(per.ogs.data$ogs))), stringsAsFactors=F)
ogs.list<-as.character(unique(per.ogs.data$ogs))
for (i in 1:length(ogs.list)) {
  ogs<-ogs.list[i]
  ogs.data<-per.ogs.data[per.ogs.data$ogs == ogs,]
  per.ogs.results[i,1]=ogs
  per.ogs.results[i,c(2,3)]=coef(summary(glm(count ~ musca + offset(log(br)), family=poisson, data=ogs.data)))[2,c(1,4)]
}
per.ogs.results$qval = p.adjust(per.ogs.results$pval, method="fdr")

#add annotations
per.ogs.results.annot <- merge(per.ogs.results, unique(all.pois[,c("ogs", "hmm", "dmel.imm", "immune.narrow", "immune.broad")]), by="ogs")

#fisher tests
fisher.test(table(per.ogs.results.annot$qval < 0.05 & per.ogs.results.annot$rate > 0, per.ogs.results.annot$immune.broad))
table(per.ogs.results.annot$qval < 0.05 & per.ogs.results.annot$rate > 0, per.ogs.results.annot$immune.broad)

per.ogs.results.annot[per.ogs.results.annot$qval < 0.05 & per.ogs.results.annot$rate > 0 & per.ogs.results.annot$immune.broad==T,]


