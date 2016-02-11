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

#now run models on all.pois
#first, look at main results -- recognition, modulation, signaling, and effector differences on musca branch and overall
require(multcomp)

all.pois$immune.narrow = as.factor(all.pois$dmel.imm!="none")
all.pois$immune.broad = all.pois$immune.narrow
all.pois$immune.broad[all.pois$hmm != "none"] = TRUE

#result 1 -- musca vs other dipts for immune vs non-immune
musca.imm.test<-glm(count ~ type*musca*immune.narrow, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))
musca.imm.test.dup<-glm(count ~ musca*immune.narrow, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"))
musca.imm.test.loss<-glm(count ~ musca*immune.narrow, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="loss" & nodeclass != "root"))
musca.imm.test.tot<-glm(count ~ musca*immune.narrow, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))

musca.imm.test.br<-glm(count ~ type*musca*immune.broad, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))
musca.imm.test.dup.br<-glm(count ~ musca*immune.broad, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="dup" & nodeclass != "root"))
musca.imm.test.loss.br<-glm(count ~ musca*immune.broad, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="loss" & nodeclass != "root"))
musca.imm.test.tot.br<-glm(count ~ musca*immune.broad, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))

#first, type == dup, immune == T, difference between musca and other_dipt
musca_imm_dup<-matrix(c(1,0,1,1,0,0,1,0), nrow=1)
dipt_imm_dup<-matrix(c(1,0,0,1,0,0,0,0), nrow=1)
diff1=musca_imm_dup-dipt_imm_dup
#next, type == dup, immune == F, difference between musca and other_dipt
musca_ni_dup<-matrix(c(1,0,1,0,0,0,0,0), nrow=1)
dipt_ni_dup<-matrix(c(1,0,0,0,0,0,0,0), nrow=1)
diff2 = musca_ni_dup - dipt_ni_dup
#then difference between 1a and 1b
diff3 = diff1 - diff2
dup_test=rbind(diff1, diff2, diff3)
colnames(dup_test)=names(coef(musca.imm.test))
rownames(dup_test)=c("dup_musVsdipt_IMM", "dup_musVsdip_NI", "dup1-dup2")

#now for loss

#first, type == loss, immune == T, difference between musca and other_dipt
musca_imm_loss<-matrix(c(1,1,1,1,1,1,1,1), nrow=1)
dipt_imm_loss<-matrix(c(1,1,0,1,0,1,0,0), nrow=1)
diff1=musca_imm_loss-dipt_imm_loss
#next, type == loss, immune == F, difference between musca and other_dipt
musca_ni_loss<-matrix(c(1,1,1,0,1,0,0,0), nrow=1)
dipt_ni_loss<-matrix(c(1,1,0,0,0,0,0,0), nrow=1)
diff2 = musca_ni_loss - dipt_ni_loss
#then difference between 1a and 1b
diff3 = diff1 - diff2
loss_test=rbind(diff1, diff2, diff3)
colnames(loss_test)=names(coef(musca.imm.test))
rownames(loss_test)=c("loss_musVsdipt_IMM", "loss_musVsdip_NI", "loss1-loss2")

#final test
cont_test=rbind(dup_test, loss_test)
cont_res_narrow=glht(musca.imm.test, linfct=cont_test)
cont_res_broad=glht(musca.imm.test.br, linfct=cont_test)
summary(cont_res_narrow)
summary(cont_res_broad)


#now test by immune class
musca.immclass.test<-glm(count ~ type*musca*dmel.imm, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))
muscaVdipt_rec_dup=matrix(c(1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),nrow=1)-matrix(c(1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_rec_dup")))
muscaVdipt_sig_dup=matrix(c(1,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),nrow=1)-matrix(c(1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_sig_dup")))
muscaVdipt_mod_dup=matrix(c(1,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0),nrow=1)-matrix(c(1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_mod_dup")))
muscaVdipt_eff_dup=matrix(c(1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1)-matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=1, dimnames=list(c("muscaVdipt_eff_dup")))
byclass_cont=rbind(muscaVdipt_rec_dup,muscaVdipt_eff_dup, muscaVdipt_mod_dup, muscaVdipt_sig_dup)
colnames(byclass_cont)=names(coef(musca.immclass.test))
summary(glht(musca.immclass.test, linfct=byclass_cont))


#for each gene familiy

musca.immhmm.test<-glm(count ~ type*musca*hmm, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))
#make tests
fams=c("BGBP", "CEC", "CLIPA", "CLIPB", "CLIPC", "CLIPD", "CTL", "FREP", "GALE", "HPX", "IGSF", "LYS", "MD2L", "NFKB", "NIM", "PGRP", "PPO", "SPRN", "SRCA", "SRCB", "SRCC", "TEP", "TLL", "TPX", "TSF")
hmm_conts<-matrix(0, ncol=length(names(coef(musca.immhmm.test))), nrow=length(fams), dimnames=list(fams,names(coef(musca.immhmm.test))))
hmm_conts[,3]=1
for (i in 1:25) {
  hmm_conts[i,i+54]=1  
}

summary(glht(musca.immhmm.test, linfct=hmm_conts))

#confirm with total, no type
musca.immhmmtot.test<-glm(count ~ musca*hmm, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))
hmm_conts2<-matrix(0, ncol=length(names(coef(musca.immhmmtot.test))), nrow=length(fams), dimnames=list(fams,names(coef(musca.immhmmtot.test))))
hmm_conts2[,2]=1
for (i in 1:25) {
  hmm_conts2[i,i+27]=1  
}
summary(glht(musca.immhmmtot.test, linfct=hmm_conts2))

#look at per-ogs rate
cons.pois$musca=factor(cons.pois$musca,levels=c("other_dipt", "musca"))
cons.rate.musca<-ddply(cons.pois, .(ogs, type), summarize, rate=coef(summary(glm(count ~ musca + offset(log(br)), family=poisson)))[2,1], pval=coef(summary(glm(count ~ musca + offset(log(br)), family=poisson)))[2,4])

cons.rate.musca.tot = cons.rate.musca[cons.rate.musca$type=="total",]
cons.rate.musca.tot$qval = p.adjust(cons.rate.musca.tot$pval, method="fdr")

#add annotations
cons.rate.musca.tot <- merge(cons.rate.musca.tot, unique(all.pois[,c("ogs", "hmm", "dmel.imm", "immune.narrow", "immune.broad")]), by="ogs")

