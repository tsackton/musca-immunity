#code to analyis gene family evolution using Poisson linear models, and also export data subsets for CAFE analysis

##organize data
#all.pois has the main dataset, need to add information about differential expression and immune function

#first get immune annotations
immune.ogs<-unique(mdom.all[,c("ogsid", "hmm", "dmel.imm")])
immune.ogs$hmm<-as.character(immune.ogs$hmm)
immune.ogs$hmm[is.na(immune.ogs$hmm)]="none"
immune.ogs$dmel.imm[is.na(immune.ogs$dmel.imm)]="none"
immune.ogs$dmel.imm[immune.ogs$dmel.imm=="modulation;signaling"] = "signaling"
immune.ogs$dmel.imm[immune.ogs$dmel.imm=="effector;recognition"] = "recognition"

#add to all.pois
all.pois<-merge(all.pois, immune.ogs, by.x="ogs", by.y="ogsid", all.x=T, all.y=F)
all.pois$dmel.imm = factor(all.pois$dmel.imm, levels=c("none", "recognition", "signaling", "modulation", "effector"))
all.pois$musca = factor(all.pois$musca, levels=c("other_dipt", "musca"))

#now run models on all.pois
#first, look at main results -- recognition, modulation, signaling, and effector differences on musca branch and overall
require(multcomp)

#rates for immune classes in whole phylogeny, assuming dup = loss; need to recode this as a ddply call to test many different options
imm.allclades<-glm(count ~ dmel.imm*type + offset(log(br)), family=poisson, data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass == "tip" & events>0))
summary(imm.allclades)

imm.byclades<-glm(count ~ type+dmel.imm*family + offset(log(br)), family=poisson, data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))
summary(imm.byclades)

imm.musca<-glm(count ~ dmel.imm*musca + offset(log(br)), family=poisson, data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root" & events>0))
summary(imm.musca)

#need to add contrasts to generate plots and test specific models


#now generate cafe subsets for each immune class, filtered and unfiltered, based on the root class
root.cafe.imm<-merge(root.cafe, immune.ogs, by.x="ogs", by.y="ogsid", all.x=T, all.y=F)
root.cafe.imm.filt<-subset(root.cafe.imm, events>0)

for (class in unique(root.cafe.imm$dmel.imm)) {
  write.table(file=paste0("../cafe/",class,"_filt",".tab"), subset(root.cafe.imm.filt, dmel.imm==class, select=c(2,1,3:16)), row.names=F, quote=F, sep="\t")
  write.table(file=paste0("../cafe/",class,"_all",".tab"), subset(root.cafe.imm, dmel.imm==class, select=c(2,1,3:16)), row.names=F, quote=F, sep="\t")  
}

