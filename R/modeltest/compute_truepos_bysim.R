#test false positive rate from overdispersion with random samples
library(lme4)
library(MASS)

#load data
sim.root.pois<-read.table("sim.root.pois.data", header=T, sep="\t", stringsAsFactors=F)

#setup data
sim.ct.tot <- subset(sim.root.pois, type == "total")
sim.ct.dup <- subset(sim.root.pois, type == "dup")
sim.ct.loss <- subset(sim.root.pois, type=="loss")
sim.ogs<-data.frame(ogs=c(unique(sim.ct.tot$ogs)))
sim.ogs$rate = as.numeric(sub("sim\\.\\d+", "", sim.ogs$ogs, perl=T))
sim.ogs$lograte = log10(sim.ogs$rate)

SIMNUM<-1000
TRUECT<-100
results2<-list()

#set up list
for (method in c("poisson", "quasipoisson", "negbin", "mm1", "mm2")) {
  results2[[method]]<-data.frame(iter=seq(1,SIMNUM), real.coef=numeric(SIMNUM), est.coef=numeric(SIMNUM), stderr=numeric(SIMNUM), pval=numeric(SIMNUM))
  
}

for (i in seq(1,SIMNUM)) {
  #shuffle for regular
  sim.ogs$case=FALSE
  sim.ogs$case[sample(which(sim.ogs$rate>mean(sim.ogs$rate)), TRUECT)]=TRUE
  sim.input <- merge(sim.ct.tot, sim.ogs)
  sim.realcoef <- log(mean(sim.ogs$rate[sim.ogs$case==T])/mean(sim.ogs$rate[sim.ogs$case==F]))
  #glm poisson
  sim.pois <- summary(glm(count ~ case, offset=log(br), family="poisson", data=sim.input))
  #glmer multilevel model (random ogs)
  sim.ml1 <- summary(glmer(count ~ case + (1|ogs), offset=log(br), family="poisson", data=sim.input))
  sim.ml2 <- summary(glmer(count ~ case + (1|ogs) + (1|branch), offset=log(br), family="poisson", data=sim.input))
  #quaispoisson
  sim.qp <- summary(glm(count ~ case, offset=log(br), family="quasipoisson", data=sim.input))
  #negative binomial
  sim.nb <- summary(glm.nb(count ~ case + offset(log(br)), data=sim.input))
  
  #assign 
  results2[["poisson"]][i,]=c(i, sim.realcoef, coef(sim.pois)[2,c(1,2,4)])
  results2[["mm1"]][i,]=c(i, sim.realcoef, coef(sim.ml1)[2,c(1,2,4)])
  results2[["mm2"]][i,]=c(i, sim.realcoef, coef(sim.ml2)[2,c(1,2,4)])
  results2[["negbin"]][i,]=c(i, sim.realcoef, coef(sim.nb)[2,c(1,2,4)])
  results2[["quasipoisson"]][i,]=c(i, sim.realcoef, coef(sim.qp)[2,c(1,2,4)])
  
}

saveRDS(results2, file="perms.results2")
