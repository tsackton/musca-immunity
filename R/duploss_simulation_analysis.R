#process simulated tree data
#requires objects from clean_duploss_data.R so run that first
#focus on just the root subset for now

sim.ct<-read.table("../input_data/simulations/sim_dup_ct.txt", header=T, stringsAsFactors=F)

#clean up for analysis
require(reshape2)

#look at three different data sets

#the first set is all groups inferrred to be present at all four major internal branches (edges 1, 18, 2, and 5 -- so in the ancestors of each family, roughly)
#define based on total counts
sim.classes<-sim.ct[sim.ct$type=="total",]
sim.classes$cons = F
sim.classes$cons[sim.classes$X1 != "NP" & sim.classes$X6 != "NP" & sim.classes$X7 != "NP" & sim.classes$X8 != "NP" & sim.classes$X2 != "NP"] = T

#second is just all genes present at the root, so this will include some mos/dros only genes or mos/mus only genes
sim.classes$root = F
sim.classes$root[sim.classes$X1 != "NP"] = T

sim.root<-droplevels(sim.ct[sim.ct$ogs %in% sim.classes$ogs[sim.classes$root==T],])
sim.root<-replace(sim.root, sim.root=="NP", 0)
sim.ct.root<-melt(sim.root, id.vars=c(1,2), variable.name="treefix", value.name="count")
sim.ct.root=droplevels(sim.ct.root[sim.ct.root$treefix != "total" & sim.ct.root$treefix!="X1",])

sim.root.pois<-merge(sim.ct.root, ultra.edges, by.x="treefix", by.y="treefix")
sim.root.pois$count<-as.numeric(sim.root.pois$count)
sim.root.pois$rate = as.numeric(sub("sim\\.\\d+", "", sim.root.pois$ogs, perl=T))

#counts of analyzed genes in each simulation bin
table(sim.root.pois$rate[sim.root.pois$type=="dup"])/length(unique(sim.root.pois$edge))

#now look at each rate subset, in a few ways
require(plyr)
sim.root.pois.tot = subset(sim.root.pois, type == "total")
sim.root.pois.notot = subset(sim.root.pois, type != "total")
duploss_fit <- function(dat) {
  fit <- glm(count ~ type + offset(log(br)), family="poisson", data=dat)
  setNames(coef(fit), c("dup.rate", "loss.rate"))
}
tot_fit <- function(dat) {
  fit <- glm(count ~ offset(log(br)), family="poisson", data=dat)
  setNames(coef(fit), c("tot.rate"))
}
sim.rate.est.all<-ddply(sim.root.pois.notot, ~ rate, duploss_fit)
sim.rate.est.all$loss.rate = sim.rate.est.all$dup.rate + sim.rate.est.all$loss.rate
sim.rate.est.tot<-ddply(sim.root.pois.tot, ~ rate, tot_fit)

#make final simulation dataset
sim.rates<-merge(sim.rate.est.all, sim.rate.est.tot, by="rate")
sim.rates$log = log(sim.rates$rate)
sim.rates$tot.rmse = sqrt((exp(sim.rates$tot.rate) - (sim.rates$rate*2))^2)
sim.rates$dup.rmse = sqrt((exp(sim.rates$dup.rate) - (sim.rates$rate))^2)
sim.rates$loss.rmse = sqrt((exp(sim.rates$loss.rate) - (sim.rates$rate))^2)

cafe.rates<-read.table("../input_data/simulations/sim_cafe_results.txt", header=F, stringsAsFactors=F)
names(cafe.rates)<-c("rate", "cafe.dup", "cafe.loss", "cafe.tot", "sim.rate")

sim.rates<-merge(sim.rates, cafe.rates)
sim.rates$pois.dup=exp(sim.rates$dup.rate)
sim.rates$pois.los=exp(sim.rates$loss.rate)
sim.rates$pois.tot=exp(sim.rates$tot.rate)/2
