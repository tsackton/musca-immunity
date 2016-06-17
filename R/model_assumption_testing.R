#prep
setwd("~/Projects/immunity/musca/musca-immunity/R/")
source("load_data.R")
source("clean_duploss_data.R")
source("duploss_simulation_analysis.R")
all.sim.results<-read.table("model_check_results.200iter.out", sep="\t", header=T)

#packages
library(MASS)
library(vcd)
library(AER)
library(pscl)
library(plyr)
library(lme4)

#model testing

#fit several models:
sim.mod.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$rate==0.114 & sim.root.pois$nodeclass != "root",])
sim.mod.type<-glm(count ~ type, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type!="total" & sim.root.pois$rate==0.002 & sim.root.pois$nodeclass != "root",])

sim.mod2.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$nodeclass != "root",])
sim.mod2.type<-glm(count ~ type, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type!="total" & sim.root.pois$nodeclass != "root",])

real.mod.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))
real.mod.type<-glm(count ~ type, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))

real.mod.ogs<-glm(count ~ 1, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root" & ogs=="23876"))

#testing 
model<-sim.mod2.tot
data<-model$data

#dispersion
dispersiontest(model)
fit<-goodfit(data$count)
rootogram(fit)
summary(fit)

#basic issue: for an individual OGS, poisson is a good fit, and in the simulation data for a single rate class or with little rate variation the same is true, but for the full simulation set or for the full real set lots of dispersion.

#two possible approaches: 
#mixed-effects model with a random OGS effect
#negative binomial

#to test, first fit to simulated full dataset
sim.full.rand<-glmer(count ~ 1 + (1|ogs), offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$nodeclass != "root",])
summary(sim.full.rand)
sim.rand.est<-ranef(sim.full.rand)$ogs
sim.rand.est<-merge(sim.rand.est, unique(sim.root.pois[,c("ogs", "rate")]), by.x="row.names", by.y="ogs")
head(sim.rand.est)
#hmm seems pretty good
plot(-5.238+sim.rand.est[,2] ~ log(sim.rand.est[,3]*2), xlab="Real Rate (log scale)", ylab="Estimated Rate (log scale)")
abline(a=0, b=1)
cor.test((-5.238+sim.rand.est[,2]), log(sim.rand.est[,3]*2), method="ke")

sim.type.rand<-glmer(count ~ 1 + (1+type|ogs), offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type!="total" & sim.root.pois$nodeclass != "root",])
summary(sim.type.rand)
