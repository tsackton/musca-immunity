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
library(MASS)
library(vcd)
library(AER)
library(pscl)
library(plyr)

model<-real.mod.ogs
data<-model$data

#check the fit for poisson for all branches
#fit.stats<-ddply(data, .(branch), summarize, pois.gof=summary(goodfit(count))[1,3])
#Ord_plot(model$data$count[model$data$branch=="musDom"])

#dispersion
dispersiontest(model)
fit<-goodfit(data$count)
rootogram(fit)
summary(fit)

#mus only model
mus.mod<-glm(model$formula, data[data$branch=="musDom",], family="poisson")
dispersiontest(mus.mod)$p.value
disptest.branch<-ddply(data, .(branch), summarize, p.value=dispersiontest(glm(model$formula, family="poisson", data[data$branch==branch,]))$p.value)


#mixed model?
library(lme4)

sim.mixed<-glmer(count ~ 1 + (1|rate), offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$nodeclass != "root",])
real.mixed<-glmer(count ~ 1 + (1|ogs), offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))

real.mixed.immune<-glmer(count ~ 1 + type*musca*immune.narrow + (1|ogs), offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))

real.mixed.immune.nb<-glm.nb()
