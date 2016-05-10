#model testing

#fit several models:
sim.mod.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$rate==0.002 & sim.root.pois$nodeclass != "root",])
sim.mod.type<-glm(count ~ type, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type!="total" & sim.root.pois$rate==0.002 & sim.root.pois$nodeclass != "root",])

sim.mod2.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type=="total" & sim.root.pois$nodeclass != "root",])
sim.mod2.type<-glm(count ~ type, offset=log(br), family="poisson", data=sim.root.pois[sim.root.pois$type!="total" & sim.root.pois$nodeclass != "root",])

real.mod.tot<-glm(count ~ 1, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root"))
real.mod.type<-glm(count ~ type, offset=log(br), family="poisson", data=subset(all.pois, subset=="conserved" & type!="total" & nodeclass != "root"))

#testing 
model<-sim.mod.tot
data<-model$data

library(MASS)
library(vcd)

#check the fit for poisson for all branches
fit.stats<-ddply(data, .(branch), summarize, pois.gof=summary(goodfit(count))[1,3])
Ord_plot(model$data$count[model$data$branch=="musDom"])

#dispersion
library(AER)
