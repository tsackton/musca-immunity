#generate estimates of the duplication and loss rate for each ogs separately

#rate1 is the calibrated tree
#rate2 is the uncalibrated tree, normalized to be on the same scale

require(plyr)
cons.rate<-ddply(cons.pois, .(ogs, type), summarize, rate1=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1], rate2=coef(summary(glm(count ~ offset(log(I(br*306.2357))), family=poisson)))[1,1])
root.rate<-ddply(root.pois, .(ogs, type), summarize, rate1=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1], rate1=coef(summary(glm(count ~ offset(log(I(br*306.2357))), family=poisson)))[1,1])
sch.rate<-ddply(sch.pois, .(ogs, type), summarize, rate1=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1], rate1=coef(summary(glm(count ~ offset(log(I(br*306.2357))), family=poisson)))[1,1])

cons.rate<-merge(cons.rate, cons.event)
root.rate<-merge(root.rate, root.event)
sch.rate<-merge(sch.rate, sch.event)

#estimates of per-unit-branchlength dup/loss rates
table(round(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]),3)/2)
table(round(exp(cons.rate$rate1[cons.rate$type=="total" & cons.rate$events > 0]),3)/2)
table(round(exp(sch.rate$rate1[sch.rate$type=="total" & sch.rate$events > 0]),3)/2)

#rates for simulation
sim.rates<-c(median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))/4, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))/2, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0])), median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*2, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*3, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*4, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*6, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*10, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*15, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*25, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*50, median(exp(root.rate$rate1[root.rate$type=="total" & root.rate$events > 0]))*150)

write.table(round(sim.rates,5), file="../data/rates.for.simulations", sep="\t", row.names=F, col.names=F, quote=F)
