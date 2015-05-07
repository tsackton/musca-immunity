#generate estimates of the duplication and loss rate for each ogs separately

require(plyr)
cons.rate<-ddply(cons.pois, .(ogs, type), summarize, rate=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1])
root.rate<-ddply(root.pois, .(ogs, type), summarize, rate=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1])
sch.rate<-ddply(sch.pois, .(ogs, type), summarize, rate=coef(summary(glm(count ~ offset(log(br)), family=poisson)))[1,1])

cons.rate<-merge(cons.rate, cons.event)
root.rate<-merge(root.rate, root.event)
sch.rate<-merge(sch.rate, sch.event)
