#test false positive rate from overdispersion

sim.ct.tot <- subset(sim.root.pois, type == "total")
sim.ct.sep <- subset(sim.root.pois, type != "total")

imm.broad.count<-table(unique(subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root", select=c("ogs", "immune.broad")))[,2])
imm.narrow.count<-table(unique(subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root", select=c("ogs", "immune.narrow")))[,2])
imm.rec.count<-table(unique(subset(all.pois, subset=="conserved" & type=="total" & nodeclass != "root", select=c("ogs", "dmel.imm")))[,2])

sim.ogs<-data.frame(ogs=c(unique(sim.ct.tot$ogs)))
narrow.ct<-round(length(sim.ogs$ogs)*(imm.narrow.count[2]/sum(imm.narrow.count)),0)
rec.ct<-round(length(sim.ogs$ogs)*(imm.rec.count[2]/sum(imm.rec.count)),0)
broad.ct<-round(length(sim.ogs$ogs)*(imm.broad.count[2]/sum(imm.broad.count)),0)
sim.ogs$rate = as.numeric(sub("sim\\.\\d+", "", sim.ogs$ogs, perl=T))

sim.results.imm = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
sim.results.rec = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
for (i in seq(1,200)) {
  sim.ogs$immune.narrow = FALSE
  sim.ogs$rec = FALSE
  sim.ogs$immune.narrow[sample(nrow(sim.ogs), narrow.ct)] = TRUE
  sim.ogs$rec[sample(nrow(sim.ogs), rec.ct)] = TRUE
  sim.ct.test = merge(sim.ct.sep, sim.ogs)
  sim.test<-summary(glm(count ~ type*musca*immune.narrow, offset=log(br), family="poisson", data=sim.ct.test))
  rec.test<-summary(glm(count ~ type*musca*rec, offset=log(br), family="poisson", data=sim.ct.test))
  sim.results.imm[i,2]=sim.test$coef[4,4]
  sim.results.imm[i,3]=sim.test$coef[7,4]
  sim.results.imm[i,4]=sim.test$coef[4,1]
  sim.results.imm[i,5]=sim.test$coef[7,1]
  sim.results.rec[i,2]=rec.test$coef[4,4]
  sim.results.rec[i,3]=rec.test$coef[7,4]
  sim.results.rec[i,4]=rec.test$coef[4,1]
  sim.results.rec[i,5]=rec.test$coef[7,1]
}

sim.onerate.results.imm = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
sim.onerate.results.rec = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
sim.ogs.onerate = subset(sim.ogs, rate <= 0.002)
narrow.ct<-round(length(sim.ogs.onerate$ogs)*(imm.narrow.count[2]/sum(imm.narrow.count)),0)
rec.ct<-round(length(sim.ogs.onerate$ogs)*(imm.rec.count[2]/sum(imm.rec.count)),0)
broad.ct<-round(length(sim.ogs.onerate$ogs)*(imm.broad.count[2]/sum(imm.broad.count)),0)

for (i in seq(1,200)) {
  sim.ogs.onerate$immune.narrow = FALSE
  sim.ogs.onerate$rec = FALSE
  sim.ogs.onerate$immune.narrow[sample(nrow(sim.ogs.onerate), narrow.ct)] = TRUE
  sim.ogs.onerate$rec[sample(nrow(sim.ogs.onerate), rec.ct)] = TRUE
  sim.ct.test = merge(sim.ct.sep, sim.ogs.onerate)
  sim.test<-summary(glm(count ~ type*musca*immune.narrow, offset=log(br), family="poisson", data=sim.ct.test))
  rec.test<-summary(glm(count ~ type*musca*rec, offset=log(br), family="poisson", data=sim.ct.test))
  sim.onerate.results.imm[i,2]=sim.test$coef[4,4]
  sim.onerate.results.imm[i,3]=sim.test$coef[7,4]
  sim.onerate.results.imm[i,4]=sim.test$coef[4,1]
  sim.onerate.results.imm[i,5]=sim.test$coef[7,1]
  sim.onerate.results.rec[i,2]=rec.test$coef[4,4]
  sim.onerate.results.rec[i,3]=rec.test$coef[7,4]
  sim.onerate.results.rec[i,4]=rec.test$coef[4,1]
  sim.onerate.results.rec[i,5]=rec.test$coef[7,1]
}


##try quasipoisson as well

sim.results.qp.imm = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
sim.results.qp.rec = data.frame(iter=seq(1,200), immune.pval=numeric(200), interaction.pval=numeric(200), immune.coef = numeric(200), interaction.coef = numeric(200))
for (i in seq(1,200)) {
  sim.ogs$immune.narrow = FALSE
  sim.ogs$rec = FALSE
  sim.ogs$immune.narrow[sample(nrow(sim.ogs), narrow.ct)] = TRUE
  sim.ogs$rec[sample(nrow(sim.ogs), rec.ct)] = TRUE
  sim.ct.test = merge(sim.ct.sep, sim.ogs)
  sim.test<-summary(glm(count ~ type*musca*immune.narrow, offset=log(br), family="quasipoisson", data=sim.ct.test))
  rec.test<-summary(glm(count ~ type*musca*rec, offset=log(br), family="quasipoisson", data=sim.ct.test))
  sim.results.qp.imm[i,2]=sim.test$coef[4,4]
  sim.results.qp.imm[i,3]=sim.test$coef[7,4]
  sim.results.qp.imm[i,4]=sim.test$coef[4,1]
  sim.results.qp.imm[i,5]=sim.test$coef[7,1]
  sim.results.qp.rec[i,2]=rec.test$coef[4,4]
  sim.results.qp.rec[i,3]=rec.test$coef[7,4]
  sim.results.qp.rec[i,4]=rec.test$coef[4,1]
  sim.results.qp.rec[i,5]=rec.test$coef[7,1]
}

sim.results.qp.imm$type = "quasipoission-imm"
sim.results.qp.rec$type = "quasipoission-rec"
sim.onerate.results.imm$type = "poisson-lowratevar-imm"
sim.onerate.results.rec$type = "poisson-lowratevar-rec"
sim.results$type = "poisson-imm"
sim.results.rec$type = "poisson-rec"

all.sim.results <- rbind(sim.results.qp.imm, sim.results.qp.rec, sim.onerate.results.imm, sim.onerate.results.rec, sim.results, sim.results.rec)
write.table(all.sim.results, file="model_check_results.200iter.out", sep="\t", row.names=F, quote=F)
