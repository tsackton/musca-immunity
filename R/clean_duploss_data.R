#construct ultrametric tree
require(ape)

#load RAxML tree into R
dipt.tree<-read.tree(file="../data/ultra.reroot.nwk")
spec.tree<-read.tree(file="../data/dipt.tree")
#compute ultrametric tree using chronos function in ape package

#make list of known dates -- this is not meant to be 100% accurate
#rather it is intended to loosely calibrate to MYA
#data for Drosophila from Tamura et al (2004) http://www.ncbi.nlm.nih.gov/pubmed/12949132
#data for musca from Wiegmann et al 2003 http://www.ncbi.nlm.nih.gov/pubmed/14668115
#data for mosquitos from  http://www.ncbi.nlm.nih.gov/pubmed/16473530
#data for deepest split (mos/dros) from http://rstb.royalsocietypublishing.org/content/363/1496/1435.long

dipt.cal<-data.frame(node=c(15,24,26,16,22,18,20), age.min=c(250,145,70,39,10,49,44), age.max=c(300,200,70,121,17.4,75,65), soft.bounds=F)

dipt.ultra<-chronos(dipt.tree, calibration=dipt.cal)
dipt.ultra.nocal<-chronos(dipt.tree)
dipt.cafe=dipt.ultra
dipt.cafe.nocal=dipt.ultra.nocal
dipt.cafe$edge.length = round(dipt.cafe$edge.length,0)
dipt.cafe.nocal$edge.length = round(dipt.cafe.nocal$edge.length,0)

#write ultrametric tree as newick
write.tree(dipt.ultra, file="../data/ultra.final.nwk")
write.tree(dipt.ultra.nocal, file="../data/ultra.final.nocal.nwk")
write.tree(dipt.cafe, file="../data/ultra.forcafe.nwk")
write.tree(dipt.cafe.nocal, file="../data/ultra.forcafe.nocal.nwk")

#done with ultrametric tree construction

#now make edge key - this is a bit ugly because there is no obvious way to automate it
edge.key=data.frame(treefix=c(names(mdom.all)[28:54]),ultra=c(3,4,2,13,12,15,11,14,10,9,8,7,6,5,17,16,1,26,25,24,23,22,21,20,19,18,NA))
ultra.edges<-data.frame(edge=seq(1,length(edge.key$ultra)-1,1), br=dipt.ultra$edge.length)
ultra.edges=merge(ultra.edges, edge.key, by.x="edge", by.y="ultra")
ultra.edges$nodeclass = "tip"
ultra.edges$nodeclass[grepl("^X", ultra.edges$treefix, perl=T)]="internal"
ultra.edges$nodeclass[ultra.edges$edge=="1"] = "root"
ultra.edges$nodeclass[ultra.edges$edge=="18"] = "root"
ultra.edges$branch = ultra.edges$edge
ultra.edges$branch[ultra.edges$nodeclass=="tip"] = as.character(ultra.edges$treefix[ultra.edges$nodeclass=="tip"])

#add clade information
ultra.edges$family="Diptera" #interal edges 1, 18, 2, 5 are generic
ultra.edges$family[grepl("^dro", ultra.edges$branch, perl=T)]="Drosophilidae" #Drosophilidae
ultra.edges$family[ultra.edges$branch == 15 | ultra.edges$branch == 10 | ultra.edges$branch == 8 | ultra.edges$branch == 6 | ultra.edges$branch == 11] = "Drosophilidae"
ultra.edges$family[ultra.edges$branch == "musDom"] = "Muscidae"
ultra.edges$family[ultra.edges$branch == "gloMor"] = "Glossinidae"
ultra.edges$family[ultra.edges$branch == 24 | ultra.edges$branch == 22 | ultra.edges$branch == 19 | ultra.edges$branch == "aedAeg" | ultra.edges$branch == "culQui"] = "Culicidae"
ultra.edges$family[grepl("^ano", ultra.edges$branch, perl=T)]="Culicidae" 

#add uncalibrated branch lengths to ultra.edges
ultra.edges$br.nocal = dipt.ultra.nocal$edge.length

#musca is just a shorthand for family=="Muscidae
ultra.edges$musca="other_dipt"
ultra.edges$musca[ultra.edges$branch=="musDom"]="musca"

#clean up count data for analysis and model fitting

#clean up for analysis
require(reshape2)

#look at three different data sets

#the first set is all groups inferrred to be present at all four major internal branches (edges 1, 18, 2, and 5 -- so in the ancestors of each family, roughly)
#define based on total counts
ogs.classes<-dup.ct[dup.ct$type=="total",]
ogs.classes$cons = F
ogs.classes$cons[ogs.classes$X1 != "NP" & ogs.classes$X6 != "NP" & ogs.classes$X7 != "NP" & ogs.classes$X8 != "NP" & ogs.classes$X2 != "NP"] = T

#second is just all genes present at the root, so this will include some mos/dros only genes or mos/mus only genes
ogs.classes$root = F
ogs.classes$root[ogs.classes$X1 != "NP"] = T

#final is the set of genes not present at the root, but present in the common ancestory of drosophila and musca (Schizophora)
ogs.classes$schiz = F
ogs.classes$schiz[ogs.classes$X1 == "NP" & ogs.classes$X6 != "NP" & ogs.classes$X8 != "NP" & ogs.classes$X7 != "NP"] = T

dup.cons<-droplevels(dup.ct[dup.ct$ogs %in% ogs.classes$ogs[ogs.classes$cons==T],])
dup.cons<-replace(dup.cons, dup.cons=="NP", 0)
dup.ct.cons<-melt(dup.cons, id.vars=c(1,2), variable.name="treefix", value.name="count")
dup.ct.cons=droplevels(dup.ct.cons[dup.ct.cons$treefix != "total" & dup.ct.cons$treefix!="X1",])

dup.root<-droplevels(dup.ct[dup.ct$ogs %in% ogs.classes$ogs[ogs.classes$root==T],])
dup.root<-replace(dup.root, dup.root=="NP", 0)
dup.ct.root<-melt(dup.root, id.vars=c(1,2), variable.name="treefix", value.name="count")
dup.ct.root=droplevels(dup.ct.root[dup.ct.root$treefix != "total" & dup.ct.root$treefix!="X1",])

dup.sch<-droplevels(dup.ct[dup.ct$ogs %in% ogs.classes$ogs[ogs.classes$schiz==T],])
dup.sch<-replace(dup.sch, dup.sch=="NP", 0)
dup.ct.sch<-melt(dup.sch, id.vars=c(1,2), variable.name="treefix", value.name="count")
dup.ct.sch=droplevels(dup.ct.sch[dup.ct.sch$treefix != "total" & dup.ct.sch$treefix!="X1" & dup.ct.sch$treefix != "X2" & dup.ct.sch$treefix != "X3" & dup.ct.sch$treefix != "X4" & dup.ct.sch$treefix != "X5" & dup.ct.sch$treefix != "aedAeg" & dup.ct.sch$treefix != "culQui" & dup.ct.sch$treefix != "anoDar" & dup.ct.sch$treefix != "anoGam" & dup.ct.sch$treefix != "anoSte",])

#dup.ct{cons,root,sch} contain the conserved set, the root set, and the schi only set

#now merge ultrametric tree info for each set
cons.pois<-merge(dup.ct.cons, ultra.edges, by.x="treefix", by.y="treefix")
cons.pois$count<-as.numeric(cons.pois$count)

root.pois<-merge(dup.ct.root, ultra.edges, by.x="treefix", by.y="treefix")
root.pois$count<-as.numeric(root.pois$count)

sch.pois<-merge(dup.ct.sch, ultra.edges, by.x="treefix", by.y="treefix")
sch.pois$count<-as.numeric(sch.pois$count)

#add event counts
require(plyr)

cons.event<-ddply(cons.pois[cons.pois$type=="total",], .(ogs), summarize, events=sum(count))
root.event<-ddply(root.pois[root.pois$type=="total",], .(ogs), summarize, events=sum(count))
sch.event<-ddply(sch.pois[root.pois$type=="total",], .(ogs), summarize, events=sum(count))

cons.pois<-merge(cons.pois, cons.event, by="ogs")
root.pois<-merge(root.pois, root.event, by="ogs")
sch.pois<-merge(sch.pois, sch.event, by="ogs")

cons.pois$subset="conserved"
root.pois$subset="root"
sch.pois$subset="sch"

all.pois<-rbind(cons.pois, root.pois, sch.pois)

#write cleaned output
write.table(all.pois, "../data/musca_counts_for_poisson.tsv", sep="\t", row.names=F, quote=F)

#now make the same set of outputs but with species counts for CAFE
#in this case, include full set (not filtered for events) and filtered set (events>=1), except for sch only which by definition has more than one event (loss in mosquitoes)

ogs.for.cafe<-ogs.ct[,c(1:15)]
ogs.for.cafe$desc="none"
ogs.for.cafe=ogs.for.cafe[,c(1,16,2:15)]

cons.cafe<-merge(ogs.for.cafe, cons.event, by="ogs")
root.cafe<-merge(ogs.for.cafe, root.event, by="ogs")
sch.cafe<-merge(ogs.for.cafe, sch.event, by="ogs")

write.table(cons.cafe[cons.cafe$events>0,c(1:16)], "../data/ogs_cafe_filt_cons.tsv", sep="\t", row.names=F, quote=F)
write.table(cons.cafe[,c(1:16)], "../data/ogs_cafe_cons.tsv", sep="\t", row.names=F, quote=F)
write.table(root.cafe[cons.cafe$events>0,c(1:16)], "../data/ogs_cafe_filt_root.tsv", sep="\t", row.names=F, quote=F)
write.table(root.cafe[,c(1:16)], "../data/ogs_cafe_root.tsv", sep="\t", row.names=F, quote=F)
write.table(sch.cafe[,c(1:16)], "../data/ogs_cafe_sch.tsv", sep="\t", row.names=F, quote=F)
