#construct ultrametric tree
require(ape)

#load RAxML tree into R
dipt.tree<-read.tree(file="./data/ultra.reroot.nwk")
spec.tree<-read.tree(file="./data/dipt.tree")
#compute ultrametric tree (primitively)
dipt.ultra<-chronos(dipt.tree)

#write ultrametric tree as newick
write.tree(dipt.ultra, file="./data/ultra.final.nwk")

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

#write cleaned output
write.table(cons.pois, "./data/musca_counts_for_poisson_cons.tsv", sep="\t", row.names=F, quote=F)
write.table(root.pois, "./data/musca_counts_for_poisson_root.tsv", sep="\t", row.names=F, quote=F)
write.table(sch.pois, "./data/musca_counts_for_poisson_sch.tsv", sep="\t", row.names=F, quote=F)
