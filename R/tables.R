#TABLES

#Tables 1-3: just text

#Table 4: full immune list
write.table(mdom.dif[mdom.dif$immune | (mdom.dif$difexp > 0 & !is.na(mdom.dif$difexp)),c("gene", "hmm", "dmel", "dmel.imm", "dmel.sym", "log2FoldChange", "difexp")], file="full_immune_list.tsv", sep="\t", quote=F, row.names=F)

#Table 5: GO terms
#written out at the end of the GOanalysis script

#Table 6 - HMM results
#need to clean up to get one per gene
hmm.results.pergene<-hmm.results
hmm.results.pergene$gene<-sub("-\\w+", "", hmm.results.pergene$prot, perl=T)
dros.protkey<-read.table("../input_data/annotations//fbgn_fbtr_fbpp_fb_2016_02.tsv.gz", sep="\t", fill=T)
dros.protkey<-dros.protkey[,c(1,3)]
names(dros.protkey)<-c("gene", "protein")
full.protkey<-rbind(protkey, dros.protkey)
hmm.results.pergene<-merge(hmm.results.pergene, full.protkey, by.x="prot", by.y="protein", all.x=T, all.y=F, suffixes=c(".1", ".2"))
hmm.results.pergene$gene = ifelse(is.na(hmm.results.pergene$gene.2), hmm.results.pergene$gene.1, hmm.results.pergene$gene.2)
hmm.results.pergene<-hmm.results.pergene[!duplicated(hmm.results.pergene$gene),]
hmm.results.table<-as.data.frame(table(hmm.results.pergene$hmm, hmm.results.pergene$species))
hmm.results.table$Var1 = factor(hmm.results.table$Var1, levels=c("ATT","DEF","DIPT","CEC","LYS","TPX","PPO","GPX","HPX","TSF","NIM","PGRP","BGBP","TEP","CTL","FREP","GALE","IGSF","MD2L","SRCA","SRCB","SRCC","NFKB","SPRN","TLL","CLIPA","CLIPB","CLIPC","CLIPD", "CLIPE"), ordered=T)
hmm.results.table<-hmm.results.table[order(hmm.results.table$Var1),]
hmm.results.table<-hmm.results.table[order(hmm.results.table$Var2),]

write.table(hmm.results.table, "hmm_results_allspecies.tsv", sep="\t", quote=F, row.names=F)

#table 7 -text

#table 8, 9, 10: in duploss_analysis.R script

