#TABLES

#Tables 1-3: just text

#Table 4: full immune list
write.table(mdom.dif[mdom.dif$immune | (mdom.dif$difexp > 0 & !is.na(mdom.dif$difexp)),c("gene", "hmm", "dmel", "dmel.imm", "dmel.sym", "log2FoldChange", "difexp")], file="full_immune_list.tsv", sep="\t", quote=F, row.names=F)

#Table 5: GO terms
#written out at the end of the GOanalysis script

