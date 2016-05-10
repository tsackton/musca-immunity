#Load differential expression data and analyze with DESeq2
require(DESeq2)
require(gplots)
require(RColorBrewer)

#load rsem results
rsem.inf1.g <- read.table("../input_data/rsem/dmel_inf1.genes.results", header=T)
rsem.inf1.i <- read.table("../input_data/rsem/dmel_inf1.isoforms.results", header=T)
rsem.inf2.g <- read.table("../input_data/rsem/dmel_inf2.genes.results", header=T)
rsem.inf2.i <- read.table("../input_data/rsem/dmel_inf2.isoforms.results", header=T)
rsem.inf3.g <- read.table("../input_data/rsem/dmel_inf3.genes.results", header=T)
rsem.inf3.i <- read.table("../input_data/rsem/dmel_inf3.isoforms.results", header=T)
rsem.unf1.g <- read.table("../input_data/rsem/dmel_unf1.genes.results", header=T)
rsem.unf1.i <- read.table("../input_data/rsem/dmel_unf1.isoforms.results", header=T)
rsem.unf2.g <- read.table("../input_data/rsem/dmel_unf2.genes.results", header=T)
rsem.unf2.i <- read.table("../input_data/rsem/dmel_unf2.isoforms.results", header=T)
rsem.unf3.g <- read.table("../input_data/rsem/dmel_unf3.genes.results", header=T)
rsem.unf3.i <- read.table("../input_data/rsem/dmel_unf3.isoforms.results", header=T)

#set up data
rsem.gene.counts=cbind(rsem.unf1.g$exp, rsem.unf2.g$exp, rsem.unf3.g$exp, rsem.inf1.g$exp, rsem.inf2.g$exp, rsem.inf3.g$exp)
rsem.gene.counts<-round(rsem.gene.counts)
rownames(rsem.gene.counts)=rsem.unf1.g$gene_id
colnames(rsem.gene.counts)=c("Unf1", "Unf2", "Unf3", "Inf1", "Inf2", "Inf3")
conditions=data.frame(col=c("Unf1", "Unf2", "Unf3", "Inf1", "Inf2", "Inf3"), cond=c("unf", "unf", "unf", "inf", "inf", "inf"), stringsAsFactors=F)
conditions$cond=factor(conditions$cond, levels=c("unf", "inf"))

#Run DESeq2 with default options
des.rsem.input.dmel<-DESeqDataSetFromMatrix(countData=rsem.gene.counts, colData=conditions, design= ~ cond)
des.rsem.analysis.dmel<-DESeq(des.rsem.input.dmel)
des.rsem.results.dmel<-results(des.rsem.analysis.dmel)

#check results
summary(des.rsem.results.dmel)

#MA-plot, comparing shrunken vs unshrunked results
par(mfcol=c(1,2))
DESeq2::plotMA(des.rsem.results.dmel, ylim=c(-3,3), las=1, main="Shruken")
resMLE<-results(des.rsem.analysis.dmel, addMLE=TRUE)
resMLE.plot<-data.frame(baseMean=resMLE$baseMean, log2FoldChange=resMLE$lfcMLE, sig=resMLE$padj<0.1)
resMLE.plot$sig[is.na(resMLE.plot$sig)]<-FALSE
plotMA(resMLE.plot, ylim=c(-3,3), las=1, main="Raw")
par(mfcol=c(1,1))

#create regularized log transformation for heatmaps, etc
rld<-rlog(des.rsem.analysis.dmel)
distsRL<-as.matrix(dist(t(assay(rld))))
rownames(distsRL) <- colnames(distsRL) <- with(colData(des.rsem.analysis.dmel), cond)
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(distsRL, col=rev(hmcol), symm=T, trace="none", key=F)

#PCA
print(plotPCA(rld, intgroup=c("cond")))

#Dispersion plot
plotDispEsts(des.rsem.analysis.dmel)

#difexp has rsem analysis
difexp_dmel=as.data.frame(des.rsem.results.dmel)

#write table
write.table(difexp_dmel, file="../results/dmel.difexp.tsv", sep="\t", quote=F, row.names=T)


