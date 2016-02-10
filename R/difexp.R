#Load differential expression data and analyze with DESeq2
require(DESeq2)
require(gplots)
require(RColorBrewer)

#load rsem results
rsem.inf1.g <- read.table("../input_data/rsem/fem_inf1.genes.results", header=T)
rsem.inf1.i <- read.table("../input_data/rsem/fem_inf1.isoforms.results", header=T)
rsem.inf2.g <- read.table("../input_data/rsem/fem_inf2.genes.results", header=T)
rsem.inf2.i <- read.table("../input_data/rsem/fem_inf2.isoforms.results", header=T)
rsem.inf3.g <- read.table("../input_data/rsem/fem_inf3.genes.results", header=T)
rsem.inf3.i <- read.table("../input_data/rsem/fem_inf3.isoforms.results", header=T)
rsem.unf1.g <- read.table("../input_data/rsem/fem_st1.genes.results", header=T)
rsem.unf1.i <- read.table("../input_data/rsem/fem_st1.isoforms.results", header=T)
rsem.unf2.g <- read.table("../input_data/rsem/fem_st2.genes.results", header=T)
rsem.unf2.i <- read.table("../input_data/rsem/fem_st2.isoforms.results", header=T)
rsem.unf3.g <- read.table("../input_data/rsem/fem_st3.genes.results", header=T)
rsem.unf3.i <- read.table("../input_data/rsem/fem_st3.isoforms.results", header=T)

#set up data
rsem.gene.counts=cbind(rsem.unf1.g$exp, rsem.unf2.g$exp, rsem.unf3.g$exp, rsem.inf1.g$exp, rsem.inf2.g$exp, rsem.inf3.g$exp)
rsem.gene.counts<-round(rsem.gene.counts)
rownames(rsem.gene.counts)=rsem.unf1.g$gene_id
colnames(rsem.gene.counts)=c("Unf1", "Unf2", "Unf3", "Inf1", "Inf2", "Inf3")
conditions=data.frame(col=c("Unf1", "Unf2", "Unf3", "Inf1", "Inf2", "Inf3"), cond=c("unf", "unf", "unf", "inf", "inf", "inf"), stringsAsFactors=F)
conditions$cond=factor(conditions$cond, levels=c("unf", "inf"))

#Run DESeq2 with default options
des.rsem.input<-DESeqDataSetFromMatrix(countData=rsem.gene.counts, colData=conditions, design= ~ cond)
des.rsem.analysis<-DESeq(des.rsem.input)
des.rsem.results<-results(des.rsem.analysis)

#check results
summary(des.rsem.results)

#MA-plot, comparing shrunken vs unshrunked results
par(mfcol=c(1,2))
plotMA(des.rsem.results, ylim=c(-3,3), las=1, main="Shruken")
resMLE<-results(des.rsem.analysis, addMLE=TRUE)
resMLE.plot<-data.frame(baseMean=resMLE$baseMean, log2FoldChange=resMLE$lfcMLE, sig=resMLE$padj<0.1)
resMLE.plot$sig[is.na(resMLE.plot$sig)]<-FALSE
plotMA(resMLE.plot, ylim=c(-3,3), las=1, main="Raw")
par(mfcol=c(1,1))

#create regularized log transformation for heatmaps, etc
rld<-rlog(des.rsem.analysis)
distsRL<-as.matrix(dist(t(assay(rld))))
rownames(distsRL) <- colnames(distsRL) <- with(colData(des.rsem.analysis), cond)
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(distsRL, col=rev(hmcol), symm=T, trace="none", key=F)

#PCA
print(plotPCA(rld, intgroup=c("cond")))

#Dispersion plot
plotDispEsts(des.rsem.analysis)

#difexp has rsem analysis
difexp=as.data.frame(des.rsem.results)

#write table
write.table(difexp, file="../results/mdom.difexp.tsv", sep="\t", quote=F, row.names=T)


