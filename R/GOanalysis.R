#GO analysis

library(GOstats)
library(AnnotationForge)
library(GSEABase)
library(plyr)

#create a GO annotation in R from the blast2GO results from the genome paper

musca.go<-read.table("../input_data/GO/mdom_go_list", sep="\t", header=F)
names(musca.go)<-c("go_id", "evidence", "prot_id")

#replace protein ids with gene ids
prot.to.gene<-protkey
prot.to.gene$protein<-sub("\\.\\d+$", "", prot.to.gene$protein)
musca.go<-merge(musca.go, prot.to.gene, by.x="prot_id", by.y="protein", all.x=T, all.y=F)
musca.go<-musca.go[,c("go_id", "evidence", "gene")]
musca.go.all<-GOAllFrame(GOFrame(musca.go))
musca.gsc<-GeneSetCollection(musca.go.all, setType=GOCollection())

#define the background set as all genes not filtered in the difexp analysis, excluding those with an annotated immune function

mdom.imm.list<-unique(mdom.all[,c("gene", "hmm", "dmel.imm")])
mdom.imm.list<-mdom.imm.list[!is.na(mdom.imm.list$hmm) | !is.na(mdom.imm.list$dmel.imm),]

universe<-rownames(difexp[!is.na(difexp$padj),]) 
#universe<-universe[!(universe %in% mdom.imm.list$gene)]


#define upreg and downreg gene sets with 10% FDR

upreg<-rownames(difexp[!is.na(difexp$padj) & difexp$padj<0.05 & difexp$log2FoldChange > 0,])
#upreg<-upreg[!(upreg %in% mdom.imm.list$gene)]
downreg<-rownames(difexp[!is.na(difexp$padj) & difexp$padj<0.05 & difexp$log2FoldChange < 0,])
#downreg<-downreg[!(downreg %in% mdom.imm.list$gene)]
  
#define parameters for analysis

up.params<-GSEAGOHyperGParams(name="Musca upreg GSEA",
                              geneSetCollection=musca.gsc,
                              geneIds = upreg,
                              universeGeneIds = universe,
                              ontology = "BP",
                              pvalueCutoff = 0.05,
                              conditional = T,
                              testDirection = "over")

down.params<-GSEAGOHyperGParams(name="Musca downreg GSEA",
                              geneSetCollection=musca.gsc,
                              geneIds = downreg,
                              universeGeneIds = universe,
                              ontology = "BP",
                              pvalueCutoff = 0.05,
                              conditional = T,
                              testDirection = "over")

#perform tests
upreg.over.res<-hyperGTest(up.params)
downreg.over.res<-hyperGTest(down.params)

#verify similar result with limma-voom
library(edgeR)
library(limma)

md.dge<-DGEList(counts=counts(des.rsem.input, normalized=F))
md.dge<-calcNormFactors(md.dge)
md.design<-data.frame(unf=c(1,1,1,1,1,1), inf=c(0,0,0,1,1,1))
md.v<-voomWithQualityWeights(md.dge,md.design,plot=T,normalization="none")
md.fit<-lmFit(md.v,md.design)
md.fit<-eBayes(md.fit)

#construct indicies for limma-based gene set tests
musca.golist<-dlply(musca.go[,c("go_id", "gene")], .(go_id), .fun=function(x) unique(c(as.character(x$gene))))
limma.go<-ids2indices(musca.golist, rownames(md.v))

#run the ROMER test from limma to get GO categories with evidence for upregulation, downregulation or regulation overall ('mixed')
romer.gse<-romer(y=md.v, index=limma.go, design=md.design, contrast=2)


