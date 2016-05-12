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


#define upreg and downreg gene sets with 5% FDR

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

#multiple test correction
upreg.df<-summary(upreg.over.res, pvalue=1, categorySize=5)
upreg.df$padj<-p.adjust(upreg.df$Pvalue, method="holm")
downreg.df<-summary(downreg.over.res, pvalue=1, categorySize=5)
downreg.df$padj<-p.adjust(downreg.df$Pvalue, method="holm")

#molecular function
up.params.m<-GSEAGOHyperGParams(name="Musca upreg GSEA",
                              geneSetCollection=musca.gsc,
                              geneIds = upreg,
                              universeGeneIds = universe,
                              ontology = "MF",
                              pvalueCutoff = 0.05,
                              conditional = T,
                              testDirection = "over")

down.params.m<-GSEAGOHyperGParams(name="Musca downreg GSEA",
                                geneSetCollection=musca.gsc,
                                geneIds = downreg,
                                universeGeneIds = universe,
                                ontology = "MF",
                                pvalueCutoff = 0.05,
                                conditional = T,
                                testDirection = "over")

#perform tests
upreg.over.res.m<-hyperGTest(up.params.m)
downreg.over.res.m<-hyperGTest(down.params.m)

#multiple test correction
upreg.df.m<-summary(upreg.over.res.m, pvalue=1, categorySize=5)
upreg.df.m$padj<-p.adjust(upreg.df.m$Pvalue, method="holm")
downreg.df.m<-summary(downreg.over.res.m, pvalue=1, categorySize=5)
downreg.df.m$padj<-p.adjust(downreg.df.m$Pvalue, method="holm")


#write results out
write.table(upreg.df.m, file="../results/GO_MF_upreg.tsv", sep="\t", quote=F, row.names=F)
write.table(upreg.df, file="../results/GO_BP_upreg.tsv", sep="\t", quote=F, row.names=F)
write.table(downreg.df.m, file="../results/GO_MF_downreg.tsv", sep="\t", quote=F, row.names=F)
write.table(downreg.df, file="../results/GO_BP_downreg.tsv", sep="\t", quote=F, row.names=F)

##DMEL ANALYSIS##

#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Dm.eg.db")
library(org.Dm.eg.db)

#define the background set as all genes not filtered in the difexp analysis, excluding those with an annotated immune function
universe<-rownames(difexp_dmel[!is.na(difexp_dmel$padj),]) 

#define upreg and downreg gene sets with 5% FDR

upreg<-rownames(difexp_dmel[!is.na(difexp_dmel$padj) & difexp_dmel$padj<0.05 & difexp_dmel$log2FoldChange > 0,])
downreg<-rownames(difexp_dmel[!is.na(difexp_dmel$padj) & difexp_dmel$padj<0.05 & difexp_dmel$log2FoldChange < 0,])

#define parameters for analysis

up.params<-new("GOHyperGParams",
                              geneIds = unlist(mget(upreg, org.Dm.egFLYBASE2EG, ifnotfound=NA))[!is.na(unlist(mget(upreg, org.Dm.egFLYBASE2EG, ifnotfound=NA)))],
                              universeGeneIds = unlist(mget(universe, org.Dm.egFLYBASE2EG, ifnotfound=NA))[!is.na(unlist(mget(universe, org.Dm.egFLYBASE2EG, ifnotfound=NA)))],
                              ontology = "BP",
                              pvalueCutoff = 0.05,
                              conditional = T,
                              testDirection = "over",
                              annotation = "org.Dm.eg.db")

down.params<-new("GOHyperGParams",
               geneIds = unlist(mget(downreg, org.Dm.egFLYBASE2EG, ifnotfound=NA))[!is.na(unlist(mget(downreg, org.Dm.egFLYBASE2EG, ifnotfound=NA)))],
               universeGeneIds = unlist(mget(universe, org.Dm.egFLYBASE2EG, ifnotfound=NA))[!is.na(unlist(mget(universe, org.Dm.egFLYBASE2EG, ifnotfound=NA)))],
               ontology = "BP",
               pvalueCutoff = 0.05,
               conditional = T,
               testDirection = "over",
               annotation = "org.Dm.eg.db")

#perform tests
upreg.over.res.dmel<-hyperGTest(up.params)
downreg.over.res.dmel<-hyperGTest(down.params)

#multiple test correction
upreg.df.dmel<-summary(upreg.over.res.dmel, pvalue=1, categorySize=5)
upreg.df.dmel$padj<-p.adjust(upreg.df.dmel$Pvalue, method="holm")
downreg.df.dmel<-summary(downreg.over.res.dmel, pvalue=1, categorySize=5)
downreg.df.dmel$padj<-p.adjust(downreg.df.dmel$Pvalue, method="holm")


#write results out
write.table(upreg.df.dmel, file="../results/GO_BP_upreg_dmel.tsv", sep="\t", quote=F, row.names=F)
write.table(downreg.df.dmel, file="../results/GO_BP_downreg_dmel.tsv", sep="\t", quote=F, row.names=F)

