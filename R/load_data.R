#Final Version. April 2015

#load data into R and check for consistency#
#This code loads transcript and protein ids, protein lengthts, phylostratigraphic data, immune hmm calls, signal peptide inference, and orthogroup data from the OMA/Treefix pipeline

#get transcripts
isokey<-read.table("../data/mdom.final.rnakey", header=F)
names(isokey) = c("gene", "transcript")

#get proteins
protkey<-read.table("../data/mdom.final.protkey", header=F)
names(protkey) = c("gene", "protein")

#get longest isoform
mdom.longest=read.table("../data/longest.musDom.txt", header=F, sep="|")
names(mdom.longest) = c("sp", "gene", "isoform")

#load strata
mdom.strat<-read.table("../data/mdom.strata", header=T)

#get signal peptide calls
sigpep<-read.table("../data/mdom.sigpep")
sigpep=sigpep[,c(1,10)]
names(sigpep) = c("id", "sigp")

#load hmm
mdom.imm<-read.table("../data/musca.immune.hmm", header=F)
names(mdom.imm) = c("prot", "hmm", "eval")
mdom.imm = mdom.imm[order(mdom.imm$eval),]
mdom.imm = mdom.imm[!duplicated(mdom.imm$prot),]

#dmel immune annotation
dmel.imm<-read.table("../data/dmel-immune.txt", sep="\t", header=T)

#basic ogs stats
ogs.ct<-read.table("../data/ogs_stats.txt", header=T, stringsAsFactors=F)
mdom.ogs.key<-read.table("../data/mdom_to_ogs.txt", header=T)
mdom.ogs.key$mdomid=sub("_XP_\\d+.*$", "", mdom.ogs.key$mdomid, perl=T)
mdom.ogs.key$mdomid=sub("_PASA\\d+_model.*", "", mdom.ogs.key$mdomid, perl=T)
dmel.ogs.key<-read.table("../data/dmel_to_ogs.txt", header=T)
dmel.ogs.key$dmelid=sub("_FBpp\\d+", "", dmel.ogs.key$dmelid, perl=T)
ogs.ct$total = apply(ogs.ct[,2:15], 1, sum)
ogs.ct$max = apply(ogs.ct[,2:15], 1, max)

#merge strata and proteins
mdom.prot<-merge(protkey, mdom.strat, by.x="protein", by.y="id", all.x=T, all.y=F)

#merge with ogs key
mdom.genes<-merge(mdom.prot, mdom.ogs.key, by.x="gene", by.y="mdomid", all=T)
mdom.genes = mdom.genes[!is.na(mdom.genes$gene),]

#remove duplicates by keeping result for longest protein only
mdom.genes=merge(mdom.longest, mdom.genes, all.x=T, all.y=F, by.x="isoform", by.y="protein")
mdom.genes = subset(mdom.genes, select=-c(gene.y, sp))
names(mdom.genes)[2]="gene"

#mdom.genes now has strata and OGS info for all protein-coding genes in the genome

#read treefix data
dup.ct<-read.table("../data/ogs_dup_ct.txt", header=T, stringsAsFactors=F)
dmel.mdom.hom<-read.table("../data/dmel_mdom_homology.txt", header=T)
dmel.mdom.hom$mdom=sub("_XP_\\d+.*$", "", dmel.mdom.hom$mdom, perl=T)
dmel.mdom.hom$mdom=sub("_PASA\\d+_model.*", "", dmel.mdom.hom$mdom, perl=T)
dmel.mdom.hom$dmel=gsub("_FBpp\\d+", "", dmel.mdom.hom$dmel, perl=T)

#merge dmel immune annotation
parse.class <- function(x) {
  x=as.character(x[3])
  com<-unlist(strsplit(x, split=c(","), fixed=T))
  class<-paste(as.character(unique(sort(dmel.imm[dmel.imm$fbgn %in% com, "class"]))), collapse=";")
  if(class==""){
    class<-NA
  }
 return(as.character(class))
}

parse.sym <- function(x) {
  x=as.character(x[3])
  com<-unlist(strsplit(x, split=c(","), fixed=T))
  class<-paste(as.character(unique(sort(dmel.imm[dmel.imm$fbgn %in% com, "sym"]))), collapse=";")
  if(class==""){
    class<-NA
  }
  return(as.character(class))
}

dmel.mdom.hom$dmel.imm<-apply(dmel.mdom.hom, 1, parse.class)
dmel.mdom.hom$dmel.sym<-apply(dmel.mdom.hom, 1, parse.sym)

#merge everything
mdom.all<-merge(mdom.genes, mdom.imm, all.x=T, all.y=F, by.x="isoform", by.y="prot")
mdom.all<-merge(mdom.all, sigpep, all.x=T, all.y=F, by.x="isoform", by.y="id")
mdom.all = merge(mdom.all, dmel.mdom.hom, all.x=T, all.y=F, by.x="gene", by.y="mdom")
mdom.all = merge(mdom.all, dup.ct, by.x="ogsid", by.y="ogs", all.x=T, all.y=T)

##mdom.all has all information for protein-coding genes##
##write out mdom.all##

write.table(mdom.all, file="../data/mdom.info.tsv", sep="\t", quote=F, row.names=F)

#get list of genes to use for ultrametric tree construction / species tree estimation
conserved.trees<-ogs.ct[ogs.ct$total==14 & ogs.ct$max==1 & !is.na(ogs.ct$total),]
#write out conserved tree IDs for processing
write.table(file="../data/ogs_for_ultrametric.txt", conserved.trees[,1], quote=F, sep="\t", row.names=F, col.names=F)

#get strata ages
strata.age<-read.table('../data/strata_age.txt', header=F)
names(strata.age)=c("strata", "min.age")

#get protein length
protlen<-read.table("../data/mdom.final.proteins.fa.lengths", header=F)
names(protlen)=c("isoform", "length")

#run difexp analysis,saving diagnostic plots to a PDF file
pdf(file="difexp_plots.pdf")
source("difexp.R")
dev.off()

##END##
