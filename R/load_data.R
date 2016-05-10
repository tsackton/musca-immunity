#Final Version. April 2016

#load data into R and check for consistency#
#This code loads transcript and protein ids, protein lengthts, phylostratigraphic data, immune hmm calls, signal peptide inference, and orthogroup data from the OMA/Treefix pipeline

#set working directory
setwd("~/Projects/immunity/musca/musca-immunity/R/")

#get transcripts
isokey<-read.table("../input_data/annotations/mdom.final.rnakey", header=F)
names(isokey) = c("gene", "transcript")

#get proteins
protkey<-read.table("../input_data/annotations/mdom.final.protkey", header=F)
names(protkey) = c("gene", "protein")
protkey$protein = sub("\\.\\d+", "", protkey$protein, perl=T)

#added: dmel protein key
dmel.protkey<-read.table("../input_data/annotations/dmel.protkey", header=F)
names(dmel.protkey) = c("gene", "protein")

#get protein length
protlen<-read.table("../input_data/annotations/mdom.final.proteins.fa.lengths", header=F)
names(protlen)=c("isoform", "length")
protlen$isoform = sub("\\.\\d+", "", protlen$isoform, perl=T)


#get longest isoform
mdom.longest=merge(protkey, protlen, by.x="protein", by.y="isoform", all.x=T)
mdom.longest=mdom.longest[order(mdom.longest$gene, mdom.longest$length, decreasing=T),]
mdom.longest=mdom.longest[!duplicated(mdom.longest$gene),]

#load strata
mdom.strat<-read.table("../input_data/annotations/mdom_strata.txt", header=T)
mdom.strat$id =  sub("\\.\\d+", "", mdom.strat$id, perl=T)


#get strata ages
strata.age<-read.table('../supplemental_methods/strata/strata_age.txt', header=F)
names(strata.age)=c("strata", "min.age")


#get signal peptide calls
sigpep<-read.table("../input_data/annotations/mdom.sigpep")
sigpep=sigpep[,c(1,10)]
names(sigpep) = c("id", "sigp")
sigpep$id = sub("\\.\\d+", "", sigpep$id, perl=T)

#load hmm
#edited April 2016 to load all species HMM annotations
hmm.results<-read.table("../supplemental_methods/hmm/parsed_hmm_results.tsv", header=F)
names(hmm.results) = c("species", "prot", "hmm", "eval")
hmm.results = hmm.results[order(hmm.results$eval),]
hmm.results = hmm.results[!duplicated(hmm.results$prot),]
mdom.imm = hmm.results[hmm.results$species=="mdom.final.proteins", c("prot", "hmm", "eval")]
dmel.hmm.imm = hmm.results[hmm.results$species=="dmel-all-translation-r6.02", c("prot", "hmm", "eval")]
dmel.hmm.imm = merge(dmel.hmm.imm, dmel.protkey, by.x="prot", by.y="protein", all.x=T, all.y=F)
dmel.hmm.imm = dmel.hmm.imm[order(dmel.hmm.imm$eval),]
dmel.hmm.imm = dmel.hmm.imm[!duplicated(dmel.hmm.imm$gene),]

#dmel immune annotation
dmel.imm<-read.table("../input_data/annotations/dmel-immune.txt", sep="\t", header=T)

#basic ogs stats
ogs.ct<-read.table("../input_data/homology/ogs_stats.txt", header=T, stringsAsFactors=F)
mdom.ogs.key<-read.table("../input_data/homology/mdom_to_ogs.txt", header=T)
mdom.ogs.key$mdomid=sub("_XP_\\d+.*$", "", mdom.ogs.key$mdomid, perl=T)
mdom.ogs.key$mdomid=sub("_PASA\\d+_model.*", "", mdom.ogs.key$mdomid, perl=T)
dmel.ogs.key<-read.table("../input_data/homology/dmel_to_ogs.txt", header=T)
dmel.ogs.key$dmelid=sub("_FBpp\\d+", "", dmel.ogs.key$dmelid, perl=T)
ogs.ct$total = apply(ogs.ct[,2:15], 1, sum)
ogs.ct$max = apply(ogs.ct[,2:15], 1, max)

#merge strata and proteins
mdom.prot<-merge(protkey, mdom.strat, by.x="protein", by.y="id", all.x=T, all.y=F)

#merge with ogs key
mdom.genes<-merge(mdom.prot, mdom.ogs.key, by.x="gene", by.y="mdomid", all=T)
mdom.genes = mdom.genes[!is.na(mdom.genes$gene),]

#remove duplicates by keeping result for longest protein only
mdom.genes=merge(mdom.longest, mdom.genes, all.x=T, all.y=F, by.x="protein", by.y="protein")
mdom.genes = subset(mdom.genes, select=-c(gene.y))
names(mdom.genes)[2]="gene"

#mdom.genes now has strata and OGS info for all protein-coding genes in the genome

#read treefix data
dup.ct<-read.table("../input_data/homology/ogs_dup_ct.txt", header=T, stringsAsFactors=F)
dmel.mdom.hom<-read.table("../input_data/homology/dmel_mdom_homology.txt", header=T)
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
mdom.all<-merge(mdom.genes, mdom.imm, all.x=T, all.y=F, by.x="protein", by.y="prot")
mdom.all<-merge(mdom.all, sigpep, all.x=T, all.y=F, by.x="protein", by.y="id")
mdom.all = merge(mdom.all, dmel.mdom.hom, all.x=T, all.y=F, by.x="gene", by.y="mdom")
mdom.all = merge(mdom.all, dup.ct, by.x="ogsid", by.y="ogs", all.x=T, all.y=T)

##mdom.all has all information for protein-coding genes##
##write out mdom.all##

write.table(mdom.all, file="../results/mdom.info.tsv", sep="\t", quote=F, row.names=F)

#get list of genes to use for ultrametric tree construction / species tree estimation
conserved.trees<-ogs.ct[ogs.ct$total==14 & ogs.ct$max==1 & !is.na(ogs.ct$total),]
#write out conserved tree IDs for processing
#write.table(file="../results/ogs_for_ultrametric.txt", conserved.trees[,1], quote=F, sep="\t", row.names=F, col.names=F)


#run difexp analysis,saving diagnostic plots to a PDF file
pdf(file="difexp_plots_musca.pdf")
source("difexp.R")
dev.off()

#run difexp analysis on dmel, saving diagnostic plots
pdf(file="difexp_plots_dmel.pdf")
source("difexp_dmel.R")
dev.off()

##END##
