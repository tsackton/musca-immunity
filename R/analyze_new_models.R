#analyze PASA genes

pasa<-subset(mdom.all, class=="novel")
pasa<-unique(pasa[,c(1:6)])
known<-subset(mdom.all, class=="known")
known<-unique(known[,c(1:6)])

difexp.pasa<-difexp
difexp.pasa$class="known"
difexp.pasa$class[grep("PASA", rownames(difexp.pasa))]="novel"
boxplot(difexp.pasa$baseMean ~ difexp.pasa$class, outline=F)