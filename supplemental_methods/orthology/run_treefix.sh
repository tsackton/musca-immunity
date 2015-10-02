#!/bin/bash

#define variables
OG=$1
RUN=$2
TREEFIXHOME=/n/home12/tsackton/regal/musca/treefix
PYTHONPATH=$PYTHONPATH:/n/home12/tsackton/sw/lib/python2.7/site-packages

#remove all gap columns from aligned fasta back to fasta
trimal -in $TREEFIXHOME/$RUN/inAlign/$OG.aligned -out $OG.in.fa -fasta -noallgaps

#copy tree
cp $TREEFIXHOME/$RUN/inTree/$OG.nwk $OG.in.nwk

#fix bars
perl -p -i -e 's/\|/_/g;' $OG.in.fa
perl -p -i -e 's/U/X/g;' $OG.in.fa #replace selenium with X
perl -p -i -e 's\./-/g;' $OG.in.fa #replace stop with gap
perl -p -i -e 's/-(\d+)/\.$1/g;' $OG.in.fa #fix Musca deflines
perl -p -i -e 's/\|/_/g;' $OG.in.nwk

#get protein model
OGMOD=${OG%%.*}
PROTMAT=$(grep $OGMOD $TREEFIXHOME/ogs_to_model.key | cut -f2,2)
if [ -z "$PROTMAT" ]
then
	PROTMAT="JTT"
fi
MODEL=PROTGAMMA$PROTMAT

#run treefix
treefix -s $TREEFIXHOME/EstimatedSpeciesTree.nwk -S $TREEFIXHOME/dipt.smap -A .in.fa -o .in.nwk -n .out.nwk -e "-m $MODEL" --niter=1000 -V1 -l - $OG.in.nwk
tree-annotate -s $TREEFIXHOME/EstimatedSpeciesTree.nwk -S $TREEFIXHOME/dipt.smap -T .out.nwk $OG.out.nwk
