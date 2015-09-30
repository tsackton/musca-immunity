#!/bin/bash

LAUNCHDIR=$(pwd)
PASADIR=/home/tim/bin/PASA
SPEC=$1

#move to correct dir	
cd $SPEC	

#PASA cmds
TIME=$(date)
echo "[$TIME] started" > $LAUNCHDIR/$SPEC.timing

#make config files
cp $LAUNCHDIR/alignAssembly.config .
cp $LAUNCHDIR/annotCompare.config .

sed -i -e "s/REPLACETHIS/${SPEC}_new_db/" alignAssembly.config
sed -i -e "s/REPLACETHIS/${SPEC}_new_db/" annotCompare.config

#cleanup old runs
$PASADIR/scripts/drop_mysql_db_if_exists.dbi -c alignAssembly.config

#seqclean
$PASADIR/seqclean/seqclean/seqclean Trinity.$SPEC.fa
TIME=$(date)
echo "[$TIME] finished seqclean" >> $LAUNCHDIR/$SPEC.timing

#get de novo accessions
$PASADIR/misc_utilities/accession_extractor.pl < Trinity.$SPEC.fa > tdn.accs

#map transcripts
$PASADIR/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $SPEC.fa --ALIGNERS gmap -N 2 --TRANSDECODER -T -u Trinity.$SPEC.fa -t Trinity.$SPEC.fa.clean --TDN tdn.accs -L --annots_gff3 $SPEC.coding.gff --gene_overlap 50.0 --CPU 4
TIME=$(date)
echo "[$TIME] finished transcript mapping" >> $LAUNCHDIR/$SPEC.timing

#compare transcripts to existing annotation:

#annotation compare
$PASADIR/scripts/Launch_PASA_pipeline.pl -c annotCompare.config -A -L -g $SPEC.fa -t Trinity.$SPEC.fa.clean --annots_gff3 $SPEC.coding.gff --CPU 4
TIME=$(date)
echo "[$TIME] finished transcript compare round 1" >> $LAUNCHDIR/$SPEC.timing

#get new gff
NEWGFF=$(ls -t *.gff3 | head -n 1)

#rename ids and check for dups
perl -p -i.bak1 -e "s/novel_gene/PASA1_gene/g" $NEWGFF
perl -p -i.bak2 -e "s/(temp)|(novel)_model_/PASA1_model_/g" $NEWGFF
../fix_dup_ids.pl $NEWGFF
cp $NEWGFF $NEWGFF.bak3
mv $NEWGFF.updated $NEWGFF

#second round annotation
$PASADIR/scripts/Launch_PASA_pipeline.pl -c annotCompare.config -A -L -g $SPEC.fa -t Trinity.$SPEC.fa.clean --annots_gff3 $NEWGFF --CPU 4
TIME=$(date)
echo "[$TIME] finished transcript compare round 2" >> $LAUNCHDIR/$SPEC.timing

#get final gff
FINALGFF=$(ls -t *.gff3 | head -n 1)

#rename ids and check for dups
perl -p -i.bak1 -e "s/novel_gene/PASA2_gene/g" $FINALGFF
perl -p -i.bak2 -e "s/(temp)|(novel)_model_/PASA2_model_/g" $FINALGFF
../fix_dup_ids.pl $FINALGFF
cp $FINALGFF $FINALGFF.bak3
mv $FINALGFF.updated $FINALGFF

#comprehensive transcriptome build
$PASADIR/scripts/build_comprehensive_transcriptome.dbi -c alignAssembly.config -t Trinity.$SPEC.fa.clean --min_per_ID 95 --min_per_aligned 30
TIME=$(date)
echo "[$TIME] finished transcriptome build" >> $LAUNCHDIR/$SPEC.timing

#reset dir
cd $LAUNCHDIR

