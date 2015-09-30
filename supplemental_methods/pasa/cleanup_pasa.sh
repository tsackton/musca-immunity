#!/bin/bash

SPE=$1
LAUNCHDIR=$(pwd)

#first move to right dir
cd $SPE
mkdir -p pasa_out

#original unmodified data
OLDGFF=$SPE.coding.gff
gffread $OLDGFF -g $SPE.fa -w pasa_out/${SPE}.orig.transcripts.fa -y pasa_out/${SPE}.orig.proteins.fa
../make_isoform_file_from_GFF.pl $OLDGFF > pasa_out/${SPE}.orig.isokey.txt

#new data - only problem is that I don't check for duplicates
NEWGFF=$(ls -t *.gff3 | head -n1)
gffread $NEWGFF -g $SPE.fa -w pasa_out/${SPE}.new.transcripts.fa -y pasa_out/${SPE}.new.proteins.fa
../make_isoform_file_from_GFF.pl $NEWGFF > pasa_out/${SPE}.new.isokey.txt

#unmapped data
../parse_unmapped.pl $SPE

#copy to dropbox
cp ${SPE}_new_db.gene*.*.gff3 pasa_out
cp ${SPE}_new_db.gene*.*.bak? pasa_out
cp $OLDGFF pasa_out
tar -cf $SPE.pasa.tar pasa_out
gzip $SPE.pasa.tar
mv $SPE.pasa.tar.gz ~/Dropbox/Work

#reset
cd $LAUNCHDIR
