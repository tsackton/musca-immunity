#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 2000
#SBATCH -p general
#SBATCH -o makedb_%j.out
#SBATCH -e makedb_%j.err
#SBATCH -J makedb
#SBATCH -t 3:00:00

#blast input
SEQDB=$1

#make blast db
zcat $SEQDB.fa.gz | segmasker -in - -parse_seqids -infmt fasta -outfmt maskinfo_asn1_bin -out $SEQDB.mask.asnb
zcat $SEQDB.fa.gz | makeblastdb -in - -parse_seqids -input_type fasta -dbtype prot -mask_data $SEQDB.mask.asnb -title $SEQDB -out $SEQDB
blastdbcmd -db $SEQDB -info
