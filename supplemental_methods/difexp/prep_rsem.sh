#!/bin/bash

#SBATCH -p general
#SBATCH --mem 1000
#SBATCH -n 1
#SBATCH -J rsem
#SBATCH -o prep_rsem_%j.out
#SBATCH -e prep_rsem_%j.err
#SBATCH --time=1:00:00

module load bowtie2
module load rsem

SPE=$1

rsem-prepare-reference --transcript-to-gene-map ../annot/$SPE.final.rnakey --bowtie2 ../annot/$SPE.final.transcripts.fa $SPE

