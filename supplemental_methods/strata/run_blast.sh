#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 1000
#SBATCH -p general
#SBATCH -o blast_%j.out
#SBATCH -e blast_%j.err
#SBATCH -J blast
#SBATCH --time-min=24:00:00

SP=$1
DB=$2
cat ../annot/$SP.final.proteins.fa | blastp -query - -task blastp -db $DB -out ./$SP/${SP}_results_$DB.out -evalue 0.001 -outfmt 6 -num_threads 16 -max_target_seqs 1 -db_soft_mask 21 -seg "yes" -soft_masking T


