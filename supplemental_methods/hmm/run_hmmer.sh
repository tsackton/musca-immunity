#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 4000
#SBATCH -p general
#SBATCH -o hmmr_imm.out
#SBATCH -e hmmr_imm.err
#SBATCH -J hmmrimm
#SBATCH -t 24:00:00

module load centos6/hmmer-3.1b1

hmmbuild immune.hmm immune_final.sto
hmmsearch -o imm_full.out --domtblout dom_out.tab --tblout target_out.tab --noali --notextw -E 1 --domE 1 --incE 0.001 --cpu 16 immune.hmm ../compara2/all_prot.fa

