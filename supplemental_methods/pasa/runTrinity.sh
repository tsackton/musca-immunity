#!/bin/sh
#SBATCH -p general
#SBATCH --mem=200000
#SBATCH -e mdom_trin.err
#SBATCH -o mdom_trin.out
#SBATCH -n 6
#SBATCH -J trin_mdom

module load centos6/perl5mods
module load centos6/trinityrnaseq_r20131110
SEQLIST=$(ls fem_*.fastq.gz)

Trinity.pl --seqType fq --JM 100G --single $SEQLIST --output trinout --CPU 6 --min_contig_length 150 --min_kmer_cov 2
