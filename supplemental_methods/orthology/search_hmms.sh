#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 32000
#SBATCH --time 12:00:00
#SBATCH -e search_hmm_%j.err
#SBATCH -o search_hmm_%j.out
#SBATCH -J search_hmm
#SBATCH --array=1-2

INDEX=$SLURM_ARRAY_TASK_ID

PROFILE_NUM=$(ls hmms/*/*.hmm | wc -l)
SEQ_NUM=$(grep -c ">" allprot.fa)
DB_SIZE=$(($PROFILE_NUM * $SEQ_NUM))

mkdir -p hits/${INDEX}r

for INPUT in $(ls hmms/${INDEX}*/*.hmm)
do
	FILE=${INPUT##*/}
	ALN=${FILE%%.*}
	DONE=$(grep -c "$ALN" finished.hits.2) 
	if [ $DONE -eq 0 ]
	then
		hmmsearch -o /dev/null --tblout hits/$INDEX/$ALN.out -E 0.01 -Z $DB_SIZE --cpu 7 $INPUT allprot.fa
	fi
done
