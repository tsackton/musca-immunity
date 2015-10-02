#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 8000
#SBATCH --time 12:00:00
#SBATCH -e aln_ogs_%j.err
#SBATCH -o aln_ogs_%j.out
#SBATCH -J aln_ogs
#SBATCH --array=100-300

INDEX=$SLURM_ARRAY_TASK_ID
export LD_PRELOAD=

echo -e "Aligning with mafft.\n"
echo -e "This in part $INDEX and will process new_ogs/$INDEX/*.fa\n"
echo -e "This script runs the following commands:\n"
echo -e "mafft --genafpair --maxiterate 1000 --thread 8 INPUT > OUTPUT\n"
echo -e "mafft --maxiterate 1000 --thread 8 INPUT > OUTPUT for >200 sequences"

mkdir -p aligned/${INDEX}
mkdir -p logs/${INDEX}

for INPUT in $(ls newogs_seq/${INDEX}/*.fa)
do
	FILE=${INPUT##*/}
	ALN=${FILE%%.*}
	SEQS=$(grep -c ">" $INPUT)
	#one sequence: don't align, just copy file
	if [ ! -s "aligned/$INDEX/$ALN.aligned" ]
	then
		if [ $SEQS -eq 1 ]
		#one seq: copy
		then
			cp $INPUT aligned/$INDEX/$ALN.aligned
		#two seqs: align
		else 
			if [ $SEQS -lt 200 ]
			then
				mafft --genafpair --maxiterate 1000 --thread 8 $INPUT > aligned/$INDEX/$ALN.aligned	
			else
				mafft --maxiterate 1000 --thread 8 $INPUT > aligned/$INDEX/$ALN.aligned
			fi
		fi
	fi
done
