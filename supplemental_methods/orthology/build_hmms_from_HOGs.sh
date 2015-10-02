#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 8000
#SBATCH --time 4:00:00
#SBATCH -e aln_hmm_%j.err
#SBATCH -o aln_hmm_%j.out
#SBATCH -J aln_hmm
#SBATCH --array=100-300

INDEX=$SLURM_ARRAY_TASK_ID
export LD_PRELOAD=

echo -e "Aligning with mafft and making HMMs.\n"
echo -e "This in part $INDEX and will process HOGs/$INDEX/*.fa\n"
echo -e "This script runs the following commands:\n"
echo -e "mafft --globalpair --maxiterate 1000 --thread 8 INPUT > OUTPUT"
echo -e "hmmbuild --amino --cpu 7 -n HOG -o LOG OUTPUT INPUT"

mkdir -p global_align/${INDEX}
mkdir -p hmms/${INDEX}
mkdir -p logs/${INDEX}

for INPUT in $(ls HOGs/${INDEX}/*.fa)
do
	FILE=${INPUT##*/}
	ALN=${FILE%%.*}
	SEQS=$(grep -c ">" $INPUT)
	#one sequence: don't align, just copy file
	if [ $SEQS -eq 1 ]
	then
		if [ ! -s "global_align/$INDEX/$ALN.aligned"]
		then
			cp $INPUT global_align/$INDEX/$ALN.aligned
		fi
	fi
	#two seqs: align and make hmm
	if [ $SEQS -gt 1 ]
	then
		if [ ! -s "global_align/$INDEX/$ALN.aligned" ]
		then
			mafft --globalpair --maxiterate 1000 --thread 8 $INPUT > global_align/$INDEX/$ALN.aligned	
		fi
		if [ ! -s "hmms/$INDEX/$ALN.*" ]
		then
			hmmbuild --amino --cpu 7 -n $ALN -o logs/${INDEX}/$ALN.summary hmms/$INDEX/$ALN.hmm global_align/$INDEX/$ALN.aligned
		fi	
	fi
done
