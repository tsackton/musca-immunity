#!/bin/bash

#SBATCH -p general
#SBATCH --mem 16000
#SBATCH --time 24:00:00
#SBATCH -J maketrees
#SBATCH -o raxml_%j.out
#SBATCH -e raxml_%j.err
#SBATCH --array=196
#SBATCH -N 1
#SBATCH -n 12

INDEX=$SLURM_ARRAY_TASK_ID
mkdir -p phylip/$INDEX
mkdir -p raxml/$INDEX
for INPUT in $(ls aligned/$INDEX/*.aligned)
do
	FILE=${INPUT##*/}
	ALN=${FILE%%.*}
	SEQNUM=$(grep -c ">" $INPUT)
	#first convert to phylip for RAxML
	if [ ! -s "phylip/$INDEX/$ALN.phy" ]
	then
		trimal -in aligned/$INDEX/$ALN.aligned -out phylip/$INDEX/$ALN.phy -phylip
	fi
	
	#now run RAxML if at least three species
	
	if [ -f raxml/$INDEX/RAXML_bestTree.$ALN ]
	then
		mv raxml/$INDEX/RAXML_bestTree.$ALN raxml/$INDEX/RAxML_bestTree.$ALN
	fi

	if [ ! -s raxml/$INDEX/RAxML_bestTree.$ALN ]
	then
		rm -f raxml/$INDEX/RAxML_info.$ALN
		rm -f raxml/$INDEX/RAxML_*.$ALN.RUN.*
		if [ $SEQNUM -gt 349 ]
		then
			#big tree, use -D option to speed up run
			raxmlHPC-PTHREADS-AVX -n $ALN -s phylip/$INDEX/$ALN.phy -m PROTGAMMAAUTO -p 12345 -N 10 -T 12 -D -w $HOME/regal/musca/orthopipe2/raxml/$INDEX/
		elif [ $SEQNUM -gt 3 ]
		then
			#normal tree
			raxmlHPC-PTHREADS-AVX -n $ALN -s phylip/$INDEX/$ALN.phy -m PROTGAMMAAUTO -p 12345 -N 10 -T 12 -w $HOME/regal/musca/orthopipe2/raxml/$INDEX/
		else
			#make unrooted tree from all species names
			grep ">" $INPUT | perl -pe 's/>//g' | perl -pe 's/\n/,/g' | perl -pe 's/^(.*),$/\($1\)/' > raxml/$INDEX/RAxML_bestTree.$ALN
		fi
	fi
done
