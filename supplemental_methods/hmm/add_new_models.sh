#!/bin/bash

#SBATCH -p general
#SBATCH -t 72:00:00
#SBATCH -n 1
#SBATCH --mem 8000
#SBATCH -e fsa_%j.err
#SBATCH -o fsa_%j.out

i=$1
fsa --stockholm --fast --noindel2 --logtime $i.fa > $i.fast.aligned
sed '1a\
#=GF ID '"${i%%.*}"'
' $i.fast.aligned > $i.fixed
