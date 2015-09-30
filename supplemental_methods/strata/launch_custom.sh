#!/bin/bash

for F in $(ls *.fa.gz)
do
	SPE=${F%%.*}
	sbatch ./makeblastdb.sh $SPE
done
