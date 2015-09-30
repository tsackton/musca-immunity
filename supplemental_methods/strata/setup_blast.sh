#!/bin/bash

SPEC=$1

mkdir -p $SPEC
SPECFILE=$1

for DBFILE in $(ls *.fa.gz)
do
	DB=${DBFILE%%.*}
	if [ ! -f "$SPEC/${SPEC}_results_${DB}.out" ]
	then
		sbatch run_blast.sh $SPECFILE $DB
	fi
done
