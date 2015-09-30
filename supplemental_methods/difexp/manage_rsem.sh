#!/usr/bin/env bash

#
# The sbatch option `--dependency' can be used to coordinate when jobs run.  
# This script shows a simple example of requiring two jobs to run in 
# non-overlapping, sequential order.  See the sbatch(1) man page for more 
# details.
#
# The way this is coded, stderr still goes to the screen or whatever is calling 
# this, and the script quits if any commands fail.  The command
#
#    awk '{print $NF}'
#
# prints the last field of a string, e.g. "JOBID" from 
# "Submitted batch job JOBID".
#

# stop if any commands fail
set -e 

#define species to use
SPE=$1

jobid1=$(sbatch prep_rsem.sh $SPE | awk '{print $NF}')
#echo "job1 is: $jobid1"

for F in $(ls ../trimseq/*.fastq.gz)
do
	SEQ=${F##*/}
	sbatch --partition=general --nodes=1 --dependency=afterok:$jobid1 --ntasks-per-node=8 --mem=16000 --time=12:00:00 --output=rsemexp.%j.out --error=rsemexp.%j.err -J "$SPE.rsem" run_rsemSE.sh $SEQ $SPE
done

