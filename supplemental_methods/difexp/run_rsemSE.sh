#!/bin/bash

module load rsem
module load bowtie2

SEQ=$1
SPE=$2
SAMP=${SEQ%%.*}

rsem-calculate-expression -p 8 --calc-ci --bowtie2 --estimate-rspd --ci-memory 6000 --time <(zcat ../trimseq/$SEQ) $SPE $SAMP
