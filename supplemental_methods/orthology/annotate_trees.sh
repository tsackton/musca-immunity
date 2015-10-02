#!/bin/bash

export PATH=$HOME/Applications/python/bin:$PATH
export PYTHONPATH=$HOME/Applications/python/lib/python:$PYTHONPATH

for OG in $(cat ogs_list)
do
	if [[ ! -s duplication/$OG.rel.txt ]]
	then
		cp new/$OG.out.nwk .
		tree-annotate -s dipt.tree -S dipt.smap -T .out.nwk $OG.out.nwk
		mv $OG.paralog.txt homology/para
		mv $OG.orth.txt homology/orth
		mv $OG.rel.txt duplication
		mv $OG.mpr.recon annotate
		mv $OG.nhx.tree annotate
		rm $OG.out.nwk
	fi
done