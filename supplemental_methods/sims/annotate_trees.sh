#!/bin/bash

for SAMP in 0.00057 0.00114
do
	echo "Working on $SAMP";
	for PERM in $(seq 0 999)
	do
		if [ -s run2/$SAMP/sim.$PERM.pruned.leafmap ]
		then
			tree-annotate -s ultra.final.nwk -S run2/$SAMP/sim.$PERM.pruned.leafmap -T .tree run2/$SAMP/sim.$PERM.pruned.tree
		fi
	done
done

