#!/bin/bash

for RATE in 00057 00114 002 005 007 009 014 023 034 057 114 341
do
	head -1 sim_ogs_ct.txt > sim
	egrep "^sim\.\d+\.$RATE" sim_ogs_ct.txt >> sim
	perl -p -i -e 's/OGS/OGS\tdesc/' sim
	perl -p -i -e "s/$RATE/$RATE\tnone/" sim
	cp sim sim$RATE.txt
	./run_sim.cafe
	mv sim.log sim$RATE.log
	rm sim
done