#!/bin/bash

for SAMPLE in effector_all effector_filt modulation_all modulation_filt none_all none_filt recognition_all recognition_filt signaling_all signaling_filt
do
	mkdir -p $SAMPLE
	cp run.cafe $SAMPLE/$SAMPLE.cafe
	cp $SAMPLE.tab $SAMPLE/
	cd $SAMPLE
	perl -p -i -e "s/SAMPLE_REPLACE/$SAMPLE/g" $SAMPLE.cafe
	chmod +x $SAMPLE.cafe
	nohup ./$SAMPLE.cafe 1> $SAMPLE.stdout 2> $SAMPLE.stderr &
	cd ..
done
