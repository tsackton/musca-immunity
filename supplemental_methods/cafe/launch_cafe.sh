#!/bin/bash

for SAMPLE in effector_all effector_filt modulation_all modulation_filt none_all none_filt recognition_all recognition_filt signaling_all signaling_filt
do
	cp run.cafe $SAMPLE.cafe
	perl -p -i -e "s/SAMPLE_REPLACE/$SAMPLE/g" $SAMPLE.cafe
	mkdir ${SAMPLE}_sim
	chmod +x $SAMPLE.cafe
	./$SAMPLE.cafe &
done