#!/bin/bash

for SPEC in mdom
do
	STTIME=$(date)
	LOGTIME=$(date +"%F")
	echo "[$STTIME] Starting $SPEC";
	#fix gff
	./remove_id_from_feat.pl $SPEC/$SPEC.coding.gff
	mv $SPEC/$SPEC.coding.gff.new $SPEC/$SPEC.coding.gff
	#run pasa commands
	./run_pasa.sh $SPEC > ${SPEC}.$LOGTIME.stdout.log 2> ${SPEC}.$LOGTIME.stderr.log
	#cleanup and copy necessary files to Dropbox
	./cleanup_pasa.sh $SPEC	
	ENDTIME=$(date)
	echo "[$ENDTIME] Finishing $SPEC";
	tar -cf ${SPEC}_log.tar ${SPEC}.$LOGTIME.*.log 
	gzip ${SPEC}_log.tar
	mv ${SPEC}_log.tar.gz ~/Dropbox/Work
	cp ${SPEC}.timing ~/Dropbox/Work
done
