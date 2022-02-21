#!/bin/bash
cd /home/dw73x
TODAY=$(TZ="Europe/London" date +%F)
REMOTE_PATH="/cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard/"
if ssh climb-covid19-wrightd@bham.covid19.climb.ac.uk  "[ -f ${REMOTE_PATH}${TODAY}.tar.gz ]"
then
    	echo "Copying remote file ${REMOTE_PATH}${TODAY}.tar.gz"
	scp "climb-covid19-wrightd@bham.covid19.climb.ac.uk:${REMOTE_PATH}${TODAY}.tar.gz" .
	echo "Extracting files"
	tar -xvzf "${TODAY}.tar.gz"
	echo "Removing archive"
	rm "${TODAY}.tar.gz"
	echo "Removing unused RDS files"
	cd "/home/dw73x/COG-UK/${TODAY}"
	rm mutations_uk.rds insertions_mapping.rds deletions_mapping.rds insertions.rds
else
 	echo "No remote file (or ssh failure)"
	exit 1
fi
