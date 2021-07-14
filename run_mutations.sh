#!/bin/bash
TZ="Europe/London" date
cd /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard
source /cephfs/covid/bham/climb-covid19-wrightd/miniconda3/etc/profile.d/conda.sh
conda activate r_env
Rscript mutations.R
if [ $? -eq 0 ]
	then
		TODAY=$(TZ="Europe/London" date +%F)
		echo "Creating archive ${TODAY}.tar.gz"
		tar -czvf "${TODAY}.tar.gz" "COG-UK/${TODAY}" COG-UK/vui_voc.rds
		echo "Done"
fi
