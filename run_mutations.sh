#!/bin/bash
TZ="Europe/London" date
cd /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard
source /cephfs/covid/bham/climb-covid19-wrightd/miniconda3/etc/profile.d/conda.sh
conda activate r_env
Rscript mutations.R
if [ $? -eq 0 ]
	then
   		
		TODAY=$(TZ="Europe/London" date +%F)
		
		echo "Sorting and indexing GFF3 file"
		conda activate gt_env
		cd "COG-UK/${TODAY}"
		gt gff3 -sortlines -tidy spike_epitopes.gff3 >  spike_epitopes_sorted.gff3
		bgzip spike_epitopes_sorted.gff3
		tabix -p gff spike_epitopes_sorted.gff3.gz
		rm spike_epitopes.gff3

		echo "Creating archive ${TODAY}.tar.gz"
		cd /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard
		tar -czvf "${TODAY}.tar.gz" "COG-UK/${TODAY}" COG-UK/vui_voc.rds
		TZ="Europe/London" date	
		echo "Done"
fi
