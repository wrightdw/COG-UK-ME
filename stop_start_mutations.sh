#!/bin/bash
cd /opt/shiny-server/my-apps/COG-UK/
TZ="Europe/London" date
TODAY=$(TZ="Europe/London" date +%F)
if [ -d ./${TODAY} ];  then
    echo "$PWD/${TODAY}" exists
else
	if sudo -u dw73x /home/dw73x/update_mutations.sh;  then
		echo "Moving GFF3 and index to /srv/shiny-server/jbrowser/"
		mv -f /home/dw73x/COG-UK/${TODAY}/spike_epitopes_sorted.gff3.gz* /srv/shiny-server/jbrowser/
		
		echo "Moving /home/dw73x/COG-UK/${TODAY} to ${PWD}"
	    	mv /home/dw73x/COG-UK/${TODAY} .

		echo "Moving /home/dw73x/COG-UK/vui_voc.rds to ${PWD}"
	    	mv /home/dw73x/COG-UK/vui_voc.rds .

		TZ="Europe/London" date
		echo "Done"
	else
		exit 1
	fi
fi
