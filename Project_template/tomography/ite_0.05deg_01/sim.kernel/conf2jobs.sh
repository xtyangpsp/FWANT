#!/bin/bash
conflist='station_conf_list'
template='submit_kernel_mpi_template.sh'
assem_template='submit_assem_kernel_template.sh'

for sta in `awk '{print $1}' $conflist`
do
	echo $sta
	stalist=${sta}.list
	stasubmit=${sta}.submit
	stasubmit_assem=${sta}.assem
	
	echo $sta ${sta}_conf > ${stalist}
	
	sed -e 's/KNJOBTEMPLATE/'${sta}'.kn/g' $template | sed -e 's/STATION_CONF_LIST_TEMPLATE/'$stalist'/g' > ${stasubmit}
	chmod a+x ${stasubmit}

	#assembly script
	sed -e 's/KAJOBTEMPLATE/'${sta}'.as/g' ${assem_template} | sed -e 's/STATION_CONF_LIST_TEMPLATE/'$stalist'/g' > ${stasubmit_assem}
	chmod a+x ${stasubmit_assem}
done
