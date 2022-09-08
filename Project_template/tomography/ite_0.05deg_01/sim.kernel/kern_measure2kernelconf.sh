#! /bin/bash
############################################################
### Script to generate the configuration files before the
### kernels calculation
############################################################

ite=ite_0.05deg_01
wkdir=/depot/xtyang/data/projects/xtyang/craton/$ite
#wkdir=/scratch/bell/xtyang/FWANT_craton/${ite}
listdir=$wkdir/measure
echo $listdir

# raw delay measurement file
listfile=$listdir/ite_0.05deg_01_measure_result_craton_dt7cc0.75snr7.dat

###################################################################
################# End of global parameters.
################# Normally, no need to change the parameters below.
##################################################################
# adjust for time shift in the source time function (see SeisSource.conf in sim.station)
sourceconf=$wkdir/sim.station/skel/fx/SeisSource.conf

tshift=`cat $sourceconf | grep force_stf_timefactor | awk '{print $3}'`
echo "tshift = " $tshift

statmess=`awk -F_ '{print $3}' $listfile | sort | uniq`

for statmes in $statmess
do
conffile=$statmes'_conf'
#if [ -f "$conffile" ]
#then
#    continue
#fi

echo "# pick for "$statmes > $conffile

awk -F_ '{if($3 == "'$statmes'") print $0}' $listfile | awk '{print $1":"$2":"$3":"$4":"$5":"$6":"$7}' > .recs
nrec=`wc -l .recs | awk '{print $1}'`

echo "Number of records for station" $statmes ":" $nrec

if [ "$nrec" -ge 1 ] 
	then
	statrefs=`awk -F/ '{print $1}' .recs | sort | uniq`

	for statref in $statrefs
	do
		#echo $statmes $statref
		echo "# station correlation " $statref >> $conffile
		snfile='./'$statmes'/'$statref'/BHN/synthetic.vel.disp.dat'
		sefile='./'$statmes'/'$statref'/BHE/synthetic.vel.disp.dat'
		szfile='./'$statmes'/'$statref'/BHZ/synthetic.vel.disp.dat'

		awk -F/ '{if($1 == "'$statref'") print $0}' .recs > .subrecs
		nsubrec=`wc -l .subrecs | awk '{print $1}'`

		echo $statref "0 0 "$nsubrec "0 0 0 0 0 0 0 0 0 0 0 0" >> $conffile
		echo $snfile >> $conffile
		echo $sefile >> $conffile
		echo $szfile >> $conffile

		if [ "$nsubrec" -ge 1 ]
			then
			for subrec in `awk '{print $0}' .subrecs`
			do
				filter=`echo $subrec | awk -F: '{print $7}' | sed 's/f//g'`
				t1=`echo $subrec | awk -F: '{print $5 + "'$tshift'"}'`
				t2=`echo $subrec | awk -F: '{print $6 + "'$tshift'"}'`
				pf=`echo $subrec | awk -F/ '{print $3}' | awk -F_ '{print $3}' | awk -F. '{print $2}'`
				tt=`echo $subrec | awk -F/ '{print $3}' | awk -F_ '{print $3}' | awk -F. '{print $4}'`
				suffixe='./'$statmes'/'$statref'/BHZ/'$filter'/'$tt'.'$pf
				echo $filter $pf $t1 $t2 $suffixe | awk '{printf "%s %s %7.3f %7.3f %s\n", $1, $2, $3, $4, $5}' >> $conffile
			done
		fi

	done
	# for statref

fi
# nrec >=1

done
# for statmes

rm -f .recs
rm -f .subrecs
