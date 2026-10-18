#! /bin/bash
############################################################
### Love-wave (TT) version of kern_measure2kernelconf.sh
### Script to generate the configuration files before the
### kernels calculation.
###
### Differences from the ZZ version:
###  * the picks are made on the ROTATED component 1 (rotated x', which is T):
###    picknum = 0 0 0 <npick> 0 0, followed by rot_mat
###  * rot_mat (3x3, listed column by column) holds the unit vectors, in the
###    code frame (x = colatitude/south, y = east, z = up), of
###        column 1 = T at the receiver = (-sin b, -cos b, 0)
###        column 2 = R at the receiver = ( cos b, -sin b, 0)
###        column 3 = Z                  = (0, 0, 1)
###    with b = azimuth receiver -> source (back azimuth)
###  * two numbers are appended: (sin a, cos a) with a = azimuth source ->
###    receiver.  They are the weights of the fx and fy EVENT runs that form
###    the source-side transverse force.  The current kernel code ignores them
###    (it stops reading after rot_mat); the patched code reads them.
###  * azimuths come from <syn.seismograms>/<statref>/<statref>.azimuth.txt,
###    written by snap2sac_sphere_fxfy.m, so geometry and synthetics agree.
###
### The synthetic files listed (BHN/BHE/BHZ slots) are the x,y,z responses
### to the source-side transverse force (see kern_synthetic4kernel_fxfy.m).
############################################################

ite=ite_0.025deg_01
wkdir=/depot/xtyang/data/projects/xtyang/RioGrande/FWANT_TT/$ite    # CHANGE for the TT project
listdir=$wkdir/measure
echo $listdir

# raw delay measurement file
listfile=$listdir/${ite}_measure_result_riogrande.dat

# azimuth files written by snap2sac_sphere_fxfy.m
synthdir=$wkdir/syn.seismograms

# directory name used for the kernel (pick) output.  The ZZ workflow uses BHZ;
# BHT avoids mixing up Love and Rayleigh kernels.  Set it to BHZ if downstream
# scripts (assembly, inversion lists such as inv_cmp_list) hard-code that name.
pickcmp=BHT

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
if [ -f "$conffile" ]
then
    continue
fi

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

		# azimuths for this pair: statref is the (virtual) source, statmes the receiver
		azfile=$synthdir/$statref/$statref.azimuth.txt
		if [ ! -f "$azfile" ]
		then
			echo "WARNING: $azfile not found, skipping pair $statref -> $statmes"
			continue
		fi
		read alpha beta <<< "$(awk -v r="$statmes" '$1==r {print $2, $3; exit}' "$azfile")"
		if [ -z "$alpha" ] || [ -z "$beta" ]
		then
			echo "WARNING: no azimuth for $statmes in $azfile, skipping pair $statref -> $statmes"
			continue
		fi
		# rot_mat (column by column: T, R, Z at the receiver) and source weights (sin a, cos a)
		rotline=`awk -v a="$alpha" -v b="$beta" 'BEGIN{pi=atan2(0,-1); d=pi/180;
			sa=sin(a*d); ca=cos(a*d); sb=sin(b*d); cb=cos(b*d);
			printf "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f", -sb,-cb,0, cb,-sb,0, 0,0,1, sa,ca}'`

		echo "# station correlation " $statref >> $conffile
		snfile='./'$statmes'/'$statref'/BHN/synthetic.vel.disp.dat'
		sefile='./'$statmes'/'$statref'/BHE/synthetic.vel.disp.dat'
		szfile='./'$statmes'/'$statref'/BHZ/synthetic.vel.disp.dat'

		awk -F/ '{if($1 == "'$statref'") print $0}' .recs > .subrecs
		nsubrec=`wc -l .subrecs | awk '{print $1}'`

		# picknum(1:6) = 0 0 0 <npick> 0 0 : picks on rotated component 1 (T), then rot_mat, then source weights
		echo $statref "0 0 0" $nsubrec "0 0" $rotline >> $conffile
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
				suffixe='./'$statmes'/'$statref'/'$pickcmp'/'$filter'/'$tt'.'$pf
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
