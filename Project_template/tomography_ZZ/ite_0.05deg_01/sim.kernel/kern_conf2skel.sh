#!/bin/bash
###################################################################
### Script to create directories where to put the kernel files
### to create the list of stations which we would like to calculate
### the kernels
###################################################################

stlist=station_conf_list

conffiles=`ls -f *_conf`
if [ "${#conffiles}" -eq 0 ]
then
echo "*_conf files do not exist"
exit
else
ls -f *_conf | awk -F/ '{print $NF}' | awk -F_ '{print $1, $1"_"$2}' > $stlist
fi

for conffile in $conffiles
do
echo $conffile

echo "make synthetic files direcotry"
for pnm in `grep / $conffile | awk '{if(NF==1) print $1}' | sed 's/synthetic.vel.disp.dat//g' | sed 's/\.\///g' | sort | uniq`
do
if [ ! -d $pnm ] 
then
mkdir -p $pnm
echo $pnm
fi
done

echo "make kernel files direcotry"
# find in the conf file a line with five "words" to identify the time series, like the following
#  4 P2  53.000 153.000 ./TA.E07A/TA.B04A/BHZ/4/T1T2.P2
for pnm in `grep / $conffile  | awk '{if(NF==5) print $5}' | sed 's/\.\///g' | sort | uniq`
do
if [ ! -d $pnm ]
then
mkdir -p $pnm
echo $pnm
fi
done

done
