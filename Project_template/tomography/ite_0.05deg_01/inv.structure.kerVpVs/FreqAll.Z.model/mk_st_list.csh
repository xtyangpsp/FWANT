#!/bin/csh
# generate the inv_st_list and inv_et_list used in ambient noise tomography
#

set ite = ite_0.05deg_04
set wkdir = /home/xtyang/depot-projects/xtyang/WNACraton/$ite
set kendir = $wkdir/sim.kernel
set targetdir = $wkdir/inv.structure.kerVpVs/FreqAll.Z.model
set stnlst = $targetdir/inv_st_list
set evtlst = $targetdir/inv_ev_list
set tmplst = $targetdir/tmplst
set uniqsrclst = $targetdir/tmplst2

cat /dev/null > $stnlst
cat /dev/null > $tmplst

cd $kendir
@ ii=0
foreach conf (`/bin/ls *_conf`)
@ ii = $ii + 1
set stn = `echo $conf | awk -F_ '{print $1}'`
echo $stn $ii
echo $stn $ii >> $stnlst

#get event list. will sort to unique later
grep "synthetic" $conf | grep "/BHZ/" | sed -e 's/\// /g'|awk '{print $3}' >> $tmplst
end

cat $tmplst | sort | uniq > $uniqsrclst

cat /dev/null > $evtlst
@ ii = 0
foreach srcstn (`cat $uniqsrclst | awk '{print $1}'`)
@ ii = $ii + 1
echo $srcstn $ii >> $evtlst
end

/bin/rm $tmplst $uniqsrclst


