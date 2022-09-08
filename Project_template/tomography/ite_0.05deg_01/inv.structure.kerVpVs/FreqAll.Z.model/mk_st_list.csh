#!/bin/csh
# generate the inv_st_list and inv_et_list used in ambient noise tomography
#

set ite = ite_0.05deg_01
set wkdir = /home/xtyang/depot-projects/xtyang/craton/$ite
set kendir = $wkdir/sim.kernel
set targetdir = $wkdir/inv.structure.kerVpVs/FreqAll.Z.model
set stnlst = $targetdir/inv_st_list
set evtlst = $targetdir/inv_ev_list
set tmplst = $targetdir/tmplst
set uniqsrclst = $targetdir/tmplst2

cat /dev/null > $stnlst
cd $kendir
@ ii=0
foreach conf (`/bin/ls *_conf`)
@ ii = $ii + 1
set stn = `echo $conf | awk -F_ '{print $1}'`
echo $stn $ii >> $stnlst
end

cat /dev/null > $tmplst
foreach rcvstn (`cat $stnlst | awk '{print $1}'` )
cd $rcvstn
foreach srcstn (`/bin/ls -d *.*`)
echo $srcstn
echo $srcstn >> $tmplst
end
cd ..
end

cat $tmplst | sort | uniq > $uniqsrclst

cat /dev/null > $evtlst
@ ii = 0
foreach srcstn (`cat $uniqsrclst | awk '{print $1}'`)
@ ii = $ii + 1
echo $srcstn $ii >> $evtlst
end

/bin/rm $tmplst $uniqsrclst


