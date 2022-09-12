#!/bin/bash
#if [ $# -lt 1 ]
#	then
echo 'USAGE: make_inv_Gd_list.sh [ftag] [snrmin]'
#	exit 1
#fi
snrmin=0
if [ $# -eq 2 ]
	then
	snrmin=$2
fi
hline="-----------------------------------------------"
fnm_meas_master='../../measure/ite_0.001deg_01_measure_result_sitkin.dat'
# threshold of travel time delays (see delay_vs_dist.png)
thres_dat=15

if [ $# -eq 1 ]
        then
        ftag=$1
	fnm_meas=".measure_temp_"${ftag}
	awk '{ if ( $7 == "'${ftag}'" && $10 >= '$snrmin' ) print }' ${fnm_meas_master} > ${fnm_meas}
	list_Gd="inv_Gd_list_"${ftag}
else
	fnm_meas=".measure_temp"
	awk '{ if ( $10 >= '$snrmin' ) print }' ${fnm_meas_master} > ${fnm_meas}
	list_Gd="inv_Gd_list"
fi

############################################################
# black list: Do not use stations in this list
list_black='inv_black_list'
# list of used station ($rcvstn)
list_st='inv_st_list'
# list of source (stations $srcstn)
list_ev='inv_ev_list'
nev=`wc -l $list_ev | awk '{print $1}'`
# components (e.g., BHZ BHR BHT in separate lines)
list_cmp='inv_cmp_list'
# frequency bands
list_freq='inv_freq_list'
# time window (e.g., T1T2, T3T4, T5T6 in separate lines)
list_twin='inv_twin_list'
# arrival phases (e.g., PD #direct P; SD #direct S; RL #surface wave)
list_phase='inv_phase_list'

list_G_raw='list_G_raw'

weig0=0

if [ -e $list_Gd ] ;then
   echo "mv $list_Gd ${list_Gd}.bak"
   mv -f $list_Gd ${list_Gd}.bak
fi
count=0
# loop all dat
cat ${fnm_meas} | while read rcd; do
(( count=count+1 ))
echo $count
  #echo $rcd

# Here f1, f2 stand for first and second strings separated by "/"
# for example, for the following line in measure_results.dat 
# CD.BJI/bp0.00667_0.01333/CD.BJI_CD.SSE_BHZ.P2.CORR.T1T2.SAC  5.799  1 RL   45.00 1477.74 f5  2.62  0.93  13.1
# KEVNM == CD.BJI
# KFREQ == bp0.00667_0.01333
# KSTNM == CD.SSE (second string separated by "_' in the third string separated by "/")
  KEVNM=`echo "$rcd" | cut -d/ -f1` 
  KFREQ=`echo "$rcd" | cut -d/ -f2`
  KSTNM=`echo "$rcd" | cut -d/ -f3 | cut -d\_ -f2`
  KCMPNM=`echo "$rcd" | cut -d/ -f3 | cut -d\_ -f3 | cut -d\. -f1` 
  KFILT=`echo "$rcd" | cut -d/ -f3 | cut -d\_ -f3 | cut -d\. -f2`
  KTWIN=`echo "$rcd" | cut -d/ -f3 | cut -d\_ -f3 | cut -d\. -f4`
  dval=`echo "$rcd" | awk '{print $2}'`

  tstd=`echo "$rcd" | awk '{print $8}'`
  weig=`echo $tstd | awk '{print 1/$1}'`
  KPHNM=`echo "$rcd" | awk '{print $4}'`
  T1=`echo "$rcd" | awk '{print $5}'`
  T2=`echo "$rcd" | awk '{print $6}'`
  echo $KEVNM $KFREQ $KSTNM $KCMPNM $KFILT $KTWIN $dval $weig $KPHNM $T1 $T2
# exit

  # check black list
  #if [ `grep ^"$KEVNM" $list_black|grep "$KSTNM"|grep "$KCMPNM"|grep "$KFREQ"|grep "$KTWIN"|wc -l` -eq 1 ]; then
  if [ `grep ^"$KEVNM.*$KSTNM.*$KCMPNM.*$KFREQ.*$KTWIN.*$KFILT" $list_black|wc -l` -eq 1 ]; then
     echo $hline
     echo "  omit $rcd"
     echo "       according to $list_black"
     echo $hline
     continue
  fi

  # check event list
  if [ `awk '{ if ( $1 == "'$KEVNM'" ) print}' $list_ev | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit event $KEVNM"
     echo "       according to $list_ev"
     echo $hline
     continue
  else
     evno=`awk '{ if ( $1 == "'$KEVNM'" ) print $2}' $list_ev`
     #evno=`grep ^"$KEVNM" $list_ev | awk '{print $2}'`
  fi

  # check station list
  if [ `awk '{ if ( $1 == "'$KSTNM'" ) print }'  $list_st | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit station $KSTNM"
     echo "       according to $list_st"
     echo $hline
     continue
  else
     stno=`awk '{ if ( $1 == "'$KSTNM'" ) print $2}' $list_st`
     #stno=`grep ^"$KSTNM" $list_st |awk '{print $2}'`
  fi

  # check freq list
  if [ `grep ^"$KFREQ" $list_freq | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit frequency $KFREQ"
     echo "       according to $list_freq"
     echo $hline
     continue
  else
     IFREQ=`grep ^"$KFREQ" $list_freq |awk '{print $2}'`
  fi

  # check cmp list
  if [ `grep ^"$KCMPNM" $list_cmp | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit component $KCMPNM"
     echo "       according to $list_cmp"
     echo $hline
     continue
  fi

  # check twin list
  if [ `grep ^"$KTWIN" $list_twin | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit time window $KTWIN"
     echo "       according to $list_twin"
     echo $hline
     continue
  fi

  # check phase list
  if [ `grep ^"$KPHNM" $list_phase | wc -l` -eq 0 ]; then
     echo $hline
     echo "  omit phase $KPHNM"
     echo "       according to $list_phase"
     echo $hline
     continue
  fi

#  # check weig threshold
#  if [ `echo "$weig < $weig0"|bc` -eq 1 ]; then
#     echo $hline
#     echo "  omit small $weig"
#     echo "       according to threshold $weig0"
#     echo $hline
#     continue
#  fi

  # check measurement threshold
  if [ `echo "$dval $thres_dat"|awk '{if ($1>$2 || -$1>$2) print 1; else print 0;}'` -eq 1 ]; then
     echo $hline
     echo "  omit too large measure $dval"
     echo "       according to threshold $thres_dat"
     echo $hline
     continue
  fi

  #echo ${KSTNM} ${KEVNM} ${KCMPNM} ${IFREQ} ${KTWIN}
  fnm_k=''
  for fnm in `grep "^[^#].*${KSTNM}_${KEVNM}_${KCMPNM}.*${IFREQ}.*${KTWIN}.*${KFILT}" $list_G_raw`; do
     fnm_k="${fnm_k} $fnm"
  done
  if [ -z "$fnm_k" ]; then
     echo $hline
     echo "no kernel"
     echo $rcd
     continue
  fi

# travel time sensivitity to variation in locations
#  LOCNM="../../inv.srcmech/syn_sac_vol/loc_kernels/measure_${KEVNM}.lk.dat"
#  if [ ! -e $LOCNM ]; then
#     xlocv=0
#     ylocv=0
#     zlocv=0
#  else   
#    xlocv=`grep "${KFREQ}/${KEVNM}.${KSTNM}.${KCMPNM}.${KFILT}.CORR.${KTWIN}" $LOCNM | cut -d" " -f7`
#    ylocv=`grep "${KFREQ}/${KEVNM}.${KSTNM}.${KCMPNM}.${KFILT}.CORR.${KTWIN}" $LOCNM | cut -d" " -f8`
#    zlocv=`grep "${KFREQ}/${KEVNM}.${KSTNM}.${KCMPNM}.${KFILT}.CORR.${KTWIN}" $LOCNM | cut -d" " -f9`
#  fi

# for ambient noise tomography, source is known so the traveltime sensivitities to source locations are zero
  xlocv=0
  ylocv=0
  zlocv=0
# source index and source location indexes ($nev is the number of sources)
  xlocn=$evno 
  ylocn=`expr $xlocn + $nev` 
  zlocn=`expr $ylocn + $nev`
  echo $xlocn $ylocn $zlocn
  echo "$fnm_k $dval $weig $stno $evno $xlocn $xlocv $ylocn $ylocv $zlocn $zlocv $KPHNM"
  echo "$fnm_k $dval $weig $stno $evno $xlocn $xlocv $ylocn $ylocv $zlocn $zlocv $KPHNM" >> $list_Gd
done
