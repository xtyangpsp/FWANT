#!/bin/bash

list_sta="../../STinfo/craton_station_withdata.txt"
template="skel/fx"
scriptfile="submit_wave_mpi.sh"
echo $list_sta
echo $template

cat /dev/null > rest.lst0
for KSTNM in `awk '{print $1 "." $2}' ${list_sta}`; do
    KNWNM=`echo $KSTNM | awk -F. '{print $1}'`
    KEPNM=`echo $KSTNM | awk -F. '{print $2}'`
    echo $KSTNM $KNWNM $KEPNM

    if [ -e $KSTNM ]; then
      echo '$KSTNM exists'
      continue
    fi

    mkdir -p $KSTNM

    cp -R ${template} ${KSTNM}/fz

    # SeisSource.conf
    # for spherical coord.
    lon=`grep $KNWNM $list_sta | grep $KEPNM | awk '{print $3}'`
    lat=`grep $KNWNM $list_sta | grep $KEPNM | awk '{print $4}'`
	 
    colat=`echo "90 $lat" | awk '{print $1-$2}'`

    echo $colat $lon

    echo "$colat $lon 2e10   0.0  1.0e+16   0.0 0.0 1.0" >> ${KSTNM}/fz/SeisSource.conf

    cat /dev/null > ${KSTNM}/fz/${scriptfile}
    cat ${template}/${scriptfile} | while read -r strline; do
        if [ `echo $strline | grep ^"#SBATCH -J" | wc -l` -eq 1 ]; then
           strline="#SBATCH -J $KSTNM"
           echo "${strline}.fz" >> ${KSTNM}/fz/${scriptfile}
        else
           echo "$strline" >> ${KSTNM}/fz/${scriptfile}
        fi
    done
    echo "${KSTNM}/fz" >> rest.lst0
    cd ${KSTNM}/fz
    ./bin/seis3d_source
    cd ../..
done

/bin/cp rest.lst0 work.lst

