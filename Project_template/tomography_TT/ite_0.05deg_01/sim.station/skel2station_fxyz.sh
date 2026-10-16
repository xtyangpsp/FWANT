#!/bin/bash
# Prepare the simulations of every virtual source from one skeleton.
# One run folder per force component is created for each station:
#     <NET.STA>/fx   force direction (1.0 0.0 0.0)   x = colatitude (south)
#     <NET.STA>/fy   force direction (0.0 1.0 0.0)   y = east
#     <NET.STA>/fz   force direction (0.0 0.0 1.0)   z = up
# Which of them are made is chosen with the first argument (or $COMPS):
#     ./skel2station_fxyz.sh "fz"          Rayleigh (ZZ) only, same as the old skel2station.sh
#     ./skel2station_fxyz.sh "fx fy"       Love (TT): two horizontal forces   (default)
#     ./skel2station_fxyz.sh "fx fy fz"    all three
# All other source parameters are unchanged (depth 2e10, start 0.0, f0 1.0e+16).
#
# The skeleton (default skel/fx, override with $TEMPLATE) must hold the files of
# a run (SeisFD3D.conf, bin/, input/, the submit script, ...) and a
# SeisSource.conf WITHOUT any force line after <anchor_force>; the force line
# is appended here.

list_sta="../../STinfo/craton_station_withdata.txt"
template=${TEMPLATE:-skel/fx}
scriptfile="submit_wave_mpi.sh"
comps=${1:-${COMPS:-"fx fy"}}
echo $list_sta
echo $template
echo "components: $comps"

# only fx, fy and fz are known
for comp in $comps; do
    case $comp in
        fx|fy|fz) ;;
        *) echo "ERROR: unknown component '$comp' (use fx, fy and/or fz)"; exit 1 ;;
    esac
done

# check that the skeleton SeisSource.conf has no force line yet
nforce=`awk '/^<anchor_force>/{f=1;next} f && NF>0 && $1 !~ /^#/ {n++} END{print n+0}' ${template}/SeisSource.conf`
if [ "$nforce" -ne 0 ]; then
    echo "ERROR: ${template}/SeisSource.conf already has $nforce force line(s) after <anchor_force>;"
    echo "       remove them, this script appends the force line."
    exit 1
fi

cat /dev/null > rest.lst0
for KSTNM in `awk '{print $1 "." $2}' ${list_sta}`; do
    KNWNM=`echo $KSTNM | awk -F. '{print $1}'`
    KEPNM=`echo $KSTNM | awk -F. '{print $2}'`
    echo $KSTNM $KNWNM $KEPNM

    # station location (exact match on network and station name)
    lon=`awk -v n=$KNWNM -v s=$KEPNM '$1==n && $2==s {print $3; exit}' $list_sta`
    lat=`awk -v n=$KNWNM -v s=$KEPNM '$1==n && $2==s {print $4; exit}' $list_sta`
    # for spherical coord.
    colat=`echo "90 $lat" | awk '{print $1-$2}'`
    echo $colat $lon

    for comp in $comps; do
        case $comp in
            fx) fdir="1.0 0.0 0.0" ;;
            fy) fdir="0.0 1.0 0.0" ;;
            fz) fdir="0.0 0.0 1.0" ;;
        esac

        if [ -e ${KSTNM}/${comp} ]; then
            echo "${KSTNM}/${comp} exists"
            continue
        fi

        mkdir -p ${KSTNM}
        cp -R ${template} ${KSTNM}/${comp}

        # SeisSource.conf: x(colat) y(lon) z  start  f0  force direction
        echo "$colat $lon 2e10   0.0  1.0e+16   $fdir" >> ${KSTNM}/${comp}/SeisSource.conf

        # submit script with the job name <NET.STA>.<comp>
        cat /dev/null > ${KSTNM}/${comp}/${scriptfile}
        cat ${template}/${scriptfile} | while IFS= read -r strline; do
            if [ `echo "$strline" | grep ^"#SBATCH -J" | wc -l` -eq 1 ]; then
               echo "#SBATCH -J ${KSTNM}.${comp}" >> ${KSTNM}/${comp}/${scriptfile}
            else
               echo "$strline" >> ${KSTNM}/${comp}/${scriptfile}
            fi
        done
        echo "${KSTNM}/${comp}" >> rest.lst0
        cd ${KSTNM}/${comp}
        ./bin/seis3d_source
        cd ../..
    done
done

/bin/cp rest.lst0 work.lst
