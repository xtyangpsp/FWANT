#!/bin/bash
# List the unfinished simulation folders of a work list (default: work.lst).
# Works for any layout made by skel2station_fxyz.sh (fz, fx+fy or all three),
# where each line of the list is a run folder such as  NET.STA/fx , NET.STA/fy
# or NET.STA/fz , and for lists made by the old skel2station.sh.  A run is finished when its
# seis3d_wave.log contains the word "finish".
#
# usage: ./search_unfinishedjobs_fxyz.sh [worklist]
#   stdout : one unfinished folder per line, then a ">>>>>>> ..." summary line
#   Entries starting with "#!" are read with the "#!" removed, as before.
#   Exit status: 0 = everything finished, 1 = something unfinished, 2 = no list.

list=${1:-work.lst}
if [ ! -f "$list" ]; then
    echo "work list $list not found" >&2
    exit 2
fi

nbad=0
declare -A nbad_comp ntot_comp
while read -r work; do
    work=`echo "$work" | sed -e 's/#!//g'`
    [ -z "$work" ] && continue
    comp=`basename "$work"`
    ntot_comp[$comp]=$(( ${ntot_comp[$comp]:-0} + 1 ))

    count=0
    if [ -f "${work}/seis3d_wave.log" ]; then
        count=`grep -c finish "${work}/seis3d_wave.log"`
    fi
    if [ "${count}" -eq 0 ]; then
        echo "${work}"
        (( nbad = nbad + 1 ))
        nbad_comp[$comp]=$(( ${nbad_comp[$comp]:-0} + 1 ))
    fi
done < "$list"

# per component summary (fx, fy, fz)
detail=""
for comp in `echo "${!ntot_comp[@]}" | tr ' ' '\n' | sort`; do
    detail="$detail $comp:${nbad_comp[$comp]:-0}/${ntot_comp[$comp]}"
done
echo ">>>>>>>" ${nbad} "unfinished runs  (unfinished/total per component:${detail} )"

[ "$nbad" -eq 0 ]
