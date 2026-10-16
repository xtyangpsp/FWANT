#!/bin/bash
# Print the job names of the unfinished runs, one per line, in the form used
# by the submit scripts of skel2station_fxyz.sh:  NET.STA.fx  /  NET.STA.fy  /  NET.STA.fz
# (the folder NET.STA/fx becomes NET.STA.fx).
#
# usage: ./get_unfinished_joblist_fxyz.sh [worklist]     (default work.lst)
dir=`dirname "$0"`
bash "$dir"/search_unfinishedjobs_fxyz.sh "$@" | grep -v '^>>>>>>>' | sed -e 's/\//./g' | sort
exit 0
