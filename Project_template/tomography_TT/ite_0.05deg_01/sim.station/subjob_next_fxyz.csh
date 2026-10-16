#!/bin/csh
## submit the next job(s) in the list work.lst
## Works for any mix of run folders  NET.STA/fx  NET.STA/fy  NET.STA/fz
## (the list written by skel2station_fxyz.sh).
##
## Usage: ./subjob_next_fxyz.csh n [queue]
##   n      number of jobs to submit
##   queue  optional, appended to the "#SBATCH -A" line of each submit script
##
## Choose which force components are submitted with $COMPS (default: all three):
##   ./subjob_next_fxyz.csh 10               next 10 jobs of any component
##   setenv COMPS "fz";    ./subjob_next_fxyz.csh 10     Rayleigh (ZZ) runs only
##   setenv COMPS "fx fy"; ./subjob_next_fxyz.csh 10     Love (TT) runs only
## Jobs of other components stay in the list untouched.
## A submitted job is marked in the list by a leading "#!".

unalias cd
set cdir = `pwd`
set list = 'work.lst'
set ltmp = ".work.lst.tmp"

if ( $#argv <= "0" ) then
  echo "Too few argument"
  exit 1
else if ( $#argv > "2" ) then
  echo "Too many arguments"
  exit 1
endif

# components to submit
if ( $?COMPS ) then
  set comps = ( $COMPS )
else
  set comps = ( fx fy fz )
endif
foreach c ( $comps )
  if ( "$c" != "fx" && "$c" != "fy" && "$c" != "fz" ) then
    echo "ERROR: unknown component '$c' (use fx, fy and/or fz)"
    exit 1
  endif
end
# grep pattern for the lines  <path>/fx|fy|fz  at the end of a "number:path" line
set cre = `echo $comps | sed 's/ /|/g'`
set pat = ":.*/($cre)"'$'

@ n = $argv[1]
@ i = 1
@ cter = 0
while ( $cter < "5" )
  if ( -e $ltmp ) then
    sleep 1
    @ cter = $cter + 1
  else
    cat /dev/null > $ltmp
    break
  endif
end
if ( $cter == "5" ) then
  echo "Error: Failed to visit file $list"
  exit 1
endif
while ( $i <= $n )
  # not-yet-submitted jobs of the wanted components
  if ( `grep -n '^[^#]' $list | grep -E "$pat" | wc -l` == 0 ) then
    echo "Jobs finished (components: $comps)"
    /bin/rm $ltmp
    exit 0
  endif
  set wdir = `grep -n '^[^#]' $list | grep -E "$pat" | head -1 | awk -F: '{print $2}'`
  cd $wdir
  if ( $#argv == "2" ) then
    set queue = "$argv[2]"
    set pbs_old = `grep "^#SBATCH -A" submit_wave_mpi.sh`
    set pbs_new = "${pbs_old}:$queue"
    ex - submit_wave_mpi.sh << EOF
    1,\$s/${pbs_old}/${pbs_new}/g
    wq
EOF
  endif
  sbatch -q standby -t 0-04:00 submit_wave_mpi.sh
  cd $cdir
  # mark this job as submitted: match the whole line only (so NET.STA/fx never
  # matches inside another station name or another component)
  set wdir = `echo $wdir | sed -e 's/[.[\*^$]/\\&/g' -e 's/\//\\\//g'`
  ex - $list << EOF
  1,\$s/^\($wdir\)\$/\#\!\1/
  wq
EOF
@ i = $i + 1
sleep 1
end

/bin/rm $ltmp
exit 0
