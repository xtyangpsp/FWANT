#!/bin/csh
## submit the next job in the list
## Usage: ./subjob_next.sh n

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
  if ( `grep '^[^#]' $list | wc -l` == 0 ) then
    echo "Jobs finished"
    /bin/rm $ltmp
    exit 0
  endif
  set nrd = `grep -n '^[^#]' $list | head -1 | awk -F: '{print $1}'`
  set wdir = `grep -n '^[^#]' $list | head -1 | awk -F: '{print $2}'`
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
  #sbatch -A xtyang submit_wave_mpi.sh
  sbatch -A standby -t 0-03:00 submit_wave_mpi.sh
  cd $cdir
  set wdir=`echo $wdir | sed 's/\//\\\//g'`
  ex - $list << EOF
  1,\$s/\($wdir\)/\#\!\1
  wq
EOF
@ i = $i + 1
sleep 1
end

/bin/rm $ltmp
exit 0

