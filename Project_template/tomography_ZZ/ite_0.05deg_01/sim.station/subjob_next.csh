#!/bin/csh
## Submit the jobs listed in work.lst (one run folder per line, e.g. NET.STA/fx).
## Works for any list made by skel2station_fxyz.sh (fx, fy, fz in any mix).
##
## Usage: ./subjob_next_fxyz.csh [n] [queue]
##   n      optional number of jobs to submit; default (or 0 / all) = every job in the list
##   queue  optional, appended to the "#SBATCH -A" line of each submit script
## An all-digit argument is taken as n, any other argument as the queue.
##   ./subjob_next_fxyz.csh              submit everything not yet submitted
##   ./subjob_next_fxyz.csh 10           submit the next 10
##   ./subjob_next_fxyz.csh 10 myqueue   next 10, with queue appended to the account line
##
## A job is marked in the list by a leading "#!" ONLY after sbatch accepted it
## ("Submitted batch job <id>").  If sbatch fails the entry stays in the list and is
## retried (3 tries); an entry that still fails is reported and left unmarked, so it is
## picked up by the next run.  After 3 failed entries in a row the script stops (the
## usual reason is the submit limit or a scheduler problem).
## The marking is done by line number, so station names are never matched as patterns.

unalias cd
set cdir = `pwd`
set list = 'work.lst'
set ltmp = ".work.lst.tmp"
set sbout = "$cdir/.subjob_sbatch.$$"
set sbatch_opts = "-q standby -t 0-04:00"     # options for every sbatch call
set script = "submit_wave_mpi.sh"
set ntry_max = 3                              # sbatch tries per job
set nfail_max = 3                             # stop after this many failed entries in a row

# ---------------- arguments ----------------
if ( $#argv > 2 ) then
  echo "Too many arguments (usage: $0 [n] [queue])"
  exit 1
endif
set limit = 0
set queue = ""
foreach a ( $argv )
  if ( "$a" == "all" ) then
    set limit = 0
  else
    echo "$a" | grep -q '^[0-9][0-9]*$'
    if ( $status == 0 ) then
      set limit = $a
    else
      set queue = "$a"
    endif
  endif
end

if ( ! -f $list ) then
  echo "Error: $list not found in $cdir"
  exit 1
endif

# ---------------- lock (atomic; stale locks of dead processes are removed) ----------------
set got = 0
set lpid = ""
@ cter = 0
while ( $cter < 10 )
  set noclobber
  ( echo $$ > $ltmp ) >& /dev/null
  set st = $status
  unset noclobber
  if ( $st == 0 ) then
    set got = 1
    break
  endif
  set lpid = `cat $ltmp |& head -1`
  if ( "$lpid" != "" ) then
    ps -p $lpid >& /dev/null
    if ( $status != 0 ) then
      echo "removing stale lock $ltmp (process $lpid is not running)"
      /bin/rm -f $ltmp
      continue
    endif
  endif
  sleep 1
  @ cter = $cter + 1
end
if ( $got == 0 ) then
  echo "Error: $list is in use by another run of this script (lock $ltmp, pid $lpid)."
  exit 1
endif
onintr cleanup

# ---------------- submit ----------------
# line numbers of the entries not yet submitted (blank lines and "#..." lines are ignored)
set lines = ( `grep -n '^[^#]' $list | awk -F: '{print $1}'` )
if ( $#lines == 0 ) then
  echo "Jobs finished (no unsubmitted entry in $list)"
  goto cleanup
endif

set nsub = 0
set nbad = 0
set consec = 0
set badlist = ()
foreach ln ( $lines )
  if ( $limit > 0 && $nsub >= $limit ) break

  set wdir = `sed -n "${ln}p" $list | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
  if ( "$wdir" == "" ) continue
  if ( "$wdir" =~ /* ) then
    set full = "$wdir"
  else
    set full = "$cdir/$wdir"
  endif

  if ( ! -f "$full/$script" ) then
    echo "WARNING: $full/$script not found, entry $wdir (line $ln) left in the list"
    @ nbad = $nbad + 1
    set badlist = ( $badlist $wdir )
    continue
  endif
  cd "$full"

  # optional queue: append once to the account line
  if ( "$queue" != "" ) then
    grep -q "^#SBATCH -A.*:${queue}"'$' $script
    if ( $status != 0 ) sed -i -e "s|^#SBATCH -A.*|&:${queue}|" $script
  endif

  set ok = 0
  @ try = 1
  while ( $try <= $ntry_max )
    sbatch $sbatch_opts $script >& $sbout
    set st = $status
    if ( $st == 0 ) then
      grep -q '^Submitted batch job' $sbout
      if ( $status == 0 ) then
        set ok = 1
        break
      endif
    endif
    echo "   sbatch failed for $wdir (try $try of $ntry_max): `head -1 $sbout`"
    @ try = $try + 1
    if ( $try <= $ntry_max ) sleep 10
  end
  cd "$cdir"

  if ( $ok == 1 ) then
    # mark exactly this line (by number) as submitted
    sed -i -e "${ln}s/^/#"'!'"/" $list
    @ nsub = $nsub + 1
    set consec = 0
    echo "submitted $wdir : `head -1 $sbout`"
    sleep 1
  else
    echo "ERROR: $wdir could not be submitted; left in the list"
    @ nbad = $nbad + 1
    @ consec = $consec + 1
    set badlist = ( $badlist $wdir )
    if ( $consec >= $nfail_max ) then
      echo "Stopping: $consec failed entries in a row (submit limit or scheduler problem?)"
      break
    endif
  endif
end

set nleft = `grep -c '^[^#]' $list`
echo "---- submitted $nsub job(s); $nleft entr(y/ies) still unsubmitted in $list"
if ( $nbad > 0 ) echo "---- failed entries: $badlist"

cleanup:
/bin/rm -f $ltmp $sbout
exit 0
