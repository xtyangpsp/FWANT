#!/bin/csh

cd measure_stnpair
set checkeddir = checked
if ( ! ( -e $checkeddir ) ) mkdir $checkeddir
set nodatadir = checked_no_data
if ( ! ( -e $nodatadir ) ) mkdir $nodatadir
foreach stnpairdir (`/bin/ls -l *.dat | grep "yang   0 Feb" | awk '{print $9}'`)
echo $stnpairdir
/bin/mv $stnpairdir $nodatadir
end

foreach stnpairdat (`/bin/ls *.dat`)
set stnpair = `echo $stnpairdat | awk -F. '{print $1"."$2"."$3"."$4}'`
set psname = $stnpair.eps

kghostview --portrait $psname

echo "Edit *.dat? y or n"
set vimflag = $<

if ( $vimflag == y ) then
vim $stnpairdat
endif

mv $stnpairdat $checkeddir

end

