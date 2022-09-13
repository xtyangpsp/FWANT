#! /bin/csh
##########################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO FIND SAC FILES THAT DO NOT
# HAVE FULL-DAY LENGTH
#
# Usage:
#
# Y.Shen, 03-09-2011
# Haiying,10-11-2016
# Cong,   10-12-2016  add log file
#         07-12-2017  deleting the zero-value sac files
##########################################################################
#      SET AND MAKE DIRECTORY
##########################################################################
# set area = Virginia
set wkdir =  '/Users/xiaotaoyang/Work/FWANT/Alaska'
set sacdir = $wkdir/sac_data
foreach network ( "`cat /Users/xiaotaoyang/Work/FWANT/Alaska/alaskanetcodes.txt`" )
echo $network
cd $sacdir

foreach station ( `/bin/ls -d $network*` )
cd $station
if ( -e log.txt) rm log.txt
if ( !(-e Sampledel)) mkdir Sampledel
if ( !(-e Nptdel)) mkdir Nptdel
if ( !(-e Nulldel)) mkdir Nulldel
@ nfile = 0
foreach sacfile ( `/bin/ls -f *.SAC ` ) #
@ nfile = $nfile + 1
echo $nfile >>log.txt
set npts = `echo "r $sacfile; lh NPTS; quit" | sac | grep "NPTS =" | awk '{print $3}'`
set delta = `echo "r $sacfile; lh DELTA; quit" | sac | grep "DELTA =" | awk '{print $3}'`
set dmin = ` echo " r $sacfile; lh DEPMIN; quit" | sac | grep "DEPMIN =" | awk '{print $3}'`
set dmax = ` echo " r $sacfile; lh DEPMAX; quit" | sac | grep "DEPMAX =" | awk '{print $3}'`
set dmean = ` echo " r $sacfile; lh DEPMEN; quit" | sac | grep "DEPMEN =" | awk '{print $3}'`
echo $sacfile $npts $delta $dmin $dmax $dmean >>log.txt
############################################################################
#   CHECK THE NUMBER OF POINT AND THE SAMPLING INTERVAL, PROCESS THE NON FULL-DAY LENGTH
############################################################################
if ( ! ($delta == 5.000000e-01)  ) then
echo $sacfile $delta != 0.50 >>log.txt
mv $sacfile ./Sampledel
echo "del '$sacfile' for sample" >>log.txt
endif

if ( $npts <= 172800 & $npts > 172000 ) then # Notice: if the sampling ratio change the value here should be changed
echo "cuterr fillz; cut B 0 N 172801; r $sacfile; w $sacfile; quit" | sac
echo $sacfile $npts
set npts = `echo "r $sacfile; lh NPTS; quit" | sac | grep "NPTS =" | awk '{print $3}'`
endif

if ( $npts > 172800 ) then
echo "cuterr fillz; cut B 0 N 172801; r $sacfile; w $sacfile; quit" | sac
echo $sacfile $npts
set npts = `echo "r $sacfile; lh NPTS; quit" | sac | grep "NPTS =" | awk '{print $3}'`
endif


if ( ! ($npts == 172801)  ) then # delete the data whose length is less than 430000
echo $sacfile $npts != 172801 >>log.txt
mv $sacfile ./Nptdel
echo "del '$sacfile' for npts" >>log.txt
endif 

set dmax1 = `echo $dmax | awk 'function abs(value) { return (value <0 ? -value : value) } { print (abs($1) <= 1e-40 ? 1:0 ) }'`
set dmin1 = `echo $dmin | awk 'function abs(value) { return (value <0 ? -value : value) } {print (abs($1) <= 1e-40 ? 1:0 )}'`
set dmean1 = `echo $dmean | awk 'function abs(value) { return (value <0 ? -value : value) } {print (abs($1) <= 1e-40 ? 1:0 )}'`

if ( ($dmax1 == 1) & ($dmin1 == 1) & (($dmax == 'nan')||($dmin == 'nan')) )  then
echo $sacfile $dmax $dmin $dmean >> log.txt
mv sacfile ./Nulldel
echo "del '$sacfile' for nullvalue" >>log.txt
endif

end #foreach sacfile
cd ..
end #foreach station
end #foreach network
