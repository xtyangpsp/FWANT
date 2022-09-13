#!/bin/csh
##################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO SPLIT EGF INTO POSITIVE AND NEGATIVE
# TIME LEGS
#
#
#  Y.Shen, 01-20-2011
#  Haiying Gao, 10-11-2016
#  Cong Li, 08/15/2017
#	Xiaotao Yang, 12/04/2017
#                       1. read netsta list directly to avoid the ls argument too long.
##################################################################################
set area = Alaska 
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area
set egfdir = $wkdir/stacked_xcorr

#set cnvegfdir = $wkdir/OBS_stacked_xcorr_cnv
#if ( ! ( -e $cnvegfdir ) ) mkdir $cnvegfdir

cd $egfdir

set tlen = 1000 # duration of one sided (split) egf
set taperwid = 0.005  # taper tlen*taperwid = 30 s

foreach network ( "`cat /Users/xiaotaoyang/Work/FWANT/Alaska/alaskanetstalist_AK.txt`" )

foreach stnpair ( `/bin/ls -d  $network* ` ) #.BIRD.N4.T57A
cd $stnpair
echo $stnpair

echo "echo off" > split.m

# stack of all EGFs
foreach egf ( `/bin/ls -f *sym.SAC` )
set egfid = `echo $egf | awk -F. '{print $1"."$2"."$3"."$4"."$5"."$6}'`
echo "cut 0 $tlen" >> split.m
echo "r $egf" >> split.m
echo "rmean;rtr;taper w $taperwid" >> split.m
set egf_pos = $egfid.all.P.SAC
echo "w $egf_pos" >> split.m

echo "cut -$tlen 0" >> split.m
echo "r $egf" >> split.m
echo "reverse" >> split.m
echo "ch b 0" >> split.m
echo "rmean;rtr;taper w $taperwid" >> split.m
set egf_neg = $egfid.all.N.SAC
echo "w $egf_neg" >> split.m
end # foreach egf

# monthly stacks
foreach egf ( `/bin/ls -f *sym.m*.SAC` )
set egfid = `echo $egf | awk -F. '{print $1"."$2"."$3"."$4"."$5"."$6"."$8}'`
echo "cut 0 $tlen" >> split.m
echo "r $egf" >> split.m
echo "rmean;rtr;taper w $taperwid" >> split.m
set egf_pos = $egfid.P.SAC
echo "w $egf_pos" >> split.m

echo "cut -$tlen 0" >> split.m
echo "r $egf" >> split.m
echo "reverse" >> split.m
echo "ch b 0" >> split.m
echo "rmean;rtr;taper w $taperwid" >> split.m
set egf_neg = $egfid.N.SAC
echo "w $egf_neg" >> split.m
end # foreach egf

echo "quit" >> split.m

sac <<EOF>> sac.output
m split.m
EOF

/bin/rm split.m sac.output
cd ..
end  # foreach stnpair
end

