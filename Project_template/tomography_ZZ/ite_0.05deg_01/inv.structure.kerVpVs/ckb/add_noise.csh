#!/bin/csh 
# add random noise to traveltimes in inv_Gd_list

set Gd_list = inv_Gd_list.ori

#echo > tt_list
cat /dev/null > tt_list
foreach tt ( `cat $Gd_list | awk '{print $1}'` )
cat $Gd_list | grep $tt | awk '{print $1,$3}' >> tt_list
end

matlab -nodesktop -nosplash <<EOF
[ttid,tt] = textread('tt_list','%s%f');
nt=length(tt)
% standard deviation of travel times in second
ntstd=0.5;
noise=ntstd*randn(nt,1);
ttwn=tt+noise;

outfile = ['tt_with_noise_list'];
fidout = fopen(outfile,'w');
for i=1:nt
%ttid=['num.' num2str(i)];
fprintf(fidout,'%s %6.3f\n',char(ttid(i)),ttwn(i));
end

EOF

#@ nd = 0
#echo > inv_Gd_list
cat /dev/null > inv_Gd_list
foreach tt ( `cat $Gd_list | awk '{print $1}'` )
#@ nd = $nd + 1
#set ttid = "num"$nd
#set ttwithnoise = `cat tt_with_noise_list | grep $ttid | awk '{print $2}'`
set ttwithnoise = `cat tt_with_noise_list | grep $tt | awk '{print $2}'`
#echo $nd $ttid $ttwithnoise
#echo $tt
#cat $Gd_list | grep $tt | awk '{print $1, $2, $ttwithnoise, $4 ,$5,$6,$7,$8,$9,$10,$11,$12,$13}' >> inv_Gd_list
set c2 = `cat $Gd_list | grep $tt | awk '{print $2}'`
set c4 = `cat $Gd_list | grep $tt | awk '{print $4}'`
set c5 = `cat $Gd_list | grep $tt | awk '{print $5}'`
set c6 = `cat $Gd_list | grep $tt | awk '{print $6}'`
set c7 = `cat $Gd_list | grep $tt | awk '{print $7}'`
set c8 = `cat $Gd_list | grep $tt | awk '{print $8}'`
set c9 = `cat $Gd_list | grep $tt | awk '{print $9}'`
set c10 = `cat $Gd_list | grep $tt | awk '{print $10}'`
set c11 = `cat $Gd_list | grep $tt | awk '{print $11}'`
set c12 = `cat $Gd_list | grep $tt | awk '{print $12}'`
set c13 = `cat $Gd_list | grep $tt | awk '{print $13}'`
echo $tt $c2 $ttwithnoise $c4 $c5 $c6 $c7 $c8 $c9 $c10 $c11 $c12 $c13 >> inv_Gd_list 
#>> inv_Gd_list
end

