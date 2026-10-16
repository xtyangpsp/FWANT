#!/bin/csh

set modeldir = result.1th

cat /dev/null > tradeoff.dat
foreach damp ( 2 4 8 12 16 ) 
	foreach smooth ( 2 4 8 12 16 )
		set sumfile = `/bin/ls -f $modeldir/"sum.damp"$damp".smot"$smooth.*.dat`
		echo $sumfile
		set nkaisquared = `cat $sumfile | grep kai^2/N | awk '{print $3}'`
		set var = `cat $sumfile | grep "data_var reduction =" | awk '{print $6}'`
		set modelnorm = `cat $sumfile | grep "model_norm" | awk '{print $4}'`
		echo $damp $smooth $nkaisquared $var $modelnorm >> tradeoff.dat
	end
end
