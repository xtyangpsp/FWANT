#!/bin/bash

for f in 1 2 3 4 5 6 7 8
do
	ftag="f"$f
	sed -e 's/FTAG/'${ftag}'/g' submit_make_inv_Gd_TEMPLATE.sh > submit_make_inv_Gd_${ftag}.sh
done
