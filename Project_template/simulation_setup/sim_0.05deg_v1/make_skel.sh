#!/bin/bash

hr="-"
while [ ${#hr} -lt 70 ]; do hr=${hr}"-"; done

echo $hr
# link bin to here
dir_bin='../../code/bin'
if [ ! -e ./bin ]; then
   if [ ! -d $dir_bin ]; then
      echo "$dir_bin doesn't exist"
      exit
   fi
   echo "ln $dir_bin ./bin"
   ln -s $dir_bin ./bin
else
   echo "./bin exists"
fi

echo $hr
# make directory skel
dir_list="input output checkpoint"
for nm in $dir_list; do
    if [ ! -d $nm ]; then
       echo "mkdir $nm"
       mkdir -p $nm
    else
       echo "$nm exists"
    fi
done

echo $hr
# generate checkpoint.dat
if [ ! -e checkpoint.dat ]; then
   echo "init checkpoint.dat"
   echo "0 0 0 # checkpoint, syncpoint, new nt" > checkpoint.dat
else
  echo "checkpoint.dat exists"
fi
