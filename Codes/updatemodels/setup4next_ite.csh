#!/bin/csh
set wkdir = /home/hg49a/pathtoxiaotao/FWANT/Alaska
#set wkdir = /Volumes/MacHD320/Work/UMass_PD/NewEnglandANT/FWANT/eUS
set current_ite = ite_0.05deg_01
set next_ite = ite_0.05deg_02
set citedir = $wkdir/$current_ite
set nitedir = $wkdir/$next_ite
if ( ! ( -e $nitedir ) ) then
mkdir $nitedir
endif

cd $citedir

echo "current iteration: " $citedir
echo "next iteration: " $nitedir
echo "setup subdirectories ..."

foreach subdir ( `/bin/ls -d *` )
set nextsubdir = $wkdir/$next_ite/$subdir
if ( ! ( -e $nextsubdir ) ) then
mkdir $nextsubdir
endif
end

# input model
set updatedmodeldir = $wkdir"/model_updates/updated_input_"$current_ite
set targetdir = $wkdir/$next_ite/sim.input
echo "cp " $updatedmodeldir " to " $targetdir
cp $updatedmodeldir/* $targetdir

# sim.station
echo "setup sim.station"
set subdir = sim.station
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.sh $targetdir
cp *.csh $targetdir
cp *.m $targetdir
cp -r skel $targetdir
/bin/rm $targetdir/skel/fx/input
#ln -s ../sim.input $targetdir/skel/fx/input
set current_SeisFD3D = $wkdir"/model_updates/SeisFD3D.conf_"$current_ite
set next_SeisFD3D = $targetdir/skel/fx/SeisFD3D.conf
/bin/cp $current_SeisFD3D $next_SeisFD3D
# link bin
cd skel/fx
ln -s ../../../code/bin .
cd ../../
ln -s ../../code/bin .
cd ..

# syn.seismograms
echo "setup syn.seismograms"
set subdir = syn.seismograms
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.m $targetdir
cp *.csh $targetdir
cp README $targetdir
cd ..

# measure
echo "setup measure"
set subdir = measure
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.csh $targetdir
cp *.m $targetdir
cp README $targetdir
cd ..

# sim.kernel
echo "setup sim.kernel"
set subdir = sim.kernel
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.m $targetdir
cp *.sh $targetdir
cp *.csh $targetdir
cp README $targetdir
cp butter.bp* $targetdir
set next_SeisFD3D = $targetdir/SeisFD3D.conf
/bin/cp $current_SeisFD3D $next_SeisFD3D
cp TomoKernel.conf $targetdir
cd $targetdir 
ln -s ../../codes/bin bin
ln -s ../sim.input input
cd $citedir

echo "setup inv.structure"
set subdir = inv.structure.kerVpVs
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.m $targetdir
cp block*.dat $targetdir
cp README $targetdir
set locmodeldir = Freq12345.Z.model
set targetmodeldir = $targetdir/Freq12345.Z.model
if ( ! ( -e $targetmodeldir ) ) mkdir $targetmodeldir
cd $locmodeldir
cp *.sh $targetmodeldir
cp *.csh $targetmodeldir
cp *.m $targetmodeldir
cp *_list $targetmodeldir
cp *.dat $targetmodeldir
cd ..
cd ..
