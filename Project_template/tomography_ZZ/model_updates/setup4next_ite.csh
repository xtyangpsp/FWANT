#!/bin/csh
set wkdir = /depot/xtyang/data/projects/xtyang/RioGrande/FWANT_ZZ
set current_ite = ite_0.025deg_02
set next_ite = ite_0.025deg_03
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
rsync -avz bin $targetdir/
cp *.sh $targetdir
cp *.csh $targetdir
cp *.m $targetdir
cp -r skel $targetdir
/bin/rm $targetdir/skel/fx/input
ln -s $targetdir/../sim.input $targetdir/skel/fx/input
set current_SeisFD3D = $wkdir"/model_updates/SeisFD3D.conf_"$current_ite
set next_SeisFD3D = $targetdir/skel/fx/SeisFD3D.conf
/bin/cp $current_SeisFD3D $next_SeisFD3D
# link bin
cd skel/fx
ln -s ../../bin .
ln -s ../../../sim.input input
cd ../../
cd ..

# syn.seismograms
echo "setup syn.seismograms"
set subdir = syn.seismograms
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.m $targetdir
cp *.csh $targetdir
cp *.sh $targetdir
cp README $targetdir
#change iteration number.
sed -e 's/'${current_ite}'/'${next_ite}'/g' snap2sac_spher.m > ${targetdir}/snap2sac_spher.m

cd ..

# measure
echo "setup measure"
set subdir = measure
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.csh $targetdir
cp *.sh $targetdir
cp *.m $targetdir
cp *.py $targetdir
#change iteration number.
#sed -e 's/'${current_ite}'/'${next_ite}'/g' measure_phase_delay.m > ${targetdir}/measure_phase_delay.m
#sed -e 's/'${current_ite}'/'${next_ite}'/g' measure_phase_delay_velocityanalysis.m > ${targetdir}/measure_phase_delay_velocityanalysis.m
sed -e 's/'${current_ite}'/'${next_ite}'/g' measure_phase_delay_velocityanalysis_oneside.m > ${targetdir}/measure_phase_delay_velocityanalysis_oneside.m
sed -e 's/'${current_ite}'/'${next_ite}'/g' assemble_measure.csh > ${targetdir}/assemble_measure.csh
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
#cp *.py $targetdir
cp README $targetdir
cp butter.bp* $targetdir
set next_SeisFD3D = $targetdir/SeisFD3D.conf
/bin/cp $current_SeisFD3D $next_SeisFD3D
cp TomoKernel.conf $targetdir
#change iteration number.
sed -e 's/'${current_ite}'/'${next_ite}'/g' kern_synthetic4kernel.m > ${targetdir}/kern_synthetic4kernel.m
sed -e 's/'${current_ite}'/'${next_ite}'/g' kern_measure2kernelconf.sh > ${targetdir}/kern_measure2kernelconf.sh
cd $targetdir 
rsync -avz ${citedir}/sim.kernel/bin .
ln -s ../sim.input input
cd $citedir

echo "setup inv.structure"
set subdir = inv.structure.kerVpVs
cd $subdir
set targetdir = $wkdir/$next_ite/$subdir
cp *.m $targetdir
cp block*.dat $targetdir
cp block*.nc $targetdir
cp README $targetdir
rsync -avz Codes $targetdir/

set locmodeldir = FreqAll.Z.model
set targetmodeldir = $targetdir/FreqAll.Z.model
if ( ! ( -e $targetmodeldir ) ) mkdir $targetmodeldir
cd $locmodeldir
cp *.sh $targetmodeldir
cp *.csh $targetmodeldir
cp *.m $targetmodeldir
#change iteration number.
sed -e 's/'${current_ite}'/'${next_ite}'/g' make_inv_Gd_list.sh > ${targetmodeldir}/make_inv_Gd_list.sh
sed -e 's/'${current_ite}'/'${next_ite}'/g' mk_st_list.csh > ${targetmodeldir}/mk_st_list.csh

#needed here only when running inversion in subsampling bootstrapping method.
mkdir $targetmodeldir/subsamples
rsync -avz subsamples/*.m subsamples/*.sh subsamples/README $targetmodeldir/subsamples/

cp *_list $targetmodeldir
rm $targetmodeldir/inv_Gd_list $targetmodeldir/inv_ev_list $targetmodeldir/inv_st_list
cp *.dat $targetmodeldir
cd ..
cd ..
