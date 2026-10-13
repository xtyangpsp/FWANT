#!/bin/bash

outdir="result.1th"
fnm_smooth='../block.64x58x48.1x1x1.1x1x1.smooth1th.dat'

# input list file name
list_Gd="inv_Gd_list_dt6"

# damping and smoothing list

#smot_list="2 4 8 12 16"
#damp_list="2 4 8 12 16"
smot_list="8"
damp_list="8"

weigst_list="0"
weigev_list="0"
weiglo_list="0"

# station terms per data
nst=1
# event terms per data
nev=1
# location terms per data, dt in 3D (x-y-z)
nlo=3

# num_cmp should be 2 when both Vp and Vs are inverted
num_cmp=2
# see sim.kernel/run.kernel.assm.sh
ker_thres="5e-3"
data_thres=15

SOLVER=../Codes/InvCkb/inv.LSQR/solver_damp_smo_loc_yang

if [ ! -d $outdir ]; then
   mkdir -p $outdir
fi

# ----------------------------------------------------------------------
function inv_loop {
  for weig_lo in ${weiglo_list}; do
  for weig_st in ${weigst_list}; do
  for weig_ev in ${weigev_list}; do
  for smot in ${smot_list}; do
  for damp in ${damp_list}; do
  
  echo "loc st ev smooth damp=" $weig_lo $weig_st $weig_ev $smot $damp

$SOLVER << EOF
1
$list_Gd
$fnm_smooth
$num_cmp
$nst $weig_st
$nev $weig_ev
$nlo $weig_lo
$damp $smot
$ker_thres $data_thres
EOF

  mv try.xyz $outdir/try.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.sum $outdir/sum.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.err $outdir/err.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.sta $outdir/sta.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.evt $outdir/evt.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.loc $outdir/loc.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  
  done
  done
  done
  done
  done
}

# ----------------------------------------------------------------------

echo ""
echo "solving with different daming ..."
time inv_loop;
