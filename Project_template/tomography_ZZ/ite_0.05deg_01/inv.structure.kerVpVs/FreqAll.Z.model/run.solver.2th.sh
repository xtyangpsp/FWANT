#!/bin/bash

outdir="result.2th"
fnm_smooth='../block.33x40x16.6x6x6.smooth2th.dat'

# input list file name
list_Gd="inv_Gd_list"

# damping list
#damp_list="3  5  7  9 10 15 20 30 40 50 60 70 80 90"
#damp_list="2 4  8 16 32 64 128 256"
#damp_list="2 4 8 12 16 24 32 48 64 128 256"

#smot_list="64 8 16 32 128"
#damp_list="32 4 8 16 64 128 256"

#smot_list="2 4 8 12 16 20 24 28 32 36 40 44 48 52 64 128 256"
#damp_list="2 4 8 12 16 20 24 28 32 36 40 44 48 52 64 128 256"

smot_list="2 8 12 16 20 24 28 32 36 48 64 128 256"
damp_list="8 16 32 36 48 64 128 256"

weigst_list="10"
weigev_list="50"
weiglo_list="10"

# station terms per data
nst=1
# event terms per data
nev=1

# num_cmp should be 2 when both Vp and Vs are inverted
num_cmp=2
ker_thres="1e-9"
data_thres=10

#SOLVER='../../code/inv.LSQR/solver_damp_smooth'
SOLVER='../../code/inv.LSQR/solver_damp_smo_loc'

if [ ! -d $outdir ]; then
   mkdir -p $outdir
fi

# ----------------------------------------------------------------------
function inv_compile {
  #FC=pgf90
  FC=ifort
  SRC='../../code/inv.LSQR'
  echo "$FC -O3 -o $SOLVER $SRC/mod_LSQR.f90 $SRC/solver_damp_smooth.f90"
  $FC -O3 -o $SOLVER $SRC/mod_LSQR.f90 $SRC/solver_damp_smooth.f90
}

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
#echo "compiling the program ..."
#time inv_compile;

echo ""
echo "solving with different daming ..."
time inv_loop;
