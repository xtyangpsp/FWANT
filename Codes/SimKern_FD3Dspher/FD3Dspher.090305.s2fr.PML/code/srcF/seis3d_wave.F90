program seis3d_wave

! This is the main program to simulate seismic wave propagation
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:59:46 -0500 (Fri, 16 Jan 2009) $
! $Revision: 70 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

!#define USEPML
#define PMLNT 1

!{ -- declare module used --
use mpi
use para_mod
use io_mod
use mpi_mod
use grid_mod
use media_mod
use src_mod
use abs_mod
use macdrp_mod
#ifdef WithOMP
use omp_lib
#endif
!} -- end declare module used ---

implicit none
integer ntime,ierr

call MPI_INIT(ierr)

call get_conf_name(fnm_conf)

call swmpi_init(fnm_conf)
call swmpi_cart_creat

call para_init(fnm_conf)
call swmpi_reinit_para

call grid_fnm_init(fnm_conf)
call grid_alloc
call grid_coord_import(thisid(1),thisid(2),thisid(3))
call grid_metric_import(thisid(1),thisid(2),thisid(3))

call media_fnm_init(fnm_conf)
call media_alloc
call media_import(thisid(1),thisid(2),thisid(3))

! --- C ---
! We pass the MPI coordinates (thisid) and grid/media arrays
call check_stability(thisid(1), thisid(2), thisid(3), &
                            ni, nj, nk, stept,               &
                            gx, gy, gz, lambda, mu, rho)
! --------------------------

call src_fnm_init(fnm_conf)
call src_import(thisid(1),thisid(2),thisid(3))
!call src_choose

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_snap_locate(thisid(1),thisid(2),thisid(3))
call io_pt_import
call io_seismo_init

!call grid_dealloc(iscoord=.true.)

call macdrp_init

call abs_init(fnm_conf)

call swmpi_datatype
call macdrp_mesg_init

ntime=0

call io_rest_import(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)

call swmpi_time_init(fnm_log,ntime)

#ifdef WithOMP
call OMP_SET_NUM_THREADS(2)
#endif


! 1-1A FFF
! 8-4B FFB
! 5-1B BBB
! 4-4A BBF

! 6-2B FBB
! 7-3B BFB
! 2-2A BFF
! 3-3A FBF

loop_time: do
!-----------------------------------------------------------------------------

if ( ntime>nt ) exit

#ifdef USEPML
if ( ntime< PMLNT .or. ntime >=2*PMLNT) then
#endif
! 1-1A FFF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

! 8-4B FFB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 5-1B BBB
! { ============================= third : lddrk4 =============================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

! 4-4A BBF
! { ============================= forth : lddrk4 =============================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

#ifdef USEPML
end if
if (ntime>=PMLNT) then

! 6-2B FBB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

! 7-3B BFB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

! 2-2A BFF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

! 3-3A FBF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
call macdrp_RK_fin(firRKb(4))
call abs_RK_fin(firRKb(4))
#ifdef WITHQS
 call atten_graves
#endif

! save result
     ntime=ntime+1
call macdrp_check(ntime)
call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
! ========================================================================== }

end if
#endif

!-----------------------------------------------------------------------------
end do loop_time

call io_seismo_close
call io_wave_close
call swmpi_time_end(fnm_log)

call macdrp_destroy
call grid_dealloc
call media_destroy
!call abs_destroy
call src_destroy

call MPI_FINALIZE(ierr)
!-----------------------------------------------------------------------!
!contains
!-----------------------------------------------------------------------!
subroutine check_stability(id1, id2, id3, ni, nj, nk, stept, gx, gy, gz, lambda, mu, rho)
    use constants_mod
    use math_mod, only : dist_point2plane
    use mpi
    implicit none

    integer, intent(in) :: id1, id2, id3, ni, nj, nk
    real(SP), intent(in) :: stept
    real(SP), dimension(0:ni+1,0:nj+1,0:nk+1), intent(in) :: gx, gy, gz, lambda, mu, rho

    integer  :: i, j, k, ierr
    real(SP) :: Vp, dtLe, dtlocal, dtmax_local, dtmax_global
    
    dtmax_local = 1.0e10

    ! Each MPI rank checks its own local chunk of the model
    do k = 1, nk
    do j = 1, nj
    do i = 1, ni
        ! P-wave velocity is the limiting factor for stability
        Vp = sqrt((lambda(i,j,k) + 2.0*mu(i,j,k)) / rho(i,j,k))

        ! Calculate the shortest distance to cell faces (essential for non-uniform grids)
        dtLe = min( &
            dist_point2plane((/gx(i,j,k),gy(i,j,k),gz(i,j,k)/),             &
               (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
               (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
               (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
            dist_point2plane((/gx(i,j,k),gy(i,j,k),gz(i,j,k)/),             &
               (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
               (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
               (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /))   &
        )

        ! 1.3 is the safety factor used in SeisFD3D for stability
        dtlocal = 1.3 / Vp * dtLe
        if (dtlocal < dtmax_local) dtmax_local = dtlocal
    end do
    end do
    end do

    ! Collective communication: all CPUs share their minimum to find the global minimum
    call MPI_ALLREDUCE(dtmax_local, dtmax_global, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Only the Master Rank (0,0,0) prints the results to the log
    if (id1 == 0 .and. id2 == 0 .and. id3 == 0) then
        write(*,*) "--------------------------------------------------------"
        write(*,*) "  MPI STABILITY CHECK"
        write(*,*) "  Max P-wave allowed dt: ", dtmax_global
        write(*,*) "  User configured stept: ", stept
        
        if (stept > dtmax_global) then
            write(*,*) "  ERROR: Simulation is UNSTABLE!"
            write(*,*) "  Please reduce 'stept' in SeisFD3D.conf to < ", dtmax_global
            write(*,*) "--------------------------------------------------------"
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        else
            write(*,*) "  Stability Check Passed. Starting Simulation..."
            write(*,*) "--------------------------------------------------------"
        end if
    end if
end subroutine check_stability
end program seis3d_wave

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
