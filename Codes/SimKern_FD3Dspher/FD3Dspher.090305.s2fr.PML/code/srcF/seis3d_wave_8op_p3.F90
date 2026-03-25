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
real(DP) :: dtmax_global
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

! 5-1B BBB
! 4-4A BBF
! 7-3B BFB
! 2-2A BFF

! 1-1A FFF
! 8-4B FFB
! 3-3A FBF
! 6-2B FBB
call check_my_stability(dtmax_global)
loop_time: do
!-----------------------------------------------------------------------------

if ( ntime>nt ) exit

! 5-1B BBB
! { ============================= third : lddrk4 =============================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 8-4B FFB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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

! 1-1A FFF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxB_LyB_LzB
call abs_LxB_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxF_LyF_LzF
call abs_LxF_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 4-4A BBF
! { ============================= forth : lddrk4 =============================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxF_LyF_LzB
call abs_LxF_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxB_LyB_LzF
call abs_LxB_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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

! 7-3B BFB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 6-2B FBB
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 3-3A FBF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxB_LyF_LzB
call abs_LxB_LyF_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxF_LyB_LzF
call abs_LxF_LyB_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

! 2-2A BFF
! {==========================================================================
! prepare
call swmpi_time_write(ntime,fnm_log)
call macdrp_syn
call abs_syn
! the 1th stage
call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
  call src_surface(ntime,0.0_SP,stept)
#endif
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
  call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
call macdrp_RK_beg(firRKa(1),firRKb(1))
call abs_RK_beg(firRKa(1),firRKb(1))
! the 2th stage
call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(1),stept)
#endif
call macdrp_LxF_LyB_LzB
call abs_LxF_LyB_LzB
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
call macdrp_RK_inn(firRKa(2),firRKb(2))
call abs_RK_inn(firRKa(2),firRKb(2))
! the 3th stage
call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(2),stept)
#endif
call macdrp_LxB_LyF_LzF
call abs_LxB_LyF_LzF
  call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
  call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
call macdrp_RK_inn(firRKa(3),firRKb(3))
call abs_RK_inn(firRKa(3),firRKb(3))
! the 4th stage
call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
  call src_surface(ntime,firRKa(3),stept)
#endif
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
call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
! ========================================================================== }

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
contains
!-----------------------------------------------------------------------!
subroutine check_my_stability(dt_out)
    ! Variables are split across para_mod and grid_mod
    use para_mod,  only: ni, nj, nk, ni1, nj1, nk1, stept
    use grid_mod,  only: x, y, z
    use media_mod, only: lambda, mu, rho
    use mpi_mod,   only: thisid ! ierr is handled locally or via mpi_mod depending on setup
    
    real(DP), intent(out) :: dt_out
    integer :: i, j, k, ii, jj, kk, local_ierr
    real(DP) :: Vp, dtLe, dtlocal, dtmax_local
    real(DP) :: p_a(3), p_b(3), p_c(3)

    dtmax_local = 1.0d10
    ! TEST: Print rank info once to ensure output is working
    if (thisid(1)==0 .and. thisid(2)==0 .and. thisid(3)==0) then
        write(*,*) "DEBUG: Starting stability loop for Rank 0..."
        call flush(6)
    end if

    do k = 1, nk
    do j = 1, nj
    do i = 1, ni
        ii = i + ni1 - 1
        jj = j + nj1 - 1
        kk = k + nk1 - 1
        ! Check if any input is NaN (x /= x is the standard NaN test in Fortran)
        if (rho(i,j,k) /= rho(i,j,k) .or. x(ii) /= x(ii)) then
             write(*,*) "RANK ", thisid, " FOUND NaN at i,j,k: ", i,j,k
             call flush(6)
             call MPI_ABORT(MPI_COMM_WORLD, 1, local_ierr)
        end if

        Vp = sqrt((lambda(i,j,k) + 2.0_DP*mu(i,j,k)) / rho(i,j,k))

        ! Forward neighbors
        p_a = (/ x(ii+1), y(jj),   z(kk)   /)
        p_b = (/ x(ii),   y(jj+1), z(kk)   /)
        p_c = (/ x(ii),   y(jj),   z(kk+1) /)
        dtLe = dist_p2p(x(ii), y(jj), z(kk), p_a, p_b, p_c)
        if (dtLe <= 0.0_DP) then
            print *, "ERROR: Zero grid spacing at indices: ", ii, jj, kk
            print *, "Coords: ", x(ii), x(ii+1)
        endif
        ! Backward neighbors
        p_a = (/ x(ii-1), y(jj),   z(kk)   /)
        p_b = (/ x(ii),   y(jj-1), z(kk)   /)
        p_c = (/ x(ii),   y(jj),   z(kk-1) /)
        dtLe = min(dtLe, dist_p2p(x(ii), y(jj), z(kk), p_a, p_b, p_c))

        dtlocal = 1.3_DP / Vp * dtLe
        ! Force an abort if we hit the zero
        if (dtlocal < 1.0d-18) then
            write(*,*) "--- ZERO DETECTED ON RANK ", thisid, " ---"
            write(*,*) "i,j,k: ", i, j, k
            write(*,*) "Vp: ", Vp, " dtLe: ", dtLe
            write(*,*) "Media: ", rho(i,j,k), lambda(i,j,k), mu(i,j,k)
            call flush(6)
            call MPI_ABORT(MPI_COMM_WORLD, 1, local_ierr)
        end if

        if (dtlocal < dtmax_local) dtmax_local = dtlocal
    end do
    end do
    end do

    ! Note: Use a local ierr variable for the MPI call
    call MPI_ALLREDUCE(dtmax_local, dt_out, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, local_ierr)

    if (thisid(1) == 0 .and. thisid(2) == 0 .and. thisid(3) == 0) then
        write(*,'(A)') "--------------------------------------------------------"
        write(*,'(A,F12.8)') "  CFL STABILITY ANALYSIS"
        write(*,'(A,F12.8)') "  Max allowable dt: ", dt_out
        write(*,'(A,F12.8)') "  Current stept:    ", stept
        
        if (real(stept, DP) > dt_out) then
            write(*,'(A)') "  >>> STATUS: UNSTABLE <<<"
            call MPI_ABORT(MPI_COMM_WORLD, 1, local_ierr)
        else
            write(*,'(A)') "  >>> STATUS: STABLE <<<"
            write(*,'(A)') "--------------------------------------------------------"
        end if
    end if
end subroutine check_my_stability

function dist_p2p(x0,y0,z0,p1,p2,p3) result(d)
    real(DP), intent(in) :: x0, y0, z0
    real(DP), intent(in) :: p1(3), p2(3), p3(3)
    real(DP) :: d, a, b, c, length
    a = (p2(2)-p1(2))*(p3(3)-p1(3)) - (p2(3)-p1(3))*(p3(2)-p1(2))
    b = (p2(3)-p1(3))*(p3(1)-p1(1)) - (p2(1)-p1(1))*(p3(3)-p1(3))
    c = (p2(1)-p1(1))*(p3(2)-p1(2)) - (p2(2)-p1(2))*(p3(1)-p1(1))
    length = sqrt(a**2 + b**2 + c**2)
    d = abs(a*(x0-p1(1)) + b*(y0-p1(2)) + c*(z0-p1(3))) / (length + 1.0d-20)
end function dist_p2p
end program seis3d_wave

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
