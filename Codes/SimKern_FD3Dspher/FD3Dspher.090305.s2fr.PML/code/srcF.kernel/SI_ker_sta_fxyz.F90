program SI_ker_sta

! This program calculates finite-frequency sensitivity kernels for multi pairs
! of event and station by swapping event strain tensors when keeping station
! sgt in memory.
! The core of kernel calculation coming from the code of Z Li, C Po, Z Zhang and Y Shen.
!
! SI_ker_sta_fxyz.F90: modified copy of SI_ker_sta.F90 (same modules, same
! conf files; no new libraries).  Works with runs made from
!     fz only        (Rayleigh ZZ, pair lines as before),
!     fx + fy        (Love TT, pair lines with the two extra source weights),
!     or all three   (fz lines and fx+fy lines can be mixed; each pair line
!                     decides its own mode).
! Changes:
!  1) Love (TT) from two horizontal-force runs.  If a pair line of the
!     *_conf file carries two extra numbers after rot_mat, (sin a, cos a), the
!     event strain is wx*[x-force run] + wy*[y-force run] (source-side
!     rotation; run folders SGT_STAX / SGT_STAY under EVENT_ROOT/<source>).
!     Without them the program behaves as before (ZZ, RL, ...).
!     CONVENTION for TT pair lines: the three synthetic files listed after the
!     pair line (fnm_synx/y/z = BHN/BHE/BHZ slots) must hold the receiver x, y,
!     z velocity for a force along T_s ALREADY built from the two source runs
!     (fT.Vx, fT.Vy of snap2sac_sphere_fxfy.m), in the code frame (x = south,
!     y = east, z = up).  The program applies the receiver weights only; it does
!     not apply (sin a, cos a) to the synthetic.  Exactly one number after
!     rot_mat is rejected (give two or none).
!  2) One convolution per output channel: the receiver rotation (rot_mat) is
!     applied to the receiver strain BEFORE the convolution (the convolution is
!     linear) instead of convolving x,y,z and rotating the results.  Station
!     components with zero weight (z for T) are not read or convolved.
!  3) Convolution skips the leading zeros of both strain series and stops at
!     the last pick sample when all filters are causal (P1, N): exact.
!  4) OpenMP directives (ignored without -fopenmp), hoisted medium
!     coefficients in stress2strain, per-rank timing line KERNEL_PERF.
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2008 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 23:30:18 -0500 (Wed, 14 Jan 2009) $
! $Revision: 53 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#define VERBOSE
!#define DEBUG

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef KernelMPI
use mpi
#endif

implicit none

integer,parameter ::   &
  SIG_FILTER_P1  =100, &
  SIG_FILTER_P2  =120, &
  SIG_FILTER_NONE=130

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2

character (len=300) :: rcdstr
character (len=SEIS_STRLEN) ::           &
    fnm_main_conf,fnm_ker_conf,          &
    EVENT_ROOT,EVENT_OUTD,               &
    SGT_ROOT,SGT_STAX,SGT_STAY,SGT_STAZ, &
    evt_list,sgt_list,                   &
    fnm_synx,fnm_syny,fnm_synz,          &
    fnm_ker
character (len=SEIS_STRLEN) ::           &
    KSTANM,KEVNM
character (len=SEIS_STRLEN) :: filenm
character (len=15) :: idstr

character (len=SEIS_STRLEN),dimension(:),allocatable :: pnm_ker

integer :: staid,evtid

logical :: flag_fst
integer,dimension(SEIS_GEO) :: &
  bstart,bcount,bstride,       &
  sub_gs,sub_gc,sub_gt,        &
  subnc_s,subnc_c,subnc_t
integer :: sub_gtinv
integer,dimension(SEIS_GEO) :: blksiz,blkfst,blknum
integer :: kid,knt,nt1,nt2
real(SP) :: kdt,Tmin,Tmax

integer :: num_freq,num_pick
integer :: ncmp,npk,nfq
integer :: picknum(6)
integer,allocatable :: &
  picktidx(:,:),       &
  pickfreq(:),         &
  pickfilt(:)
real(SP),allocatable :: &
  picktwin(:,:)
real(SP),dimension(3,3) :: &
  rot_mat
character (len=2) :: filt_nm
real(SP) :: t1,t2

integer :: NOD,NOD1,NEL,NFACT
real(DP),dimension(:,:),allocatable :: filt_a,filt_b
real(DP),dimension(:),allocatable :: vecx,vecy,vecw

real(SP),dimension(:),pointer :: Vpr,Upr
real(SP),dimension(:),allocatable,target :: T,Vx,Ux,Vy,Uy,Vz,Uz
real(SP),dimension(:),allocatable :: VV,UU

real(SP) :: SGT_M0
real(SP),dimension(:,:,:,:),pointer :: Udapr,Udbpr
real(SP),dimension(:,:,:,:),allocatable :: &
    ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,         &
    ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,         &
    ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ
real(SP),dimension(:,:,:,:),allocatable,target :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS
real(SP),dimension(:,:,:,:),allocatable :: &
    Uda,Udb
logical :: flag_udx,flag_udy,flag_udz
logical :: flag_sgtx,flag_sgty,flag_sgtz
! channel results of the convolution (replace the in-place use of ExxS.. in
! the original) and work arrays
real(SP),dimension(:,:,:,:),allocatable,target :: Ca,Cb
real(SP),dimension(:,:,:,:),allocatable :: Etmp
real(SP),dimension(:),allocatable,target :: Vch,Uch
real(SP) :: wch(3),wsrc(2),rot_tmp(3,3)
integer :: pk_tmp(6),ic,nend,ierr2
logical :: need(3),flag_mix
! Cut the convolution at the last pick sample when all filters are causal
! (exact).  Set to .false. to reproduce the full-length original behaviour.
logical,parameter :: CONV_TRIM=.true.
! wall-clock per rank: load, convolution, filter+integrate+write
real(DP) :: perf(3)=0.0_DP,timer0,perf_start

real(SP),dimension(:),allocatable :: Ka,Kb
real(SP),dimension(:),allocatable :: E11R,E22R,E33R,E12R,E13R,E23R
real(SP),dimension(:),allocatable :: E11S,E22S,E33S,E12S,E13S,E23S
real(SP),dimension(:,:,:),allocatable :: Kap,Kaq,Kbp,Kbq

integer :: n_i,n_j,n_k
integer :: blki,blkj,blkk
integer :: i,ii,j,k,m,n,mt,ierr,i1,i2,j1,j2,k1,k2
integer :: p(1)
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt
!#ifdef KernelMPI
!integer,dimension(MPI_STATUS_SIZE) :: istatus
!#endif

integer :: ncidS,ncidF
integer :: TxxSid,TyySid,TzzSid,TxySid,TxzSid,TyzSid
integer :: TxxFid,TyyFid,TzzFid,TxyFid,TxzFid,TyzFid

integer,dimension(:),allocatable :: ncid
integer,dimension(:),allocatable :: kapid,kaqid,kbpid,kbqid

!----------------------------------------------------------------------

#ifdef KernelMPI
call MPI_INIT(ierr)
#endif
perf_start=kc_wtime()

call get_conf_name(fnm_conf)

! read kernel conf
fnm_ker_conf='TomoKernel.conf'
call init_kernel(fnm_ker_conf)
fnm_conf=fnm_main_conf

call para_init(fnm_conf)
call swmpi_init(fnm_conf)

#ifdef KernelMPI
call swmpi_cart_creat
call swmpi_reinit_para
call swmpi_datatype
  if (masternode) then
     read(*,"(a)") sgt_list
     write(*,*) sgt_list
  end if
  call MPI_BCAST(sgt_list,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
#else
call swmpi_set_gindx(0,0,0)
  read(*,"(a)") sgt_list
  write(*,*) sgt_list
#endif

call media_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

call init_filter(fnm_ker_conf)

!----------------------------------------------------------------------

#ifndef KernelMPI
print *, 'input mpi id:'
read *, n_i,n_j,n_k
call swmpi_change_fnm(n_i,n_j,n_k)
call swmpi_set_gindx(n_i,n_j,n_k)
thisid=(/ n_i,n_j,n_k /)
#else
n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#endif

call date_and_time(date=str_d0,time=str_t0)
write(idstr,"(a8,3(i2.2))") ' thisid=',n_i,n_j,n_k
#ifdef VERBOSE
  write(*,*) idstr, ' begins from ',str_d0,  &
                  ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
#endif

sub_gs=(bstart-1)*snap_subt(:,kid)+snap_subs(:,kid)
sub_gs(1)=out_i(sub_gs(1))
sub_gs(2)=out_j(sub_gs(2))
sub_gs(3)=out_k(sub_gs(3))
sub_gc=(snap_subc(:,kid)-bstart(:))/bstride(:)+1
do n=1,SEIS_GEO
   if (bcount(n)/=-1 .and. bcount(n)<sub_gc(n)) then
       sub_gc(n)=bcount(n)
   end if
end do
sub_gt=bstride*snap_subt(:,kid)
sub_gtinv=snap_tinv(kid)

call io_snap_locate(n_i,n_j,n_k)

!----------------------------------------------------------------------
kernode_if : if (snap_ishere(kid) .and. sgt_out(kid)) then

do i=1,sub_gc(1)
   i1=sub_gs(1)+(i-1)*sub_gt(1)
   if (i1>=snap_gsubs(1,kid)) exit
end do
do i=sub_gc(1),1,-1
   i2=sub_gs(1)+(i-1)*sub_gt(1)
   if (i2<=snap_gsube(1,kid)) exit
end do
subnc_s(1)=(i1-snap_gsubs(1,kid))/snap_subt(1,kid)+1
subnc_c(1)=(i2-i1)/sub_gt(1)+1

do j=1,sub_gc(2)
   j1=sub_gs(2)+(j-1)*sub_gt(2)
   if (j1>=snap_gsubs(2,kid)) exit
end do
do j=sub_gc(2),1,-1
   j2=sub_gs(2)+(j-1)*sub_gt(2)
   if (j2<=snap_gsube(2,kid)) exit
end do
subnc_s(2)=(j1-snap_gsubs(2,kid))/snap_subt(2,kid)+1
subnc_c(2)=(j2-j1)/sub_gt(2)+1

do k=1,sub_gc(3)
   k1=sub_gs(3)+(k-1)*sub_gt(3)
   if (k1>=snap_gsubs(3,kid)) exit
end do
do k=sub_gc(3),1,-1
   k2=sub_gs(3)+(k-1)*sub_gt(3)
   if (k2<=snap_gsube(3,kid)) exit
end do
subnc_s(3)=(k1-snap_gsubs(3,kid))/snap_subt(3,kid)+1
subnc_c(3)=(k2-k1)/sub_gt(3)+1

subnc_t=bstride

where (blksiz==-1 .or. blksiz>subnc_c)
   blksiz=subnc_c
end where
blknum=(subnc_c+blksiz-1)/blksiz

#ifdef DEBUG
write(100+myid,*) "thisid="
write(100+myid,*) thisid
write(100+myid,*) "sub_gs="
write(100+myid,*) sub_gs
write(100+myid,*) "sub_gc="
write(100+myid,*) sub_gc
write(100+myid,*) "sub_gt="
write(100+myid,*) sub_gt
write(100+myid,*) "sub_gtinv="
write(100+myid,*) sub_gtinv
write(100+myid,*) "subnc_s="
write(100+myid,*) subnc_s
write(100+myid,*) "subnc_c="
write(100+myid,*) subnc_c
write(100+myid,*) "subnc_t="
write(100+myid,*) subnc_t

write(100+myid,*) "blksiz="
write(100+myid,*) blksiz
write(100+myid,*) "blknum="
write(100+myid,*) blknum
#endif

call alloc_ker_var(blksiz(1),blksiz(2),blksiz(3),knt)
call alloc_media_local(blksiz(1),blksiz(2),blksiz(3))
call alloc_pick(num_pick)
call alloc_seismo_var(knt)

kdt=sub_gtinv*stept
T=(/(ii, ii=1,knt,1)/)*kdt
Tmin=kdt
Tmax=knt*kdt

staid=10001; evtid=10002
! loop all stations
open(staid,file=trim(sgt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("sgt_list open err:"//trim(sgt_list))

!-----------------------------------------------------------------------------
sta_loop: do

  read(staid,*,iostat=ierr) KSTANM,evt_list
  if (ierr<0) exit sta_loop
  if (KSTANM(1:1)=='#') cycle

#ifdef VERBOSE
  call date_and_time(date=str_d1,time=str_t1)
  write(*,*) idstr//trim(KSTANM)//" "//trim(evt_list), &
                    ': begins from ',str_d1,  &
                  ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

! loop each block
!-----------------------------------------------------------------------------
do blkk=blkfst(3),blknum(3)
do blkj=blkfst(2),blknum(2)
do blki=blkfst(1),blknum(1)

#ifdef VERBOSE
  call date_and_time(date=str_d1,time=str_t1)
   write(*,"(a7,3(i4,a1,i4.4),a9,3(i2.2))")   &
       ' block:',                             &
       blki,'/',blknum(1),                    &
       blkj,'/',blknum(2),                    &
       blkk,'/',blknum(3),                    &
       ', thisid=',n_i,n_j,n_k
  write(*,*) idstr, ': current time ',str_d1, &
           ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

   flag_sgtx=.false.; flag_sgty=.false.; flag_sgtz=.false.

   bsubs(1)=(blki-1)*blksiz(1)*subnc_t(1)+subnc_s(1)
   bsubs(2)=(blkj-1)*blksiz(2)*subnc_t(2)+subnc_s(2)
   bsubs(3)=(blkk-1)*blksiz(3)*subnc_t(3)+subnc_s(3)

   bsubc(1)=min(blksiz(1),subnc_c(1)-(bsubs(1)-subnc_s(1))/subnc_t(1))
   bsubc(2)=min(blksiz(2),subnc_c(2)-(bsubs(2)-subnc_s(2))/subnc_t(2))
   bsubc(3)=min(blksiz(3),subnc_c(3)-(bsubs(3)-subnc_s(3))/subnc_t(3))

   bsubt=subnc_t

   subc=bsubc;
   subt=snap_subt(:,kid)*bsubt
   subs=snap_subs(:,kid)+(bsubs-1)*snap_subt(:,kid)

#ifdef DEBUG
write(100+myid,*) "-------------------"
write(100+myid,*) "bli,j,k=",blki,blkj,blkk
write(100+myid,*) "bsubs="
write(100+myid,*) bsubs
write(100+myid,*) "bsubc="
write(100+myid,*) bsubc
write(100+myid,*) "bsubt="
write(100+myid,*) bsubt
write(100+myid,*) "subs="
write(100+myid,*) subs
write(100+myid,*) "subc="
write(100+myid,*) subc
write(100+myid,*) "subt="
write(100+myid,*) subt
#endif

if ( product(bsubc)*knt/=size(ExxS) ) then
   call realloc_ker_var(bsubc(1),bsubc(2),bsubc(3),knt)
   call alloc_media_local(bsubc(1),bsubc(2),bsubc(3))
end if

! load media
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt) 

! loop all event
!-----------------------------------------------------------------------------
open(evtid,file=trim(evt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("evt_list open err:"//trim(evt_list))

evt_loop: do
  read(evtid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  if (len_trim(rcdstr)==300) then
     print *, len_trim(rcdstr)
     print *, rcdstr
     call error_except("beyond maximum character 300 limit")
  end if
  rcdstr=adjustl(rcdstr)
  if (rcdstr(1:1)=='#') cycle

#ifdef VERBOSE
   write(*,*) idstr//trim(rcdstr)
#endif

  n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
  KEVNM=rcdstr(1:n)
  rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  read(rcdstr,*) picknum, rot_mat
  ! optional trailing pair (sin a, cos a): weights of the x- and y-force EVENT
  ! runs that build a source force along T (Love, TT).  Absent for ZZ.
  wsrc=0.0_SP
  read(rcdstr,*,iostat=ierr) pk_tmp, rot_tmp, wsrc
  flag_mix=(ierr==0)
  if (.not. flag_mix) then
     ! exactly one number after rot_mat is a typo, not a ZZ line: stop instead of silently running ZZ
     read(rcdstr,*,iostat=ierr2) pk_tmp, rot_tmp, wsrc(1)
     if (ierr2==0) call error_except('one source weight after rot_mat; give two (sin a, cos a) or none')
  end if
  if (flag_mix .and. all(wsrc==0.0_SP)) call error_except('source weights are both zero')
  if (flag_mix) then
     if (abs(wsrc(1)**2+wsrc(2)**2-1.0_SP)>1.0e-3_SP) &
        call warning_print('source weights (sin a, cos a) are not a unit vector')
  end if

  if (maxval(picknum)>num_pick) then
     num_pick=maxval(picknum)
     call alloc_pick(num_pick)
  end if

  ! components needed in this record = components with non-zero weight in any
  ! requested channel (ZZ: z; TT: x and y, because the z weight of T is 0)
  need=.false.
  do ic=1,6
     if (picknum(ic)<=0) cycle
     call channel_weights(ic,wch)
     need=need .or. (wch/=0.0_SP)
  end do
  flag_udx=need(1); flag_udy=need(2); flag_udz=need(3)
  if (any(picknum>0) .and. .not. any(need)) call error_except('rot_mat column is zero')

  read(evtid,"(a132)") fnm_synx
  read(evtid,"(a132)") fnm_syny
  read(evtid,"(a132)") fnm_synz
  if (flag_udx) call init_synthetic(fnm_synx,T,Vx,Ux)
  if (flag_udy) call init_synthetic(fnm_syny,T,Vy,Uy)
  if (flag_udz) call init_synthetic(fnm_synz,T,Vz,Uz)

  if (all(picknum<=0)) cycle

  if (abs(T(2)-T(1)-kdt)>SEIS_EQUAL) then
     print *, kdt,T(2),T(1)
     print *, T(2)-T(1),abs(T(2)-T(1))
     !call error_except('dt /= stept*tinv')
  end if

timer0=kc_wtime()
! ---- load station SGT: only the components that carry weight in this record.
!      ZZ: z only.  TT (rotated pick on T): x and y only, z is never read.
if ( flag_udx .and. (.not. flag_sgtx) ) then
   if ((.not. allocated(ExxX)) .or. (product(bsubc)*knt/=size(ExxX))) then
      call alloc_sgtx(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
   call load_sgt(trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAX), &
        ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,1.0_SP,.false.)
   flag_sgtx=.true.
   call stress2strain(ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,lambda,mu,bsubc)
end if
if ( flag_udy .and. (.not. flag_sgty) ) then
   if ((.not. allocated(ExxY)) .or. (product(bsubc)*knt/=size(ExxY))) then
      call alloc_sgty(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
   call load_sgt(trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAY), &
        ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,1.0_SP,.false.)
   flag_sgty=.true.
   call stress2strain(ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,lambda,mu,bsubc)
end if
if ( flag_udz .and. (.not. flag_sgtz) ) then
   if ((.not. allocated(ExxZ)) .or. (product(bsubc)*knt/=size(ExxZ))) then
      call alloc_sgtz(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
   call load_sgt(trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAZ), &
        ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ,1.0_SP,.false.)
   flag_sgtz=.true.
   call stress2strain(ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ,lambda,mu,bsubc)
end if
! ---- load event SGT
if (flag_mix) then
   ! TT: source force along T = wsrc(1)*[x-force run] + wsrc(2)*[y-force run].
   ! The two event runs live under the same sub-directory names as the
   ! station x/y runs (SGT_STAX, SGT_STAY); EVENT_OUTD is not used.
   call load_sgt(trim(EVENT_ROOT)//'/'//trim(KEVNM)//'/'//trim(SGT_STAX), &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,wsrc(1),.false.)
   call load_sgt(trim(EVENT_ROOT)//'/'//trim(KEVNM)//'/'//trim(SGT_STAY), &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,wsrc(2),.true.)
else
   call load_sgt(trim(EVENT_ROOT)//'/'//trim(KEVNM)//'/'//trim(EVENT_OUTD), &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,1.0_SP,.false.)
end if
call stress2strain(ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,lambda,mu,bsubc)
perf(1)=perf(1)+kc_wtime()-timer0

! ---- loop over output channels.  A channel is one unrotated component
!      (ncmp=1..3) or one rotated component (ncmp=4..6, weights = column
!      ncmp-3 of rot_mat).  Each channel needs ONE convolution, because the
!      convolution is linear in the receiver strain.
do ncmp=1,6

   if (picknum(ncmp)<=0) cycle
   call channel_weights(ncmp,wch)

do npk=1,picknum(ncmp)
   do
   read(evtid,'(a300)',iostat=ierr) rcdstr
   if (ierr<0) call error_except("read event pick error")
   if (len_trim(rcdstr)==300) then
      print *, len_trim(rcdstr)
      print *, rcdstr
      call error_except("beyond maximum character 300 limit")
   end if
   rcdstr=adjustl(rcdstr)
   if (rcdstr(1:1)=='#') cycle
   exit
   end do

   !read(evtid,*) pickfreq(npk),filt_nm,t1,t2,pnm_ker(npk)
   ! pickfreq
#ifdef VERBOSE
   write(*,*) idstr//trim(rcdstr)
#endif
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) pickfreq(npk)
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   if (pickfreq(npk)>num_freq) call error_except('frequency band id overflow')
   ! pickfilt
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   filt_nm=rcdstr(1:n)
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   select case (trim(filt_nm))
   case ('P1')
     pickfilt(npk)=SIG_FILTER_P1
   case ('P2')
     pickfilt(npk)=SIG_FILTER_P2
   case ('N')
     pickfilt(npk)=SIG_FILTER_NONE
   end select
   ! t1 t2
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t1
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t2
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   if (t1>Tmax) then
      write(*,*) "t1,Tmax=",t1,Tmax
      call warning_print("t1>Tmax, skip this rcd "//trim(rcdstr))
      cycle
   elseif (t1<=Tmin) then
      write(*,*) "t1,Tmin=",t1,Tmin
      call warning_print("t1<Tmin, reset t1 to Tmin"//trim(rcdstr))
      picktidx(1,npk)=1
      picktwin(1,npk)=T(1)
   else
      p=maxloc(T,T<=t1); picktidx(1,npk)=p(1); picktwin(1,npk)=T(p(1));
   end if
   if (t2<Tmin) then
      write(*,*) "t2,Tmin=",t2,Tmin
      call warning_print("t2<Tmin, skip this rcd"//trim(rcdstr))
      cycle
   elseif (t2>=Tmax) then
      write(*,*) "t2,Tmax=",t2,Tmax
      call warning_print("t2>Tmax, reset t2 to Tmax"//trim(rcdstr))
      picktidx(2,npk)=knt
      picktwin(2,npk)=T(knt)
   else
      p=minloc(T,T>=t2); picktidx(2,npk)=p(1); picktwin(2,npk)=T(p(1));
   end if
   ! pnm_ker
   pnm_ker(npk)=trim(rcdstr)

   ! kernel nc file
   if (flag_fst .and. (blki==1 .and. blkj==1 .and. blkk==1)) then
      fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
      call nfseis_grid3d_def(trim(pnm_ker(npk))//"/"//trim(fnm_ker),   &
          snap_subc(1,kid),snap_subc(2,kid), snap_subc(3,kid),ncid(npk), &
          "Finite frequency kernel")
      call nfseis_grid3d_defvar(ncid(npk),'phase_Vp'    ,kapid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'amplitude_Vp',kaqid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'phase_Vs'    ,kbpid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'amplitude_Vs',kbqid(npk))
      call nfseis_grid3d_enddef(ncid(npk))
   else
      fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
      call nfseis_open(trim(pnm_ker(npk))//"/"//trim(fnm_ker),ncid(npk))
      call nfseis_inq_varid(ncid(npk),'phase_Vp'    ,kapid(npk))
      call nfseis_inq_varid(ncid(npk),'amplitude_Vp',kaqid(npk))
      call nfseis_inq_varid(ncid(npk),'phase_Vs'    ,kbpid(npk))
      call nfseis_inq_varid(ncid(npk),'amplitude_Vs',kbqid(npk))
   end if
end do

! ---- synthetic of this channel (same arithmetic as the old syn_rotate)
   Vch=0.0_SP; Uch=0.0_SP
   if (wch(1)/=0.0_SP) then; Vch=Vch+wch(1)*Vx; Uch=Uch+wch(1)*Ux; end if
   if (wch(2)/=0.0_SP) then; Vch=Vch+wch(2)*Vy; Uch=Uch+wch(2)*Uy; end if
   if (wch(3)/=0.0_SP) then; Vch=Vch+wch(3)*Vz; Uch=Uch+wch(3)*Uz; end if
   Vpr=>Vch; Upr=>Uch
   Udapr=>Ca; Udbpr=>Cb

! ---- last sample the picks of this channel can need.  The filters are
!      causal for 'P1' and 'N' (output n depends on input <= n), so nothing
!      after the last pick window is needed and the convolution is cut there.
!      'P2' (zero phase) uses the whole trace, so it keeps the full length.
   nend=knt
   if (CONV_TRIM) then
      nend=maxval(picktidx(2,1:picknum(ncmp)))
      if (any(pickfilt(1:picknum(ncmp))/=SIG_FILTER_P1 .and. &
              pickfilt(1:picknum(ncmp))/=SIG_FILTER_NONE)) nend=knt
      nend=max(1,min(nend,knt))
   end if

   timer0=kc_wtime()
   call sgt2scatter_ch(wch,nend,kdt)
   perf(2)=perf(2)+kc_wtime()-timer0
   timer0=kc_wtime()
! cal kernel
do npk=1,picknum(ncmp)
   nfq=pickfreq(npk)

   ! filtering when different frq or filt
   if (npk==1  &
       .or. (nfq/=pickfreq(npk-1)) .or. (pickfilt(npk)/=pickfilt(npk-1)) ) then
   select case (pickfilt(npk))
   case (SIG_FILTER_P1)
     call scatter_filter(Udapr,Uda,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call scatter_filter(Udbpr,Udb,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call syn_filter(Vpr,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
     call syn_filter(Upr,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   case (SIG_FILTER_P2)
     call scatter_filtfilt(Udapr,Uda,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call scatter_filtfilt(Udbpr,Udb,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call syn_filtfilt(Vpr,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
     call syn_filtfilt(Upr,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   case (SIG_FILTER_NONE)
     Uda=Udapr;Udb=Udbpr;
     VV=Vpr; UU=Upr;
   case default
     print *, "npk,nfq,pickfilt=",npk,nfq,pickfilt(npk)
     call error_except('pickfilt error')
   end select
   end if

#ifdef VERBOSE
#ifdef KernelMPI
if (myid==0) &
#endif
   call export_info(pnm_ker(npk),T,VV,UU, &
        picktwin(:,npk),picktidx(:,npk),  &
        filt_a(:,nfq),filt_b(:,nfq))
#ifdef DEBUG
#ifdef KernelMPI
if (myid==0) &
#endif
   call export_debug(pnm_ker(npk),T,Vpr,Upr)
#endif
#endif

   nt1=picktidx(1,npk);nt2=picktidx(2,npk)
   call scatter2kernel(VV,UU,Uda,Udb,nt1,nt2,kdt,bsubc,Kap,Kaq,Kbp,Kbq)

   ! put
   call nfseis_put(ncid(npk),kapid(npk),Kap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kaqid(npk),Kaq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kbpid(npk),Kbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kbqid(npk),Kbq,bsubs,bsubc,bsubt)
end do ! npk

! close nc file
do npk=1,picknum(ncmp)
   call nfseis_close(ncid(npk))
end do ! npk

perf(3)=perf(3)+kc_wtime()-timer0
end do ! ncmp

end do evt_loop
close(evtid)

end do ! blki
end do ! blkj
end do ! blkk

end do sta_loop
close(staid)

!-----------------------------------------------------------------------------
end if kernode_if

#ifdef VERBOSE
  call date_and_time(date=str_d2,time=str_t2)
  write(*,*) idstr, ' finish at ',str_d2,  &
                     ',',str_t2(1:2),'/',str_t2(3:4),'/',str_t2(5:10)
#endif

!-----------------------------------------------------------------------------
#ifdef VERBOSE
write(*,'(a,3(i2.2),a,4(es12.4,1x))') 'KERNEL_PERF rank=',n_i,n_j,n_k, &
  ' load conv filter+kernel+write total(s)=',perf,kc_wtime()-perf_start
#endif
call media_destroy
call dealloc_all

#ifdef KernelMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine alloc_filt_var(nf,NOD1,NEL)
  integer,intent(in) :: nf,NOD1,NEL
  allocate(filt_b(NOD1,nf)); filt_b=0.0_DP
  allocate(filt_a(NOD1,nf)); filt_a=0.0_DP
  allocate(vecx(NEL)); vecx=0.0_DP
  allocate(vecy(NEL)); vecy=0.0_DP
  allocate(vecw(NEL)); vecw=0.0_DP
end subroutine alloc_filt_var

subroutine alloc_pick(npk)
  integer,intent(in) :: npk
  if (allocated(pnm_ker)) deallocate(pnm_ker)
  if (allocated(pickfreq)) deallocate(pickfreq)
  if (allocated(pickfilt)) deallocate(pickfilt)
  if (allocated(picktidx)) deallocate(picktidx)
  if (allocated(picktwin)) deallocate(picktwin)
  if (allocated(ncid ))    deallocate(ncid )
  if (allocated(kapid))    deallocate(kapid)
  if (allocated(kaqid))    deallocate(kaqid)
  if (allocated(kbpid))    deallocate(kbpid)
  if (allocated(kbqid))    deallocate(kbqid)

  allocate(pnm_ker(npk)); 
  allocate(pickfreq(npk)); pickfreq=0
  allocate(pickfilt(npk)); pickfilt=0
  allocate(picktidx(2,npk)); picktidx=0
  allocate(picktwin(2,npk)); picktwin=0.0
  allocate(ncid (npk)); ncid=0
  allocate(kapid(npk)); kapid=0
  allocate(kaqid(npk)); kaqid=0
  allocate(kbpid(npk)); kbpid=0
  allocate(kbqid(npk)); kbqid=0
end subroutine alloc_pick

subroutine alloc_seismo_var(nt)
  integer,intent(in) :: nt
  allocate(T(nt)); T=0.0
  allocate(VV(nt)); VV=0.0
  allocate(UU(nt)); UU=0.0
  allocate(Vx(nt)); Vx=0.0
  allocate(Ux(nt)); Ux=0.0
  allocate(Vy(nt)); Vy=0.0
  allocate(Uy(nt)); Uy=0.0
  allocate(Vz(nt)); Vz=0.0
  allocate(Uz(nt)); Uz=0.0
  allocate(Vch(nt)); Vch=0.0
  allocate(Uch(nt)); Uch=0.0
end subroutine alloc_seismo_var

subroutine alloc_ker_var(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0

  allocate(Uda(ki,kj,kk,kt)); Uda=0.0
  allocate(Udb(ki,kj,kk,kt)); Udb=0.0
  allocate(Ca(ki,kj,kk,kt)); Ca=0.0
  allocate(Cb(ki,kj,kk,kt)); Cb=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0

  allocate(Ka(kt),Kb(kt));Ka=0.0;Kb=0.0
  allocate(E11S(kt)); E11S=0.0
  allocate(E22S(kt)); E22S=0.0
  allocate(E33S(kt)); E33S=0.0
  allocate(E12S(kt)); E12S=0.0
  allocate(E13S(kt)); E13S=0.0
  allocate(E23S(kt)); E23S=0.0
  allocate(E11R(kt)); E11R=0.0
  allocate(E22R(kt)); E22R=0.0
  allocate(E33R(kt)); E33R=0.0
  allocate(E12R(kt)); E12R=0.0
  allocate(E13R(kt)); E13R=0.0
  allocate(E23R(kt)); E23R=0.0
end subroutine alloc_ker_var

subroutine realloc_ker_var(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' realloc_ker_var'
#endif
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
  if (allocated(Uda)) deallocate(Uda)
  if (allocated(Udb)) deallocate(Udb)
  if (allocated(Ca)) deallocate(Ca)
  if (allocated(Cb)) deallocate(Cb)
  if (allocated(Etmp)) deallocate(Etmp)

  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0

  allocate(Uda(ki,kj,kk,kt)); Uda=0.0
  allocate(Udb(ki,kj,kk,kt)); Udb=0.0
  allocate(Ca(ki,kj,kk,kt)); Ca=0.0
  allocate(Cb(ki,kj,kk,kt)); Cb=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0 
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0 
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0 
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0 
end subroutine realloc_ker_var

subroutine alloc_sgtx(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgtx'
#endif
  if (allocated(ExxX)) deallocate(ExxX)
  if (allocated(EyyX)) deallocate(EyyX)
  if (allocated(EzzX)) deallocate(EzzX)
  if (allocated(ExyX)) deallocate(ExyX)
  if (allocated(ExzX)) deallocate(ExzX)
  if (allocated(EyzX)) deallocate(EyzX)
  allocate(ExxX(ki,kj,kk,kt)); ExxX=0.0
  allocate(EyyX(ki,kj,kk,kt)); EyyX=0.0
  allocate(EzzX(ki,kj,kk,kt)); EzzX=0.0
  allocate(ExyX(ki,kj,kk,kt)); ExyX=0.0
  allocate(ExzX(ki,kj,kk,kt)); ExzX=0.0
  allocate(EyzX(ki,kj,kk,kt)); EyzX=0.0
end subroutine alloc_sgtx
subroutine alloc_sgty(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgty'
#endif
  if (allocated(ExxY)) deallocate(ExxY)
  if (allocated(EyyY)) deallocate(EyyY)
  if (allocated(EzzY)) deallocate(EzzY)
  if (allocated(ExyY)) deallocate(ExyY)
  if (allocated(ExzY)) deallocate(ExzY)
  if (allocated(EyzY)) deallocate(EyzY)
  allocate(ExxY(ki,kj,kk,kt)); ExxY=0.0
  allocate(EyyY(ki,kj,kk,kt)); EyyY=0.0
  allocate(EzzY(ki,kj,kk,kt)); EzzY=0.0
  allocate(ExyY(ki,kj,kk,kt)); ExyY=0.0
  allocate(ExzY(ki,kj,kk,kt)); ExzY=0.0
  allocate(EyzY(ki,kj,kk,kt)); EyzY=0.0
end subroutine alloc_sgty
subroutine alloc_sgtz(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgtz'
#endif
  if (allocated(ExxZ)) deallocate(ExxZ)
  if (allocated(EyyZ)) deallocate(EyyZ)
  if (allocated(EzzZ)) deallocate(EzzZ)
  if (allocated(ExyZ)) deallocate(ExyZ)
  if (allocated(ExzZ)) deallocate(ExzZ)
  if (allocated(EyzZ)) deallocate(EyzZ)
  allocate(ExxZ(ki,kj,kk,kt)); ExxZ=0.0
  allocate(EyyZ(ki,kj,kk,kt)); EyyZ=0.0
  allocate(EzzZ(ki,kj,kk,kt)); EzzZ=0.0
  allocate(ExyZ(ki,kj,kk,kt)); ExyZ=0.0
  allocate(ExzZ(ki,kj,kk,kt)); ExzZ=0.0
  allocate(EyzZ(ki,kj,kk,kt)); EyzZ=0.0
end subroutine alloc_sgtz

subroutine alloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  if (allocated(lambda)) deallocate(lambda)
  if (allocated(mu)) deallocate(mu)
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine alloc_media_local

subroutine dealloc_all
  if (allocated(Vx)) deallocate(Vx)
  if (allocated(Ux)) deallocate(Ux)
  if (allocated(T)) deallocate(T)
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxX)) deallocate(ExxX)
  if (allocated(EyyX)) deallocate(EyyX)
  if (allocated(EzzX)) deallocate(EzzX)
  if (allocated(ExyX)) deallocate(ExyX)
  if (allocated(ExzX)) deallocate(ExzX)
  if (allocated(EyzX)) deallocate(EyzX)
  if (allocated(Ca)) deallocate(Ca)
  if (allocated(Cb)) deallocate(Cb)
  if (allocated(Etmp)) deallocate(Etmp)
  if (allocated(Vch)) deallocate(Vch)
  if (allocated(Uch)) deallocate(Uch)
  if (allocated(Ka)) deallocate(Ka)
  if (allocated(Kb)) deallocate(Kb)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
end subroutine dealloc_all

!-----------------------------------------------------------------------------
subroutine init_kernel(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid,n

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)
call string_conf(fid,1,'snap_id',2,kid)
call string_conf(fid,1,'first_run',2,flag_fst)
do n=1,SEIS_GEO
 call string_conf(fid,1,'sub_start',n+1,bstart(n))
 call string_conf(fid,1,'sub_count',n+1,bcount(n))
 call string_conf(fid,1,'sub_stride',n+1,bstride(n))
end do
do n=1,SEIS_GEO
 call string_conf(fid,1,'block_from',n+1,blkfst(n))
 call string_conf(fid,1,'block_size',n+1,blksiz(n))
end do

call string_conf(fid,1,"SGT_TIMEDIM_SIZE",2,knt)
call string_conf(fid,1,"SGT_M0",2,SGT_M0)
call string_conf(fid,1,'SGT_ROOT',2,SGT_ROOT)
call string_conf(fid,1,'SGT_STAX',2,SGT_STAX)
call string_conf(fid,1,'SGT_STAY',2,SGT_STAY)
call string_conf(fid,1,'SGT_STAZ',2,SGT_STAZ)
!call string_conf(fid,1,'SGT_LIST',2,sgt_list)

call string_conf(fid,1,'EVENT_ROOT',2,EVENT_ROOT)
call string_conf(fid,1,'EVENT_OUTD',2,EVENT_OUTD)

! numbers
call string_conf(fid,1,'maximum_of_pick_per_component',2,num_pick)
call string_conf(fid,1,'number_of_freqband',2,num_freq)

close(fid)
end subroutine init_kernel

subroutine init_filter(fnm_conf)
character (len=*),intent(in) :: fnm_conf

character (len=SEIS_STRLEN) :: fnm_filter,str
integer :: fid,gid,n,m

fid=5002; gid=5003
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'<filter_file>',1,str)
do n=1,num_freq
   read(fid,*) fnm_filter

   if (trim(fnm_filter)/='none') then
      open(gid,file=trim(fnm_filter),status='old')
      read(gid,*) NOD1
      NOD=NOD1-1; NFACT=3*NOD
      NEL=knt+2*NFACT
      if (.not. allocated(filt_a)) then
         call alloc_filt_var(num_freq,NOD1,NEL)
      end if
      do m=1,NOD1
         !read(gid,*) filt_a(m,n),filt_b(m,n)
         read(gid,"(2(f20.16))") filt_a(m,n),filt_b(m,n)
      end do
      close(gid)
   end if
end do
close(fid)
end subroutine init_filter

subroutine init_synthetic(fnm_syn,T,V,U)
character (len=*),intent(in) :: fnm_syn
real(SP),dimension(:),intent(out) :: T,V,U
integer :: fid,n

fid=5001
open(fid,file=trim(fnm_syn),status='old')
  do n=1,knt
     read(fid,*) T(n),V(n),U(n)
  end do
close(fid)
end subroutine init_synthetic

subroutine export_info(pnm,T,V,U,twin,tidx,fcta,fctb)
character (len=*),intent(in) :: pnm
real(SP),dimension(:),intent(in) :: T,V,U
real(SP),dimension(:),intent(in) :: twin
integer,dimension(:),intent(in) :: tidx
real(DP),dimension(:),intent(in) :: fcta,fctb

character (len=SEIS_STRLEN) :: filenm
integer :: n,fid
fid=6001

filenm=trim(pnm)//'/'//'kernel.para.log'
open(fid,file=trim(filenm),status='unknown')
  write(fid,*) twin
  write(fid,*) tidx
  write(fid,"(f20.16)") fcta
  write(fid,"(f20.16)") fctb
close(fid)

filenm=trim(pnm)//'/'//'synthetic.filtered.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,knt
   write(fid,*) T(n),V(n),U(n)
end do
close(fid)
end subroutine export_info

subroutine export_debug(pnm,T,V,U)
character (len=*),intent(in) :: pnm
real(SP),dimension(:),intent(in) :: T,V,U

character (len=SEIS_STRLEN) :: filenm
integer :: n,fid
fid=6001

filenm=trim(pnm)//'/'//'kernel.debug.log'
open(fid,file=trim(filenm),status='unknown')
  write(fid,*) 'knt=',knt
  write(fid,*) 'NEL=',NEL
  write(fid,*) 'NFACT=',NFACT
  write(fid,*) 'NOD=',NOD
close(fid)

filenm=trim(pnm)//'/'//'synthetic.orig.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,knt
   write(fid,*) T(n),V(n),U(n)
end do
close(fid)

filenm=trim(pnm)//'/'//'filter.vec.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,NEL
   write(fid,*) vecx(n),vecy(n)
end do
close(fid)
end subroutine export_debug

subroutine stress2strain(Exx,Eyy,Ezz,Exy,Exz,Eyz,lam,miu,scl)
! Same arithmetic as the original; the three medium coefficients are computed
! once per grid point instead of once per grid point and time sample.
real(SP),dimension(:,:,:,:),intent(inout) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP),dimension(:,:,:),intent(in) :: lam,miu
integer,dimension(SEIS_GEO),intent(in) :: scl
real(SP) :: E0
real(SP),dimension(:,:,:),allocatable :: C1,C2,C3,C21
integer i,j,k,n,n0

allocate(C1(scl(1),scl(2),scl(3)),C2(scl(1),scl(2),scl(3)), &
         C3(scl(1),scl(2),scl(3)),C21(scl(1),scl(2),scl(3)))
do k=1,scl(3)
do j=1,scl(2)
do i=1,scl(1)
   C1(i,j,k)=(lam(i,j,k)+miu(i,j,k))/(miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   C2(i,j,k)=-lam(i,j,k)/(2.0*miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   C3(i,j,k)=0.5*(1.0/miu(i,j,k))
   C21(i,j,k)=C2(i,j,k)-C1(i,j,k)
end do
end do
end do

n0=size(Exx,4)
!$omp parallel do default(shared) private(n,k,j,i,E0)
do n=1,n0
do k=1,scl(3)
do j=1,scl(2)
do i=1,scl(1)
   E0=C2(i,j,k)*(Exx(i,j,k,n)+Eyy(i,j,k,n)+Ezz(i,j,k,n))
   Exx(i,j,k,n)=E0-C21(i,j,k)*Exx(i,j,k,n)
   Eyy(i,j,k,n)=E0-C21(i,j,k)*Eyy(i,j,k,n)
   Ezz(i,j,k,n)=E0-C21(i,j,k)*Ezz(i,j,k,n)
   Exy(i,j,k,n)=C3(i,j,k)*Exy(i,j,k,n)
   Exz(i,j,k,n)=C3(i,j,k)*Exz(i,j,k,n)
   Eyz(i,j,k,n)=C3(i,j,k)*Eyz(i,j,k,n)
end do
end do
end do
end do
!$omp end parallel do
deallocate(C1,C2,C3,C21)
end subroutine stress2strain

subroutine sgt2scatter(                   &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,nt1,nt2,dt,cblk)

integer,intent(in) :: nt1,nt2
real(SP),dimension(:,:,:,:),intent(inout) :: &
     ExxS,EyyS
real(SP),dimension(:,:,:,:),intent(in) ::    &
               EzzS,ExyS,ExzS,EyzS,          &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP),dimension(:,:,:),intent(in) :: miu,lam
integer,dimension(SEIS_GEO),intent(in) :: cblk
real(SP),intent(in) :: dt

integer :: n,i,j,k
real(kind=SP) :: x1,x2

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)
   E11R(1:nt2)=ExxR(i,j,k,1:nt2);E22R=EyyR(i,j,k,1:nt2);E33R=EzzR(i,j,k,1:nt2)
   E12R(1:nt2)=ExyR(i,j,k,1:nt2);E13R=ExzR(i,j,k,1:nt2);E23R=EyzR(i,j,k,1:nt2)
   E11S(1:nt2)=ExxS(i,j,k,1:nt2);E22S=EyyS(i,j,k,1:nt2);E33S=EzzS(i,j,k,1:nt2)
   E12S(1:nt2)=ExyS(i,j,k,1:nt2);E13S=ExzS(i,j,k,1:nt2);E23S=EyzS(i,j,k,1:nt2)

do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do

   ExxS(i,j,k,n)=x1*2.0_SP*(lam(i,j,k)+2.0_SP*miu(i,j,k))
   EyyS(i,j,k,n)=2.0_SP*(x2-x1)*2.0_SP*miu(i,j,k)

end do

end do
end do
end do
end subroutine sgt2scatter

function kc_wtime() result(t)
! wall clock in seconds
real(DP) :: t
#ifndef KernelMPI
integer(8) :: cnt,rate
#endif
#ifdef KernelMPI
t=MPI_WTIME()
#else
call system_clock(cnt,rate)
t=real(cnt,DP)/real(rate,DP)
#endif
end function kc_wtime

subroutine channel_weights(ic,w)
! Weights of the receiver components (x,y,z) that make up output channel ic:
! ic=1..3 unrotated component, ic=4..6 rotated component (column ic-3 of rot_mat)
integer,intent(in) :: ic
real(SP),intent(out) :: w(3)
if (ic<=3) then
   w=0.0_SP; w(ic)=1.0_SP
else
   w=rot_mat(:,ic-3)
end if
end subroutine channel_weights

subroutine chk_etmp(A)
real(SP),dimension(:,:,:,:),intent(in) :: A
if (allocated(Etmp)) then
   if (any(shape(Etmp)/=shape(A))) deallocate(Etmp)
end if
if (.not. allocated(Etmp)) allocate(Etmp(size(A,1),size(A,2),size(A,3),size(A,4)))
end subroutine chk_etmp

subroutine load_sgt(path,Exx,Eyy,Ezz,Exy,Exz,Eyz,wgt,acc)
! Read the six stress SGT components of one run (all time chunks) into the
! block arrays.  acc=.false.: E = wgt*data.  acc=.true.: E = E + wgt*data
! (used to add the second horizontal-force run).  Replaces the four copies of
! the same read loop of the original code; for wgt=1, acc=.false. the data
! read are identical to the original.
character (len=*),intent(in) :: path
real(SP),dimension(:,:,:,:),intent(inout) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP),intent(in) :: wgt
logical,intent(in) :: acc
character (len=SEIS_STRLEN) :: fl
integer :: ncl,nn,mm,mtt,vid(6)

if (acc) call chk_etmp(Exx)
mm=0; nn=0
mtt=knt
if (snap_tcnt(kid)<knt) mtt=snap_tcnt(kid)
do
   if (mm>=knt) exit
   nn=nn+1
   fl=get_fnm_snapnode_n(trim(path),'sgt_',kid,nn,n_i,n_j,n_k)
   call nfseis_open(fl,ncl)
   call nfseis_inq_varid(ncl, 'Txx', vid(1))
   call nfseis_inq_varid(ncl, 'Tyy', vid(2))
   call nfseis_inq_varid(ncl, 'Tzz', vid(3))
   call nfseis_inq_varid(ncl, 'Txy', vid(4))
   call nfseis_inq_varid(ncl, 'Txz', vid(5))
   call nfseis_inq_varid(ncl, 'Tyz', vid(6))
   if (mm+mtt>knt) mtt=knt-mm
   call read_var(ncl,vid(1),Exx,mm,mtt,wgt,acc)
   call read_var(ncl,vid(2),Eyy,mm,mtt,wgt,acc)
   call read_var(ncl,vid(3),Ezz,mm,mtt,wgt,acc)
   call read_var(ncl,vid(4),Exy,mm,mtt,wgt,acc)
   call read_var(ncl,vid(5),Exz,mm,mtt,wgt,acc)
   call read_var(ncl,vid(6),Eyz,mm,mtt,wgt,acc)
   mm=mm+mtt
   call nfseis_close(ncl)
end do
end subroutine load_sgt

subroutine read_var(ncl,vid,A,m0,mtt,wgt,acc)
integer,intent(in) :: ncl,vid,m0,mtt
real(SP),dimension(:,:,:,:),intent(inout) :: A
real(SP),intent(in) :: wgt
logical,intent(in) :: acc
integer :: ie
if (acc) then
   ie=nf90_get_var(ncl,vid,Etmp(:,:,:,1:mtt),(/bsubs,1/),(/bsubc,mtt/),(/bsubt,1/))
   if (ie /= nf90_noerr) call nfseis_except(ie,'tomo_kernel: read sgt fail')
   A(:,:,:,m0+1:m0+mtt)=A(:,:,:,m0+1:m0+mtt)+wgt*Etmp(:,:,:,1:mtt)
else
   ie=nf90_get_var(ncl,vid,A(:,:,:,m0+1:m0+mtt),(/bsubs,1/),(/bsubc,mtt/),(/bsubt,1/))
   if (ie /= nf90_noerr) call nfseis_except(ie,'tomo_kernel: read sgt fail')
   if (wgt/=1.0_SP) A(:,:,:,m0+1:m0+mtt)=wgt*A(:,:,:,m0+1:m0+mtt)
end if
end subroutine read_var

subroutine addtr(dst,w,src,i,j,k,ne)
! dst(1:ne) = dst(1:ne) + w*src(i,j,k,1:ne)
real(SP),dimension(:),intent(inout) :: dst
real(SP),intent(in) :: w
real(SP),dimension(:,:,:,:),intent(in) :: src
integer,intent(in) :: i,j,k,ne
integer :: n
do n=1,ne
   dst(n)=dst(n)+w*src(i,j,k,n)
end do
end subroutine addtr

subroutine sgt2scatter_ch(wr,ne,dt)
! Scattering integral for ONE output channel.
!   receiver strain  R = wr(1)*R_x + wr(2)*R_y + wr(3)*R_z   (receiver rotation)
!   event strain     S = ExxS..EyzS (already source-side combined)
!   Ca(n) = 2(lam+2mu) * dt*sum_{m=1}^{n-1} tr(R(m)) tr(S(n-m))
!   Cb(n) = 4 mu       * dt*sum_{m=1}^{n-1} [ R(m):S(n-m) - tr(R(m)) tr(S(n-m)) ]
! (same definition as the original sgt2scatter3, n=1..ne; n>ne is set to 0).
! The sums skip the leading zeros of R and of S: exact, since those terms are 0.
real(SP),intent(in) :: wr(3),dt
integer,intent(in) :: ne
integer :: i,j,k,n,m,mr0,ms0,mhi
real(SP) :: x1,x2,coefa,coefb
real(SP),dimension(:,:),allocatable :: Rw,Sw
real(SP),dimension(:),allocatable :: Rt,St,oa,ob
logical :: ux,uy,uz

ux=(wr(1)/=0.0_SP); uy=(wr(2)/=0.0_SP); uz=(wr(3)/=0.0_SP)

!$omp parallel default(shared) private(i,j,k,n,m,mr0,ms0,mhi,x1,x2,coefa,coefb,Rw,Sw,Rt,St,oa,ob)
allocate(Rw(ne,6),Sw(ne,6),Rt(ne),St(ne),oa(knt),ob(knt))
!$omp do collapse(2) schedule(dynamic)
do k=1,bsubc(3)
do j=1,bsubc(2)
do i=1,bsubc(1)
   Rw=0.0_SP
   if (ux) then
      call addtr(Rw(:,1),wr(1),ExxX,i,j,k,ne); call addtr(Rw(:,2),wr(1),EyyX,i,j,k,ne)
      call addtr(Rw(:,3),wr(1),EzzX,i,j,k,ne); call addtr(Rw(:,4),wr(1),ExyX,i,j,k,ne)
      call addtr(Rw(:,5),wr(1),ExzX,i,j,k,ne); call addtr(Rw(:,6),wr(1),EyzX,i,j,k,ne)
   end if
   if (uy) then
      call addtr(Rw(:,1),wr(2),ExxY,i,j,k,ne); call addtr(Rw(:,2),wr(2),EyyY,i,j,k,ne)
      call addtr(Rw(:,3),wr(2),EzzY,i,j,k,ne); call addtr(Rw(:,4),wr(2),ExyY,i,j,k,ne)
      call addtr(Rw(:,5),wr(2),ExzY,i,j,k,ne); call addtr(Rw(:,6),wr(2),EyzY,i,j,k,ne)
   end if
   if (uz) then
      call addtr(Rw(:,1),wr(3),ExxZ,i,j,k,ne); call addtr(Rw(:,2),wr(3),EyyZ,i,j,k,ne)
      call addtr(Rw(:,3),wr(3),EzzZ,i,j,k,ne); call addtr(Rw(:,4),wr(3),ExyZ,i,j,k,ne)
      call addtr(Rw(:,5),wr(3),ExzZ,i,j,k,ne); call addtr(Rw(:,6),wr(3),EyzZ,i,j,k,ne)
   end if
   Sw(:,1)=ExxS(i,j,k,1:ne); Sw(:,2)=EyyS(i,j,k,1:ne); Sw(:,3)=EzzS(i,j,k,1:ne)
   Sw(:,4)=ExyS(i,j,k,1:ne); Sw(:,5)=ExzS(i,j,k,1:ne); Sw(:,6)=EyzS(i,j,k,1:ne)

   Rt=Rw(:,1)+Rw(:,2)+Rw(:,3)
   St=Sw(:,1)+Sw(:,2)+Sw(:,3)
   Rw(:,4:6)=2.0_SP*Rw(:,4:6)      ! shear terms count twice

   ! first non-zero sample of R and S
   mr0=ne+1
   do m=1,ne
      if (any(Rw(m,:)/=0.0_SP)) then
         mr0=m; exit
      end if
   end do
   ms0=ne+1
   do m=1,ne
      if (any(Sw(m,:)/=0.0_SP)) then
         ms0=m; exit
      end if
   end do

   coefa=2.0_SP*(lambda(i,j,k)+2.0_SP*mu(i,j,k))
   coefb=4.0_SP*mu(i,j,k)
   oa=0.0_SP; ob=0.0_SP
   do n=1,ne
      mhi=n-ms0                    ! S(n-m) is zero for n-m < ms0
      x1=0.0_SP; x2=0.0_SP
      do m=mr0,mhi
         x1=x1+Rt(m)*St(n-m)
         x2=x2+Rw(m,1)*Sw(n-m,1)+Rw(m,2)*Sw(n-m,2)+Rw(m,3)*Sw(n-m,3) &
              +Rw(m,4)*Sw(n-m,4)+Rw(m,5)*Sw(n-m,5)+Rw(m,6)*Sw(n-m,6)
      end do
      oa(n)=dt*x1*coefa
      ob(n)=dt*(x2-x1)*coefb
   end do
   Ca(i,j,k,1:knt)=oa(1:knt)
   Cb(i,j,k,1:knt)=ob(1:knt)
end do
end do
end do
!$omp end do
deallocate(Rw,Sw,Rt,St,oa,ob)
!$omp end parallel
end subroutine sgt2scatter_ch





subroutine scatter2kernel(V,U,Sa,Sb, &
           nt1,nt2,dt,cblk,Kap,Kaq,Kbp,Kbq)

real(SP),dimension(:),intent(in) :: V,U
real(SP),dimension(:,:,:,:),intent(in) :: Sa,Sb
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt1,nt2
real(SP),intent(in) :: dt

real(SP),dimension(:,:,:),intent(out) :: Kap,Kaq,Kbp,Kbq

real(SP) :: V0,U0
integer :: n,i,j,k

! normalization
!V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
!do m=nt1+1,nt2-1
!   V0=V0+V(m)**2.0*kdt; U0=U0+U(m)**2.0*kdt
!end do
!V=V/V0/SGT_M0; U=U/U0/SGT_M0

! normalization
V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
do n=nt1+1,nt2-1
   V0=V0+V(n)**2.0*dt; U0=U0+U(n)**2.0*dt
end do
!V=V/V0/SGT_M0; U=U/U0/SGT_M0
V0=V0*SGT_M0; U0=U0*SGT_M0

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   Ka(nt1:nt2)=Sa(i,j,k,nt1:nt2);
   Kb(nt1:nt2)=Sb(i,j,k,nt1:nt2);

   ! Forming kernels for tau_p and tau_q.
   Kap(i,j,k)=0.5*dt*(V(nt1)*Ka(nt1)+V(nt2)*Ka(nt2))
   Kaq(i,j,k)=0.5*dt*(U(nt1)*Ka(nt1)+U(nt2)*Ka(nt2))
   Kbp(i,j,k)=0.5*dt*(V(nt1)*Kb(nt1)+V(nt2)*Kb(nt2))
   Kbq(i,j,k)=0.5*dt*(U(nt1)*Kb(nt1)+U(nt2)*Kb(nt2))
   do n=nt1+1,nt2-1
      Kap(i,j,k)=Kap(i,j,k)+V(n)*Ka(n)*dt
      Kaq(i,j,k)=Kaq(i,j,k)+U(n)*Ka(n)*dt
      Kbp(i,j,k)=Kbp(i,j,k)+V(n)*Kb(n)*dt
      Kbq(i,j,k)=Kbq(i,j,k)+U(n)*Kb(n)*dt
   end do
   
   Kap(i,j,k)=Kap(i,j,k)/V0
   Kaq(i,j,k)=Kaq(i,j,k)/U0
   Kbp(i,j,k)=Kbp(i,j,k)/V0
   Kbq(i,j,k)=Kbq(i,j,k)/U0

end do
end do
end do
end subroutine scatter2kernel

subroutine scatter_filtfilt(w0,w,bval,aval,nt,cblk)
real(SP),dimension(:,:,:,:),intent(in) :: w0
real(SP),dimension(:,:,:,:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt
integer :: i,j,k,n

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NFACT+1:NFACT+nt)=w0(i,j,k,1:nt)
   do n=1,NFACT
      vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(i,j,k,1:nt)=vecy(NFACT+1:NFACT+nt)
   
end do
end do
end do
end subroutine scatter_filtfilt

subroutine scatter_filter(w0,w,bval,aval,nt,cblk)
real(SP),dimension(:,:,:,:),intent(in) :: w0
real(SP),dimension(:,:,:,:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt
integer :: i,j,k
integer :: n

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NFACT+1:NFACT+nt)=w0(i,j,k,1:nt)
   do n=1,NFACT
      !vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(n)=0.0
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(i,j,k,1:nt)=vecy(NFACT+1:NFACT+nt)
   
end do
end do
end do
end subroutine scatter_filter

subroutine syn_filtfilt(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt
integer :: n

   vecx(NFACT+1:NFACT+nt)=w0(1:nt)
   do n=1,NFACT
      vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   !do n=1,NFACT
   !   vecx(n)=2*w0(1)-w0(NFACT+1-(n-1))
   !end do
   !do n=NFACT+1,NFACT+nt
   !   vecx(n)=w0(n-NFACT)
   !end do
   !do n=NFACT+nt+1,NEL
   !   vecx(n)=2*w0(nt)-w0(nt-1-(n-(NFACT+nt+1)))
   !end do
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(1:nt)=vecy(NFACT+1:NFACT+nt)

end subroutine syn_filtfilt

subroutine syn_filter(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt
integer :: n

   vecx(NFACT+1:NFACT+nt)=w0(1:nt)
   !do n=1,NFACT
   !   vecx(n)=2*vecx(1)-vecx(NFACT+1-(n-1))
   !   vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   !end do
   do n=1,NFACT
      !vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(n)=0.0
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(1:nt)=vecy(NFACT+1:NFACT+nt)

end subroutine syn_filter

subroutine filtfilt(NOD,NEL,b,a,x,y,w)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real(DP),dimension(:) :: x,y,w

call filter(NOD,NEL,b,a,x,y)
w=y(NEL:1:-1)
call filter(NOD,NEL,b,a,w,y)
y=y(NEL:1:-1)
end subroutine filtfilt

subroutine filter(NOD,NEL,b,a,x,y)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real(DP),dimension(:) :: x,y
integer :: n,m

do n=1,NEL
   y(n)=b(1)*x(n)
   do m=2,min(NOD+1,n)
      y(n)=y(n)+b(m)*x(n-m+1)-a(m)*y(n-m+1)
   end do
   y(n)=y(n)/a(1)
end do
end subroutine filter

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef KernelMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef KernelMPI
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

subroutine warning_print(msg)
  character (len=*),intent(in) :: msg
  print *, idstr//" warning :"//trim(msg)
end subroutine warning_print

end program SI_ker_sta

