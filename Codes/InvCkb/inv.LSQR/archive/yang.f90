module yang
    use LSQR_mod
    use, intrinsic :: iso_fortran_env, only : I1B => int8,&
        I2B => int16, I4B => int32,  I8B => int64

    implicit none
    private

    character (len = *),parameter ::  &
        hline="----------------------------------------------------"

    public :: yang_init, error_except

contains

    subroutine yang_main_init(fnm_list)
        ! G matrix, smoothing filename
        character (len = 132), intent(out) :: fnm_list
        character (len = 132) :: fnm_smooth

        ! mode selection
        character (len = 10) :: mode0(4), mode
        integer (I8B) :: imultis

        ! component counts and weights
        integer (I8B) :: num_cmp, nsta, nevt, nloc
        real (DP) :: weig_sta, weig_evt, weig_loc

        ! damping, smoothing, and thresholds
        real (DP) :: lamb, eta, ker_thres, data_thres

        mode0=(/'Simple    ','Multiscale','Quelling  ','HybridVQLM'/)

        ! ------------------------ input parameters ------------------------
        write(*,*) hline
        write(*,*)' mode:'

        write(*,*)'     1 = Simple damping;'
        !write(*,*)'     2 x= MultiScaled transformation:'
        !write(*,*)'     3 x= Convolution quelling.'
        !write(*,*)'     4 x= vertically quelling and laterally multiscale.'
        read(*,*) imultis
        mode=mode0(imultis)

        if (imultis/=1) call error_except("only Simple damping valid")

        print *,'  file list of G matrix, data and weighting='
        read(*,'(a)') fnm_list
        write(*,*) " ",trim(fnm_list)

        print *,'  file name of smoothing operator='
        read(*,'(a)') fnm_smooth
        write(*,*) " ",trim(fnm_smooth)

        write(6,*)'  # of inverted structure components'
        read(5,*) num_cmp
        write(*,*) num_cmp

        write(*,*)'  # of station terms per measurement and weighting'
        read(*,*) nsta,weig_sta
        write(*,*) nsta,weig_sta

        write(*,*)'  # of event terms per measurement and weighting'
        read(*,*) nevt,weig_evt
        write(*,*) nevt,weig_evt

        write(*,*)'  # of location terms per measurement and weighting'
        read(*,*) nloc,weig_loc
        write(*,*) nloc,weig_loc

        write(*,*)' dampping and smoothing'
        read(*,*) lamb,eta
        write(*,*) lamb,eta

        write(6,*)'  kernel and measurement threshold'
        read(5,*) ker_thres, data_thres
        write(*,*) " ",ker_thres, data_thres

    end subroutine yang_main_init

    subroutine yang_worker_init(fnm_smooth, num_cmp, nrow_smoth, nbd_smoth, &
        nel_smoth, nsta, nevt, nloc)

        integer (I8B) :: nrow_smoth, nbd_smoth, nel_smoth
        integer(I8B) ::                             &
            num_dat,num_mval,num_xval,num_siz1d,    &
            num_sta,num_evt,num_loc,num_row,num_bd, &
            num_nel

        call init_smoothing(fnm_smooth,num_cmp,nrow_smoth,nbd_smoth,nel_smoth)

        call init_data_d(fnm_list,num_cmp,nsta,nevt,nloc,    &
            num_dat,num_sta,num_evt,num_loc,num_mval,num_xval, &
            num_siz1d,num_row,num_bd,num_nel)
        write(*,*) 'num_nel: ',num_nel
        num_row=num_row+nrow_smoth
        num_nel=num_nel+nel_smoth
        num_bd=max(num_bd,nbd_smoth)
        num_siz1d=max(num_siz1d,num_row)
        write(*,*) 'num_siz1d: ',num_siz1d

        write(*,*) '  # of cmp, sta, evt, loc per data'
        write(*,*) num_cmp,nsta,nevt,nloc
        write(*,*) '  # of dat, mval, xval, size'
        write(*,*) num_dat,num_mval,num_xval,num_siz1d
        write(*,*) '  # of row, band, nel'
        write(*,*) num_row,num_bd,num_nel
        !!debug
        write(*,*) 'allocating variables ...'
        call alloc_var(nsta,nevt,nloc,num_xval,num_bd)
        write(*,*) 'lsqr_init ...'
        call lsqr_init(num_row,num_xval,num_siz1d,num_nel)
        write(*,*) 'lsqr_init: finished ...'
    end subroutine yang_worker_init

    subroutine yang_damping(num_xval)
        do i=1,num_xval
        !do i=1,num_mval
        !write(*,*) i
            coef_dp(1)=lamb
            indx(1)=i
            call plldrow(coef_dp,indx,1,0.0_DP)
        end do

        call lsqr_check("  ---- data, station and event + dampping ----")
    end subroutine yang_damping

subroutine init_smoothing(filenm,ncmp,max_row,max_bd,max_nel)
character (len=*),intent(in) :: filenm
integer(I8B),intent(in) :: ncmp
integer(I8B),intent(out) :: max_row,max_bd,max_nel
integer :: fid,nblk,ierr

fid=5001;
max_row=0; max_bd=0; max_nel=0
open(fid,file=trim(filenm),status='old',iostat=ierr)
if (ierr>0) call error_except("init_smoothing open err:"//trim(filenm))
read(fid,*) nblk,max_row,max_bd,max_nel
close(fid)
max_bd =max_bd *ncmp
max_row=max_row*ncmp
max_nel=max_nel*ncmp
end subroutine init_smoothing

subroutine init_data_damp(filenm,ncmp,nst,nev,nlo,  &
    lamb, eta, &
    max_dat,max_st,max_ev,max_lo,max_mval,max_xval, &
    max_siz,max_row,max_bd,max_nel)
character (len=*),intent(in) :: filenm
integer(I8B),intent(in) :: ncmp,nst,nev,nlo
integer(I8B),intent(out) ::                     &
  max_dat,max_st,max_ev,max_lo,max_mval,max_xval, &
  max_siz,max_row,max_bd,max_nel
character (len=132) :: fnm_k
integer :: idst(nst),idev(nev),idlo(nlo)
real    :: weig(nlo)
!integer :: fid,gid,k,n,nblk,nael,ierr
integer :: fid,gid,k,n,nblk,ierr
integer(I8B)::nael
real(SP) :: ker0
integer(I8B):: icount=0

character (len = 300) :: rcdstr

integer :: i, j, k

fid=5001; gid=5002
max_st=0
max_ev=0
max_lo=0
max_dat=0
max_mval=0 ! model parameterization
max_xval=0 ! m plus station and event terms
max_siz=0  ! max size of vector
max_row=0
max_bd=0
max_nel=0

! number related with kernels
! debug
write(*,*) 'file:',trim(filenm)
open(fid,file=trim(filenm),status='old',iostat=ierr)
if (ierr>0) call error_except("init_G_scale open err:"//trim(filenm))
do
  !!debug
  icount=icount+1

  read(fid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  rcdstr=adjustl(rcdstr); if (rcdstr(1:1)=='#') cycle

  max_dat=max_dat+1
  nael=0
  do k=1,ncmp
     n=index(rcdstr(1:len_trim(rcdstr))," ")
     fnm_k=rcdstr(1:n)
     !write(*,*) trim(fnm_k) !debug
     rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
     open(gid,file=trim(fnm_k),status='old',iostat=ierr,form='unformatted')
!      write(*,*) 'test1'
     if (ierr>0) call error_except("init_G_scale kfile err:"//trim(filenm))
     read(gid) nblk, ker0
     ! write(*,*) 'test2:',nblk, ker0
     read(gid) ncoef
     !!debug
     ! write(*,*) 'test3',ncoef
     close(gid)
     nael=nael+ncoef

  end do
  !!debug
  !write(*,*) 'nael,nlist: ',nael,icount
  read(rcdstr,*) rhs, weig_dat
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  do k=1,nst
      read(rcdstr,*) idst(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_st=max(max_st,abs(idst(k)))
  end do

  do k=1,nev
      read(rcdstr,*) idev(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_ev=max(max_ev,abs(idev(k)))
  end do

  do k=1,nlo
      read(rcdstr,*) idlo(k), weig(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_lo=max(max_lo,abs(idlo(k)))
  end do

  max_bd=max(max_bd,nael)
  max_nel=max_nel+nael
  !debug
  !write(*,*) max_nel,icount
end do
close(fid)

max_mval=nblk*ncmp
max_xval=max_mval+max_st+max_ev+max_lo
max_nel=max_nel+(nst+nev+nlo)*max_dat +max_xval
max_row=max_dat+max_xval
max_bd=max_bd+nst+nev+nlo
max_siz=max(max_xval,max_row)
write(*,*) 'finished init_data_damp'
end subroutine init_data_damp

end module
