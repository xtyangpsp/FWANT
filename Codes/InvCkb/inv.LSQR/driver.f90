program driver
    use mpi
    use, intrinsic :: iso_fortran_env, only : &
        I8 => int8, &
        I16 => int16, &
        I32 => int32, &
        I64 => int64, &
        SP => real32, &
        DP => real64

    implicit none

    ! openmp
    integer :: rc, rank, nprocs, nworkers, mainrank
    real(DP) :: starttime, endtime

    ! init openmp, get n of active processes
    call mpi_init(rc)
    call mpi_comm_size(mpi_comm_world, nprocs, rc)

    mainrank = 0
    nworkers = nprocs - 1

    ! current process
    call mpi_comm_rank(mpi_comm_world, rank, rc)

    if (rank .eq. mainrank) then
        ! in main process
        starttime = mpi_wtime()
        call main(nworkers)
        endtime = mpi_wtime()

        print '("--- Timing: ", f0.4, " sec on ", i0, " workers")', &
            starttime - endtime, nworkers
    else
        ! in a worker thread
        call worker(nworkers, rank)
    endif

    call mpi_finalize(rc)

contains

    subroutine main(nworkers)
        integer, intent(in) :: nworkers
        integer :: i, nlines, listid, ierr

        ! init yang module
        call yang_main_init()

        listid = 10000
        open(listid, file = trim(fnm_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call error_except("fnm_list open err in main: "//trim(fnm_list))

        ! get number of lines
        nlines = 0
        !$OMP PARALLEL DO REDUCTION(+:nlines)
        do
            read(listid, *, iostat = ierr)
            if (ierr /= 0) exit
            nlines = nlines + 1
        end do
        !$OMP END PARALLEL DO

        close(listid)

        ! unblock workers
        !$OMP PARALLEL DO
        do i = 1, nworkers

        end do
        !$OMP END PARALLEL DO
    end subroutine main

    subroutine worker(nworkers, rank)
        integer, intent(in) :: nworkers, rank
        integer, dimension(mpi_status_size) :: status
        integer :: nlines, ierr

        ! block for signal from main
        call mpi_probe(0, 0, mpi_comm_world, status, ierr)

        print "('Hello from worker ', I0)", rank
    end subroutine worker

    subroutine throw(emsg)
        character (len = *), intent(in) :: emsg

        write(*, "('Fatal: ', A") trim(emsg)
        stop 1
    end subroutine throw
end program
