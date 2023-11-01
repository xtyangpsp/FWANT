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

        print '("--- Timing: ", f6.4, " sec on ", i0, " workers")', &
            starttime - endtime, nworkers
    else
        ! in a worker thread
        call worker(nworkers, rank)
    endif

contains

    subroutine main(nworkers)
        integer, intent(in) :: nworkers
        integer :: i, data, ierr

        print *, "Running main first..."

        data = 0

        ! unblock workers
        !$OMP PARALLEL DO
        do i = 1, nworkers
            call mpi_send(data, 1, mpi_int, i, 0, mpi_comm_world, ierr)
        enddo
        !$OMP END PARALLEL DO
    end subroutine main

    subroutine worker(nworkers, rank)
        integer, intent(in) :: nworkers, rank
        integer, dimension(mpi_status_size) :: status
        integer :: ierr

        call mpi_probe(0, 0, mpi_comm_world, status, ierr)
        print "('Writing from worker ', I0)", rank
    end subroutine worker
end program
