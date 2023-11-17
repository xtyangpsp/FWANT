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
    !--------------------------------------------------------------------------
    ! solver_WS.f
    !
    ! 08/13 Y.Shen: change mdata and station location files
    !               change the definition of variance (fit in the program).
    ! 08/12/2008 W. Zhang: port to fortran90 format and only keep simple mode

    ! style
    character (len = *), parameter ::  &
        hline="----------------------------------------------------"

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
        ! openmpi
        integer, intent(in) :: nworkers

        ! files + handling
        integer :: glistio, nlines, ierr

        ! shen used line length 300 and token length 132. not sure why. 
        ! not going to change it.
        character (len = 300) :: line
        character (len = 132) :: token

        ! user input parameters
        character (len = 10) :: mode
        character (len = 132) :: gmat_list, smooth_op, gkernel_file
        integer(I8B) :: inv_struct_cmps, nstations, nevents, nlocations
        real(DP) :: station_weighting, event_weighting, location_weighting, &
            damping, smoothing, kernel_threshold, data_threshold

        ! matrix data
        integer(I8B) :: nelements
        integer :: ncoefficients
        real(DP) :: rhs, weight_data

        ! extreme value metrics
        rhs_min = 1.e7; rhs_max = -1.e7
        rhs_min_weight = 1.e7; rhs_max_weight = -1.e7

        ! get user input
        call get_input(mode, gmat_list, smooth_op, inv_struct_cmps, nstations &
            nevents, nlocations, station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold)

        ! read gmat_list file
        glistio = 10000; gkernio = 10001
        open(glistio, file = trim(gmat_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("failed to open gmat_list")

        ! loop through file line by line
        nlines = 0
        read_loop: do
            read(glistio, "(A)") line
            if (ierr /= 0) exit

            ! skip comment lines
            call advance_to(" ", line, token)
            if (token(1:1) == '#') cycle

            nlines = nlines + 1

            do i = 1, inv_struct_cmps
                ! token currently holds kernel filename
                open(gkernio, file = trim(token), status = "old", &
                    iostat = ierr, form = "unformatted")
                if (ierr /= 0) &
                    call throw("open kernel err: "//trim(token))

                ! read data from kernel file

                close(gkernio)
                ! advance to next token
                call advance_to(" ", line, token)
            end do

            ! measurement, weight, stationid
            read(rcdstr, *) rhs, weight_data
            call advance_to(" ", line, token)
            call advance_to(" ", line, token)

            ! check for extreme measurements
            if (abs(rhs) > data_threshold) cycle

            if (modulo(nlines, 1000) == 0) &
                write(*, "(A, /)") "Read 1000 lines"
        end do read_loop
    end subroutine main

    subroutine advance_to(delimiter, input, item)
        character (len = *), intent(in) :: delimiter
        character (len = *), intent(inout) :: input
        character (len = *), intent(out) :: item
        integer :: n ! next delimiter index

        ! get index of next delimiter
        n = index(input(1:len_trim(input), delimiter)

        ! store current token
        item = input(1:n)

        ! advance input
        input(1:) = input(n+1:)
        input = adjustl(input)
    end subroutine

    subroutine get_input(mode, gmat_list, smooth_op, inv_struct_cmps, &
        nstations, nevents, nlocations, station_weighting, event_weighting, &
        location_weighting, damping, smoothing, kernel_threshold, &
        data_threshold)

        ! mode
        integer(I8B) :: choice
        character (len = 10) :: mode, mode0(4)

        ! parameters
        character (len = 132), intent(out) :: gmat_list, smooth_op
        integer(I8B), intent(out) :: inv_struct_cmps, nstations, nevents, &
            nlocations
        real(DP), intent(out) :: station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold

        mode0 = (/'Simple    ','Multiscale','Quelling  ','HybridVQLM'/)

        write(*, *) hline

        write(*, "(A,/)") "Mode: "//NEW_LINE('a')// &
            "1) Simple damping"
        read(*, *) choice
        mode = mode0(choice)

        if (choice /= 1) &
            call throw("only simple damping valid")

        write(*, *) "file name of list of G matrix, data, and weighting: "
        read(*, "(A)") gmat_list
        write(*, "(/, 'Using file: ', A, /)") trim(gmat_list)

        write(*, *) "file name of smoothing operator: "
        read(*, "(A)") smooth_op
        write(*, "(/, 'Using file: ', A, /)") trim(smooth_op)

        write(*, *) "# of inverted structure components: "
        read(*, *) inv_struct_cmps

        write(*, *) "# of stations terms per measurement and weighting"
        read(*, *) nstations, station_weighting

        write(*, *) "# of event terms per measurement and weighting"
        read(*, *) nevents, event_weighting

        write(*, *) "# of location terms per measurement and weighting"
        read(*, *) nlocations, location_weighting

        write(*, *) "damping and smoothing"
        read(*, *) damping, smoothing

        write(*, *) "kernel and measurement thresholds"
        read(*, *) kernel_threshold, data_threshold
    end subroutine get_input

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
