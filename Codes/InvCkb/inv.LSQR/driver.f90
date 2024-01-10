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
        integer :: glistio, gkernio, nlines, ierr

        ! loops
        integer :: i, j

        ! shen used line length 300 and token length 132. not sure why. 
        ! not going to change it.
        character (len = 300) :: line
        character (len = 132) :: token

        ! user input parameters. taken from get_input(...)
        character (len = 10) :: mode
        character (len = 132) :: gmat_list, smooth_op, gkernel_file
        integer(I64) :: inv_struct_cmps, nstations, nevents, nlocations
        integer(I64) :: max_data, max_station, max_event, max_location, &
            max_mval, max_xval, max_size, max_row, max_bd, max_elements, &
            max_row_smooth, max_bd_smooth, max_elements_smooth
        real(DP) :: station_weighting, event_weighting, location_weighting, &
            damping, smoothing, kernel_threshold, data_threshold

        ! matrix data
        integer(I64) :: nelements
        integer :: ncoefficients, nblocks
        real(DP) :: rhs, weight_data, tmp
        integer, allocatable :: stations(:), events(:), locations(:)
        real(DP), allocatable :: location_weightings(:)
        real(SP), dimension(3) :: kernels ! kernel index in row
        integer, dimension(:), allocatable :: value_indices, indices
        real(SP), dimension(:), allocatable :: coefficients_sp
        real(DP), dimension(:), allocatable :: coefficients_dp, mvals, errxs

        ! extreme value metrics
        real(DP) :: rhs_min = 1.e7, rhs_max = -1.e7, &
            rhs_min_weight = 1.e7, rhs_max_weight = -1.e7

        ! get user input
        call get_input(mode, gmat_list, smooth_op, inv_struct_cmps, nstations, &
            nevents, nlocations, station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold)

        ! init smoothing
        call get_smoothing(smooth_op, inv_struct_cmps, max_row_smooth, &
            max_bd_smooth, max_elements_smooth)


        call checkpoint("Retrieved smoothing data")

        ! init damping
        call get_damping(gmat_list, inv_struct_cmps, nstations, nevents, &
            nlocations, max_data, max_station, max_event, max_location, &
            max_mval, max_xval, max_size, max_row, max_bd, max_elements)

        call checkpoint("Retrieved damping data")

        max_row = max_row + max_row_smooth
        max_elements = max_elements + max_elements_smooth
        max_bd = max(max_bd, max_bd_smooth)
        max_size = max(max_size, max_row)

        call alloc_vars(nstations, nevents, nlocations, max_xval, max_bd, &
            stations, events, locations, location_weightings, value_indices, &
            indices, coefficients_sp, coefficients_dp, mvals, errxs)

        ! read gmat_list file
        glistio = 10000; gkernio = 10001
        open(glistio, file = trim(gmat_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("failed to open gmat_list")

        ! loop through file line by line
        nlines = 0
        read_loop: do
            read(glistio, "(A300)", iostat = ierr) line
            if (ierr < 0 .or. nlines == 100) exit

            ! skip comment lines
            call advance_to(" ", line, token)
            if (token(1:1) == '#') cycle
            call advance_to(" ", line, token)

            nlines = nlines + 1

            nelements = 0
            do i = 1, inv_struct_cmps
                ! token currently holds kernel filename
                open(gkernio, file = trim(token), status = "old", &
                    iostat = ierr, form = "unformatted")
                if (ierr /= 0) &
                    call throw("open kernel err: "//trim(token))

                ! read data from kernel file
                read(gkernio) nblocks, kernels(i)
                read(gkernio) ncoefficients
                !read(gkernio) (value_indices(nelements+j), &
                !    coefficients_sp(nelements+j), j=1, ncoefficients)
                close(gkernio)

                ! add coefficients
                do j = 1, ncoefficients
                    value_indices(nelements+j) = &
                        value_indices(nelements+j)+nblocks*(i-1)
                end do

                nelements = nelements + 1

                ! advance to next token
                call advance_to(" ", line, token)
            end do

            ! measurement, weight, stationid
            read(line, *) rhs, weight_data
            call advance_to(" ", line, token)
            call advance_to(" ", line, token)

            ! check for extreme measurements
            if (abs(rhs) > data_threshold) cycle

            ! station elements
            do i = 1, nstations
                read(line, *) stations(i)
                call advance_to(" ", line, token)
            end do

            ! event elements
            do i = 1, nevents
                read(line, *) events(i)
                call advance_to(" ", line, token)
            end do

            ! location elements and weighting
            do i = 1, nlocations
                read(line, *) locations(i), location_weightings(i)
                call advance_to(" ", line, token)
                call advance_to(" ", line, token)
            end do

            ! check kernel value
            ncoefficients = 0
            do i = 1, nelements
                j = value_indices(i)
                tmp = real(coefficients_sp(i), DP)
                if (dabs(tmp) > kernel_threshold) then
                    ncoefficients = ncoefficients + 1
                    indices(ncoefficients) = j
                    coefficients_dp(ncoefficients) = tmp
                end if
            end do

            ! station term
            do i = 1, nstations
                ncoefficients = ncoefficients + 1
                indices(ncoefficients) = max_mval + abs(stations(i))
                coefficients_dp(ncoefficients) = sign(station_weighting, &
                    real(stations(i), DP))
            end do

            ! event term
            do i = 1, nevents
                ncoefficients = ncoefficients + 1
                indices(ncoefficients) = max_mval + nstations + abs(events(i))
                coefficients_dp(ncoefficients) = sign(event_weighting, &
                    real(events(i), DP))
            end do

            ! location term
            do i = 1, nlocations
                ncoefficients = ncoefficients + 1
                indices(ncoefficients) = max_mval + nstations + nevents + &
                    locations(i)
                coefficients_dp(ncoefficients) = location_weightings(i) * &
                    location_weighting
            end do

            ! calculate rhs min and max
            rhs_min = min(rhs_min, rhs)
            rhs_max = max(rhs_max, rhs)

            ! covariance matrix
            ! TODO: move into previous loops
            do i = 1, ncoefficients
                coefficients_dp(i) = coefficients_dp(i) * weight_data
            end do

            ! rhs min and max weighting
            rhs_min_weight = min(rhs_min_weight, rhs)
            rhs_max_weight = max(rhs_max_weight, rhs)

            ! add row to matrix
            !lldrow(coefficients_dp, indices, ncoefficients, real(rhs, DP))

            if (modulo(nlines, 1000) == 0) &
                print "('Read through line ', I0)", nlines
        end do read_loop

        close(glistio)

        call checkpoint("Built kernel matrix")
    end subroutine main

    subroutine alloc_vars(nstations, nevents, nlocations, nxvals, nbds, &
        stations, events, locations, location_weightings, value_indices, &
        indices, coefficients_sp, coefficients_dp, mvals, errxs)
        integer(I64), intent(in) :: nstations, nevents, nlocations, nxvals, nbds
        integer, dimension(:), allocatable, intent(inout) :: stations, events, &
            locations, value_indices, indices
        real(DP), dimension(:), allocatable, intent(inout) :: location_weightings, &
            coefficients_dp, mvals, errxs
        real(SP), dimension(:), allocatable, intent(inout) :: coefficients_sp

        if (nstations > 0) then
            allocate(stations(nstations))
            stations = 0.0
        end if

        if (nevents > 0) then
            allocate(events(nevents))
            events = 0.0
        end if

        if (nlocations > 0) then
            allocate(locations(nlocations))
            locations = 0.0
            allocate(location_weightings(nlocations))
            location_weightings = 0.0_DP
        end if

        allocate(value_indices(nbds))
        value_indices = 0
        allocate(indices(nbds))
        indices = 0

        allocate(coefficients_sp(nbds))
        coefficients_sp = 0.0
        allocate(coefficients_dp(nbds))
        coefficients_dp = 0.0_DP
        allocate(mvals(nxvals))
        mvals = 0.0_DP
        allocate(errxs(nxvals))
        errxs = 0.0_DP
    end subroutine alloc_vars

    subroutine advance_to(delimiter, input, item)
        character (len = *), intent(in) :: delimiter
        character (len = *), intent(inout) :: input
        character (len = *), intent(out) :: item
        integer :: n ! next delimiter index

        ! get index of next delimiter
        n = index(input(1:len_trim(input)), delimiter)

        ! store current token
        item = " "
        item(1:n) = input(1:n)

        ! advance input
        input(1:) = input(n+1:)

        ! trim leading spaces
        input = adjustl(input)
    end subroutine

    subroutine get_input(mode, gmat_list, smooth_op, inv_struct_cmps, &
        nstations, nevents, nlocations, station_weighting, event_weighting, &
        location_weighting, damping, smoothing, kernel_threshold, &
        data_threshold)

        ! mode
        integer(I64) :: choice
        character (len = 10) :: mode, mode0(4)

        ! parameters
        character (len = 132), intent(out) :: gmat_list, smooth_op
        integer(I64), intent(out) :: inv_struct_cmps, nstations, nevents, &
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

        write(*, "(A)", advance="no") "file name of list of G matrix, data, and weighting: "
        read(*, "(A)") gmat_list
        write(*, "(/, 'Using file: ', A)") trim(gmat_list)

        write(*, "(A)", advance="no") "file name of smoothing operator: "
        read(*, "(A)") smooth_op
        write(*, "(/, 'Using file: ', A)") trim(smooth_op)

        write(*, "(A)", advance="no") "# of inverted structure components: "
        read(*, *) inv_struct_cmps
        write(*, "(/, 'Using: ', I0)") inv_struct_cmps

        write(*, "(A)", advance="no") "# of stations terms per measurement and weighting: "
        read(*, *) nstations, station_weighting
        write(*, "(/, 'Using: ', I0, ' / ', F0.4)") nstations, &
            station_weighting

        write(*, "(A)", advance="no") "# of event terms per measurement and weighting: "
        read(*, *) nevents, event_weighting
        write(*, "(/, 'Using: ', I0, ' / ', F0.4)") nevents, event_weighting

        write(*, "(A)", advance="no") "# of location terms per measurement and weighting: "
        read(*, *) nlocations, location_weighting
        write(*, "(/, 'Using: ', I0, ' / ', F0.4)") nlocations, &
            location_weighting

        write(*, "(A)", advance="no") "damping and smoothing: "
        read(*, *) damping, smoothing
        write(*, "(/, 'Using: ', F0.4, ' / ', F0.4)") damping, smoothing

        write(*, "(A)", advance="no") "kernel and measurement thresholds: "
        read(*, *) kernel_threshold, data_threshold
        write(*, "(/, 'Using: ', F0.4, ' / ', F0.4)") kernel_threshold, &
            data_threshold
    end subroutine get_input

    subroutine get_smoothing(gmat_smooth, inv_struct_cmps, max_row, max_bd, &
        max_elements)
        character (len = *), intent(in) :: gmat_smooth
        integer(I64), intent(in) :: inv_struct_cmps
        integer(I64), intent(out) :: max_row, max_bd, max_elements
        integer :: fid, nblocks, ierr

        fid = 10000
        max_row = 0
        max_bd = 0
        max_elements = 0
        
        open(fid, file = trim(gmat_smooth), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("init_smoothing open err: "//trim(gmat_smooth))
        read(fid, *) nblocks, max_row, max_bd, max_elements
        close(fid)

        max_bd = max_bd * inv_struct_cmps
        max_row = max_row * inv_struct_cmps
        max_elements = max_elements * inv_struct_cmps
    end subroutine get_smoothing

    subroutine get_damping(gmat_list, inv_struct_cmps, nstations, nevents, &
        nlocations, max_data, max_station, max_event, max_location, max_mval, &
        max_xval, max_size, max_row, max_bd, max_elements)
        character (len = *), intent(in) :: gmat_list
        integer(I64), intent(in) :: inv_struct_cmps, nstations, nevents, nlocations
        integer(I64), intent(out) :: max_data, max_station, max_event, &
            max_location, max_mval, max_xval, max_size, max_row, max_bd, &
            max_elements

        character(len = 132) :: kernel_file, token
        character(len = 300) :: line
        integer :: stations(nstations), events(nevents), locations(nlocations)
        real :: weights(nlocations)
        integer :: glistio, gkernio, nlines, ierr
        integer :: i, j
        integer :: nblocks, ncoefficients
        integer(I64) :: nelements
        real(SP) :: kernel

        real(DP) :: rhs, weight_data

        max_data = 0
        max_station = 0
        max_event = 0
        max_location = 0
        max_mval = 0
        max_xval = 0
        max_size = 0
        max_row = 0
        max_bd = 0
        max_elements = 0

        glistio = 10000; gkernio = 10001
        open(glistio, file = trim(gmat_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("failed to open gmat_list")

        nlines = 0
        do
            read(glistio, "(A300)", iostat = ierr) line
            if (ierr < 0 .or. nlines == 100) exit

            ! skip comment lines
            call advance_to(" ", line, token)
            if (token(1:1) == '#') cycle
            call advance_to(" ", line, token)

            nlines = nlines + 1
            max_data = max_data + 1

            nelements = 0
            do i = 1, inv_struct_cmps
                ! token currently holds kernel filename
                open(gkernio, file = trim(token), status = "old", &
                    iostat = ierr, form = "unformatted")
                if (ierr /= 0) &
                    call throw("open kernel err: "//trim(token))

                ! read data from kernel file
                read(gkernio) nblocks, kernel
                read(gkernio) ncoefficients

                close(gkernio)

                nelements = nelements + ncoefficients

                ! advance to next token
                call advance_to(" ", line, token)
            end do

            read(line, *) rhs, weight_data
            call advance_to(" ", line, token)
            call advance_to(" ", line, token)

            do i = 1, nstations
                read(line, *) stations(i)
                call advance_to(" ", line, token)
                max_station = max(max_station, abs(stations(i)))
            end do

            do i = 1, nevents
                read(line, *) events(i)
                call advance_to(" ", line, token)
                max_event = max(max_event, abs(events(i)))
            end do

            do i = 1, nlocations
                read(line, *) locations(i)
                call advance_to(" ", line, token)
                call advance_to(" ", line, token)
                max_location = max(max_location, abs(locations(i)))
            end do

            max_bd = max(max_bd, nelements)
            max_elements = max_elements + nelements

            if (modulo(nlines, 1000) == 0) &
                print "('Read through line ', I0)", nlines
        end do

        close(glistio)

        ! compute final values
        max_mval = nblocks * inv_struct_cmps
        max_xval = max_mval + max_station + max_event + max_location
        max_elements = max_elements + (nstations + nevents + nlocations) * &
            max_data + max_xval
        max_row = max_data + max_xval
        max_bd = max_bd + nstations + nevents + nlocations
        max_size = max(max_xval, max_row)
    end subroutine get_damping

    subroutine worker(nworkers, rank)
        integer, intent(in) :: nworkers, rank
        integer, dimension(mpi_status_size) :: status
        integer :: nlines, ierr

        ! block for signal from main
        call mpi_probe(0, 0, mpi_comm_world, status, ierr)

        print "('Hello from worker ', I0)", rank
    end subroutine worker

    subroutine checkpoint(msg)
        character (len = *), intent(in) :: msg

        write(*, "(A, /, A)") trim(msg), hline
    end subroutine checkpoint

    subroutine throw(emsg)
        character (len = *), intent(in) :: emsg

        write(*, "('Fatal: ', A)") trim(emsg)
        stop 1
    end subroutine throw
end program
