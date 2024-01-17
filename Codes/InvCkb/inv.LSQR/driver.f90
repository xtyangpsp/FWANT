program driver
    use mpi
    use omp_lib
    use CompressedSparseMatrix

    use, intrinsic :: iso_fortran_env, only : &
        I8 => int8, &
        I16 => int16, &
        I32 => int32, &
        I64 => int64, &
        SP => real32, &
        DP => real64

    implicit none

    integer :: rc, rank, nproc, ierr
    real(DP) :: starttime, endtime

    character (len = *), parameter ::  &
        hline="----------------------------------------------------"

    ! init openmp
    call mpi_init(rc)
    call mpi_comm_size(mpi_comm_world, nproc, rc)

    ! current process
    call mpi_comm_rank(mpi_comm_world, rank, rc)

    call proc_info(rank)

    ! sync processes
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (ierr /= 0) call throw("mpi failed to bar")

    call main(rank, nproc)

    call mpi_finalize(rc)
contains

    subroutine main(rank, nproc)
        integer, intent(in) :: rank, nproc
        integer :: nthreads, tid
        integer :: group_size_omp = 1

        ! PARAMETERS
        character(len = 10) :: mode
        character(len = 132) :: gmat_list, smooth_op
        integer(I64)  :: inv_struct_cmps, nstations, nevents, &
            nlocations
        real(DP) :: station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold

        ! FROM SMOOTHING
        integer(I64) :: sm_max_row, sm_max_band, sm_max_elements

        ! FROM DAMPING
        integer(I64) :: max_data, max_station, max_event, &
            max_location, max_m, max_x, max_size, max_row, max_band, &
            max_elements

        !$omp parallel private(nthreads, tid)
            tid = omp_get_thread_num()
            write(*, "('OMP tid=', I0, ' on MPI rank ', I0)") tid, rank

            ! master thread
            if (rank == 0 .AND. tid == 0) then
                nthreads = omp_get_num_threads()
                write(*, "('OpenMP info: nthreads=', I0, ' on MPI rank ', I0)") nthreads, rank
                group_size_omp = nthreads * 10
            end if
        !$omp end parallel

        ! wait for threads
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        if (ierr /= 0) call throw("mpi failed to bar")

        ! load parameters
        call get_input(mode, gmat_list, smooth_op, inv_struct_cmps, nstations, &
            nevents, nlocations, station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold)
        call checkpoint("retrieved user input")

        ! get smoothing
        call get_smoothing(smooth_op, inv_struct_cmps, sm_max_row, sm_max_band, &
            sm_max_elements)
        call checkpoint("retrieved smoothing data")

        ! get damping
        call get_damping(gmat_list, inv_struct_cmps, nstations, nevents, &
            nlocations, max_data, max_station, max_event, max_location, &
            max_m, max_x, max_size, max_row, max_band, max_elements)
        max_row = max_row + sm_max_row
        max_elements = max_elements + sm_max_elements
        max_band = max(max_band, sm_max_band)
        max_size = max(max_size, max_row)
        write(*, "('max # of elements: ', I0)") max_elements
        write(*, "('max size: ', I0)") max_size
        write(*,*) '# of cmp, sta, evt, loc per data'
        write(*,*) inv_struct_cmps, nstations, nevents, nlocations
        write(*,*) '# of dat, mval, xval, size'
        write(*,*) max_data, max_m, max_x, max_size
        write(*,*) '# of row, band, nel'
        write(*,*) max_row, max_band, max_elements
        call checkpoint("retrieved damping data")

        ! allocate matrix
        call csm_init(max_elements, max_row, max_x)
        call checkpoint("allocated G, m, d")

        call load_gmat(gmat_list, smooth_op, inv_struct_cmps, nstations, nevents, &
            nlocations, station_weighting, event_weighting, location_weighting, &
            damping, smoothing, kernel_threshold, data_threshold, max_x, max_m, &
            max_band, sm_max_row, sm_max_band, sm_max_elements)

        call call_cgls(rank, nproc, max_x)

        call csm_free()
    end subroutine main

    subroutine get_input(mode, gmat_list, smooth_op, inv_struct_cmps, &
        nstations, nevents, nlocations, station_weighting, event_weighting, &
        location_weighting, damping, smoothing, kernel_threshold, &
        data_threshold)

        ! mode
        integer(I64) :: choice
        character(len = 10), intent(out) :: mode
        character(len = 10) :: mode0(4)

        ! parameters
        character(len = 132), intent(out) :: gmat_list, smooth_op
        integer(I64), intent(out) :: inv_struct_cmps, nstations, nevents, &
            nlocations
        real(DP), intent(out) :: station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold

        write(*, *) hline

        write(*, "(A, /)") "Mode: "//NEW_LINE('a')// &
            "1) Simple damping"
        read(*, *) choice
        mode = mode0(choice)

        if (choice /= 1) call throw("only simple damping valid")

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

    subroutine get_smoothing(smooth_op, inv_struct_cmps, max_row, max_band, &
        max_elements)
        character(len=*), intent(in) :: smooth_op
        integer(I64), intent(in) :: inv_struct_cmps
        integer(I64), intent(out) :: max_row, max_band, max_elements
        integer :: fid = 10000, nblock, ierr

        max_row = 0
        max_band = 0
        max_elements = 0

        open(fid, file=trim(smooth_op), status="old", iostat=ierr)
        if (ierr /= 0) &
            call throw("smooth_op open err: "//trim(smooth_op))
        read(fid, *) nblock, max_row, max_band, max_elements
        close(fid)

        max_band = max_band * inv_struct_cmps
        max_row = max_row * inv_struct_cmps
        max_elements = max_elements * inv_struct_cmps
    end subroutine get_smoothing

    subroutine get_damping(gmat_list, inv_struct_cmps, nstations, &
        nevents, nlocations, max_data, max_station, max_event, max_location, &
        max_m, max_x, max_size, max_row, max_band, max_elements)

        character(len=*), intent(in) :: gmat_list
        integer(I64), intent(in) :: inv_struct_cmps, nstations, nevents, nlocations
        integer(I64), intent(out) :: max_data, max_station, max_event, &
            max_location, max_m, max_x, max_size, max_row, max_band, &
            max_elements

        character(len=132) :: token
        character(len=300) :: line
        integer :: stations(nstations), events(nevents), locations(nlocations)
        real :: weight(nlocations)
        integer :: fid, kid, nlines, i, j, nblocks, ncoefficients, ierr
        integer(I64) :: nelements
        real(SP) :: kernel
        integer(I64) :: icount = 0

        real(DP) :: rhs, weight_data

        max_data = 0
        max_station = 0
        max_event = 0
        max_location = 0
        max_m = 0
        max_x = 0
        max_size = 0
        max_row = 0
        max_band = 0
        max_elements = 0

        fid = 10000; kid = 10001
        open(fid, file = trim(gmat_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("failed to open gmat_list")

        nlines = 0
        do
            read(fid, "(A300)", iostat = ierr) line
            if (ierr < 0) exit

            ! skip comment lines
            line = adjustl(line)
            if (line(1:1) == '#') cycle
            call advance_to(" ", line, token)

            nlines = nlines + 1
            max_data = max_data + 1

            nelements = 0
            do i = 1, inv_struct_cmps
                ! token currently holds kernel filename
                open(kid, file = trim(token), status = "old", &
                    iostat = ierr, form = "unformatted")
                if (ierr /= 0) &
                    call throw("open kernel err: "//trim(token))

                ! read data from kernel file
                read(kid) nblocks, kernel
                read(kid) ncoefficients

                close(kid)

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

            max_band = max(max_band, nelements)
            max_elements = max_elements + nelements

            if (modulo(nlines, 1000) == 0) &
                print "('Read through line ', I0)", nlines
        end do

        close(fid)

        ! compute final values
        max_m = nblocks * inv_struct_cmps
        max_x = max_m + max_station + max_event + max_location
        max_elements = max_elements + (nstations + nevents + nlocations) * &
            max_data + max_x
        max_row = max_data + max_x
        max_band = max_band + nstations + nevents + nlocations
        max_size = max(max_x, max_row)
    end subroutine get_damping

    subroutine load_gmat(gmat_list, smooth_op, inv_struct_cmps, nstations, nevents, &
        nlocations, station_weighting, event_weighting, location_weighting, &
        damping, smoothing, kernel_threshold, data_threshold, max_x, max_m, max_band, &
        sm_max_row, sm_max_band, sm_max_elements)

        character(len=*), intent(in) :: gmat_list, smooth_op
        integer(I64), intent(in) :: inv_struct_cmps, nstations, nevents, nlocations
        real(DP), intent(in) :: station_weighting, event_weighting, &
            location_weighting, damping, smoothing, kernel_threshold, &
            data_threshold
        integer(I64), intent(in) :: max_x, max_m, max_band
        integer(I64), intent(inout) :: sm_max_row, sm_max_band, sm_max_elements

        character(len=132) :: token
        character(len=300) :: line
        integer :: stations(nstations), events(nevents), locations(nlocations)
        real(DP) :: wglocations(nlocations)
        real(DP) :: rhs_min, rhs_max, rhs_min_weight, rhs_max_weight
        integer :: idx0(max_band), idx(max_band)
        real :: weight(nlocations)
        integer :: fid, kid, nlines, i, j, k, m, nblocks, ncoefficients, ierr
        integer(I64) :: nelements
        real(SP) :: kernel
        integer(I64) :: icount = 0
        real(DP) :: coef_dp(max_band)
        real(SP) :: coef_sp(max_band)
        real(DP) :: rhs, weight_data, tmp

        fid = 10000; kid = 10001
        open(fid, file = trim(gmat_list), status = "old", iostat = ierr)
        if (ierr /= 0) &
            call throw("failed to open gmat_list")

        nlines = 0
        do
            read(fid, "(A300)", iostat = ierr) line
            if (ierr < 0) exit

            ! skip comment lines
            line = adjustl(line)
            if (line(1:1) == '#') cycle
            call advance_to(" ", line, token)

            nlines = nlines + 1

            nelements = 0
            do i = 1, inv_struct_cmps
                ! token currently holds kernel filename
                open(kid, file = trim(token), status = "old", &
                    iostat = ierr, form = "unformatted")
                if (ierr /= 0) &
                    call throw("open kernel err: "//trim(token))

                ! read data from kernel file
                read(kid) nblocks, kernel
                read(kid) ncoefficients
                read(kid) (idx0(nelements+j), coef_sp(nelements+j), j=1, ncoefficients)

                close(kid)

                !$omp parallel do shared(nblocks) reduction(+:idx0)
                    do j = 1, ncoefficients
                        idx0(nelements + j) = idx0(nelements + j) + (nblocks * (i - 1))
                    end do
                !$omp end parallel do

                nelements = nelements + ncoefficients

                ! advance to next token
                call advance_to(" ", line, token)
            end do

            read(line, *) rhs, weight_data
            call advance_to(" ", line, token)
            call advance_to(" ", line, token)

            ! check measurement
            if (abs(rhs) > data_threshold) then
                write(*, "('omitting large rhs: ', F0.4)") rhs
                cycle
            end if

            do i = 1, nstations
                read(line, *) stations(i)
                call advance_to(" ", line, token)
            end do

            do i = 1, nevents
                read(line, *) events(i)
                call advance_to(" ", line, token)
            end do

            do i = 1, nlocations
                read(line, *) locations(i), wglocations(i)
                call advance_to(" ", line, token)
                call advance_to(" ", line, token)
            end do

            ! check kernel value
            m = 0
            do i = 1, nelements
                j = idx0(i)
                tmp = real(coef_sp(i), DP)
                if (dabs(tmp) > kernel_threshold)then
                    m = m + 1
                    idx(m) = j
                    coef_dp(m) = tmp
                end if
            end do

            ! station term
            do i = 1, nstations
                m = m + 1
                idx(m) = max_m + abs(stations(i))
                coef_dp(m) = sign(station_weighting, real(stations(i), DP))
            end do

            ! event term
            do i = 1, nevents
                m = m + 1
                idx(m) = max_m + nstations + abs(events(i))
                coef_dp(m) = sign(event_weighting, real(events(i), DP))
            end do

            ! location term
            do i = 1, nlocations
                m = m + 1
                idx(m) = max_m + nstations + nevents + abs(locations(i))
                coef_dp(m) = wglocations(i) * location_weighting
            end do

            rhs_min = min(rhs_min,rhs)
            rhs_max = max(rhs_max,rhs)

            ! covariance matrix
            !$omp parallel do shared(weight_data) reduction(*:coef_dp)
                do i = 1, m
                    coef_dp(i) = coef_dp(i) * weight_data
                end do
            !$omp end parallel do

            rhs = rhs * weight_data
            rhs_min_weight = min(rhs_min_weight, rhs)
            rhs_max_weight = max(rhs_max_weight, rhs)

            call csm_insert_row(m, idx, coef_dp, real(rhs, DP))

            if (modulo(nlines, 1000) == 0) &
                print "('Read through line ', I0)", nlines
        end do

        close(fid)

        write(*,*) "rhs_min, rhs_max"
        write(*,*) rhs_min, rhs_max
        write(*,*) "rhs_min_weight, rhs_max_weight"
        write(*,*) rhs_min_weight, rhs_max_weight

        call checkpoint("loaded G, d")

        ! apply damping
        do i = 1, max_x
            coef_dp(1) = damping
            idx(1) = i
            call csm_insert_row(1, idx, coef_dp, 0.0_DP)
        end do

        call checkpoint("applied damping")

        ! apply smoothing
        open(fid, file=trim(smooth_op), status="old", iostat=ierr)
        if (ierr /= 0) call throw("open smoothing error: "//trim(smooth_op))
        read(fid, *) nblocks, sm_max_row, sm_max_band, sm_max_elements
        do i = 1, sm_max_row
            read(fid, *) ncoefficients, &
                (idx0(j), j=1, ncoefficients), &
                (coef_sp(j), j=1, ncoefficients)

            do j = 1, inv_struct_cmps
                do k = 1, ncoefficients
                    idx(k) = idx0(k) + (j - 1) * inv_struct_cmps
                    coef_dp(k) = coef_sp(k) * smoothing
                end do
                call csm_insert_row(ncoefficients, idx, coef_dp, 0.0_DP)
            end do
        end do

        close(fid)

        call checkpoint("applied smoothing")
    end subroutine load_gmat

    subroutine call_cgls(rank, nproc, xlen)
        integer, intent(in) :: rank, nproc
        integer(I64), intent(in) :: xlen
        integer :: ierr, itcount
        integer :: fd = 10000
        real(DP), allocatable :: x(:)

        ! open output file
        if (rank == 0) then
            open(fd, file = "output.txt", action="write", iostat=ierr)
            if (ierr /= 0) call throw("failed to open output file")
        end if

        write(*, "('Starting PLSQR on rank ', I0)") rank

        ! wait for threads
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        if (ierr /= 0) call throw("mpi failed to bar")

        allocate(x(xlen))
        call csm_cgls(rank, nproc, x, itcount)

        ! wait for threads
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        if (ierr /= 0) call throw("mpi failed to bar")

        ! close output file
        if (rank == 0) then
            close(fd)
        end if
    end subroutine call_cgls

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

    subroutine proc_info(rank)
        integer, intent(in) :: rank
        character(len=1024) :: hostname

        call hostnm(hostname)

        !write(*, "('rank=', I0, ', pid=', I0, ' on ', A)") rank, getpid(), trim(hostname)
    end subroutine proc_info

    ! information checkpoint
    subroutine checkpoint(msg)
        character(len=*), intent(in) :: msg

        write(*, "('CHECKPOINT: ', A, /, A)") trim(msg), hline
    end subroutine checkpoint

    ! throw an error, exit program
    subroutine throw(emsg)
        character(len=*), intent(in) :: emsg

        write(*, "('Fatal: ', A)") trim(emsg)
        stop 1
    end subroutine throw

end program driver
