program trupin
    use, intrinsic :: iso_fortran_env, only : &
        i8 => int8, &
        i16 => int16, &
        i32 => int32, &
        i64 => int64, &
        sp => real32, &
        dp => real64

    implicit none

    ! kernel variables
    character (len = 132) :: kernel_list, kernel_file
    integer :: kernel_list_io, kernel_file_io
    character (len = 300) :: kernel_entry
    integer(i8) :: inv_struct_cmps, &
        nstations, nevents, nlocations, nblocks, nelements, ncoefs
    integer, allocatable :: stations(:), events(:), locations(:)
    real(dp), allocatable :: wglocations(:)
    real(dp) :: measurement, weight
    real(dp) :: kernel_threshold, data_threshold
    real(sp) :: kernel0(3)
    real(sp), allocatable :: coefs_sp(:)
    real(dp), allocatable :: coefs_dp(:)
    integer, allocatable :: index0(:)

    ! general variables
    integer :: i, j, k, n, ierr

    ! read in kernel list
    write(*, *), "File list for G matrix, data, and weighting:"
    read(*, "(a)") kernel_list
    write(*, "('using file ', a)") trim(kernel_list)

    ! inverted structure components
    write(6, *) "# of inverted structure components"
    read(5, *) inv_struct_cmps
    write(*, "('using ', i0, ' components')") inv_struct_cmps

    ! thresholds
    write(6, *) "Kernel and data thresholds (2): "
    read(5, *) kernel_threshold, data_threshold
    write(*, "('using kernel threshold ', f0.2, ' and data threshold ', f0.2)") &
        kernel_threshold, data_threshold

    ! open kernel file
    kernel_list_io = 10001
    kernel_file_io = 10002
    open(unit = kernel_list_io, file = trim(kernel_list), &
         status = "old", iostat = ierr)
    if (ierr /= 0) &
        call throw("kernel_list open err: "//trim(kernel_list))

    write(*, *) "starting loop over "//trim(kernel_list)
    i = 0
    kernel_list_loop: do

        read(kernel_list_io, "(a300)", iostat = ierr) kernel_entry
        if (ierr /= 0) exit
        write(*, "('entry ', i0)") i

        kernel_entry = adjustl(kernel_entry); if (kernel_entry(1:1) == "#") &
            cycle

        i = i + 1

        ! read kernel files
        nelements = 0
        do j = 1, inv_struct_cmps
            n = index(kernel_entry(1:len_trim(kernel_entry)), " ")
            kernel_file = kernel_entry(1:n)
            kernel_entry(1:) = kernel_entry(n+1:)
            kernel_entry = adjustl(kernel_entry)

            open(unit = kernel_file_io, file = trim(kernel_file), &
                status = "old", iostat = ierr, form = "unformatted")
            if (ierr /= 0) call throw("open kernel error: "//trim(kernel_file))
            read(kernel_file_io) nblocks, kernel0(j)
            read(kernel_file_io) ncoefs
            read(kernel_file_io) (index0(nelements + k), coefs_sp(nelements + k), k = 1, ncoefs)
            close(kernel_file_io)

            do k = 1, ncoefs
                index0(nelements + k) = index0(nelements + k) + nblocks * (j - 1)
            end do

            nelements = nelements + ncoefs
        end do

        ! read measurement and weight
        read(kernel_entry, *) measurement, weight
        call advance(kernel_entry, " ", n)
        call advance(kernel_entry, " ", n)

        ! check measurement
        if (abs(measurement) > data_threshold) then
            write(*, "('omit large measurement ', a)") measurement
            cycle
        end if

        ! stations
        do k = 1, nstations
            read(kernel_entry, *) stations(k)
            call advance(kernel_entry, " ", n)
        end do

        ! events
        do k = 1, nevents
            read(kernel_entry, *) events(k)
            call advance(kernel_entry, " ", n)
        end do

        ! locations
        do k = 1, nlocations
            read(kernel_entry, *) locations(k), wglocations(k)
            call advance(kernel_entry, " ", n)
            call advance(kernel_entry, " ", n)
        end do

    end do kernel_list_loop
    write(*, *) "completed loop over "//trim(kernel_list)

    close(kernel_list_io)

contains

    ! advance string to next delimiter
    subroutine advance(string, delimiter, n)
        character (len = *), intent(inout) :: string
        character (len = *), intent(in) :: delimiter
        integer, intent(inout), optional :: n

        n = index(string(1:len_trim(string)), delimiter)
        string(1:) = string(n+1:)
        string = adjustl(string)
    end subroutine advance

    ! print a message, stop execution
    subroutine throw(msg)
        character (len = *), intent(in) :: msg
        print *, msg
        stop
    end subroutine throw

end program trupin
