module CompressedSparseMatrix
    use, intrinsic :: iso_fortran_env, only : &
        I8 => int8, &
        I16 => int16, &
        I32 => int32, &
        I64 => int64, &
        SP => real32, &
        DP => real64

    implicit none

    private

    ! definition of variables
    !     mrow, ncol = system dimensions (mrow x ncol)
    !     nel = number of active elements in matrix
    !     values = coefficient matrix values
    !     col_index = column indices
    !     row_index = pointer to start of row
    !     b = solution vector
    integer(I64) :: mrow, ncol, nel
    real(DP), allocatable :: val(:), b(:)
    integer(I64), allocatable :: col_ind(:), row_ptr(:)
    integer(I64) :: max_row, max_nel, row_index_ptr

    public :: csm_init, csm_output_n_rows, csm_insert_row, csm_free, csm_cgls,
        csm_error, x

contains

    subroutine csm_init(max_nel_, max_row_, ncol_)
        integer, intent(in) :: max_nel_, max_row_, ncol_

        max_nel = max_nel_
        max_row = max_row_

        nel = 0
        mrow = 0
        ncol = ncol_

        ! allocate mem
        allocate(val(max_nel))
        allocate(b(max_row))
        allocate(col_ind(max_nel))
        allocate(row_ptr(max_row + 1))

        val = 0.0
        b = 0.0
        col_ind = 0.0
        row_ptr = 0.0
        row_index_ptr = 1
    end subroutine csm_init

    subroutine csm_insert_row(ncoef, index, coef, rhs)
        ! insert equation into row of matrix
        !     ncoefs = number of active elements
        !     indices = column indices
        !     coefs = coefficients
        integer, intent(in) :: ncoef, index(:)
        real(DP), intent(in) :: coef(:), rhs
        integer :: i

        ! short circuit no coefficients
        if (ncoef <= 0) return

        ! check limits
        if (nel + ncoef > max_nel) then
            write(*, *) mrow + 1, nel, ncoef, max_nel
            call throw("exceed max_nel")
        end if

        ! inc total number of rows
        mrow = mrow + 1

        ! add rhs to b vector
        b(mrow) = rhs

        ! add num of coefs to active elements vector
        row_ptr(mrow) = row_index_ptr
        row_index_ptr = row_index_ptr + ncoef

        ! copy matrix values
        val(nel+1:nel+ncoef) = coef(1:ncoef)

        ! copy column indices
        col_ind(nel+1:nel+ncoef) = index(1:ncoef)

        ! inc total number of elements
        nel = nel + ncoef
        row_ptr(max_row + 1) = nel + 1
    end subroutine csm_insert_row

    subroutine csm_free()
        deallocate(val)
        deallocate(b)
        deallocate(row_ptr)
        deallocate(col_ind)
    end subroutine csm_free

    subroutine csm_cgls(rank, nproc, x, itcount)
        !--------------------------------------------------------------
        ! CONJUGATE GRADIENT ROUTINE FOR NON-SQUARE LEAST SQUARES
        !
        integer, intent(in) :: rank, nproc
        real(DP) :: x(ncol)
        integer, intent(out) :: itcount

        real(DP), parameter :: tol = 1.0e-7_DP

        real(DP) :: alpha, beta, phibar, rhobar, phi, rho, c, s, theta, &
            tmp, anorm, test1, test2, rnorm
        real(DP) :: u(mrow), v(ncol), &
            sig(ncol), q(mrow), &
            w(ncol)
        integer :: i, ierr, maxiter = 30000

        if (rank == 0) then
            write(*, "('Conjugate gradient routine for non-square least squares.')")
        end if

        anorm = dot_product(val, val)
        anorm = dsqrt(anorm)

        ! initial values
        u = b
        call linearize(u, beta)

        ! v = A^T * u
        call tvecmul(u, v)
        call linearize(v, alpha)

        ! w = v
        w = v

        ! zero out x and sig
        x = 0.0_DP
        sig = 0.0_DP

        ! phibar rhobar
        phibar = beta
        rhobar = alpha

        ! and get the party started
        do itcount = 1, maxiter
            call vecmul(v, q)
            call scaleadd(u, q, '-', alpha)
            call linearize(u, beta)

            call tvecmul(u, q)
            call scaleadd(v, q, '-', beta)
            call linearize(v, alpha)

            ! paige and saunders (a) through (g)
            rho = dsqrt(rhobar ** 2 + beta ** 2)

            ! convergence check
            if (rho .le. tol) exit

            c = rhobar / rho
            s = beta / rho
            theta = s * alpha
            rhobar = -c * alpha
            phi = c * phibar
            phibar = s * phibar

            ! update solution vector and w vector

            !$omp parallel do private(tmp) shared(w) reduction(+:x,sig)
            do i = 1, ncol
                tmp = w(i) / rho
                w(i) = v(i) - theta * tmp
                x(i) = x(i) + phi * tmp
                sig(i) = sig(i) + tmp ** 2
            end do
            !$omp end parallel do

            ! prepare for convergence checks
            rnorm = phibar
            test1 = rnorm / anorm / dsqrt(dot_product(x, x))
            test2 = alpha * dabs(c) / anorm

            if ((test1 .le. tol) .or. (test2 .le. tol)) then
                exit
            end if
        end do

        if (rank == 0) then
            if (itcount <= maxiter) then
                print *, "Solution converged in ", itcount, " iterations."
                print "(*(X, F0.2))", x
            else
                print *, "Solution did not converge."
            end if
        end if

    end subroutine csm_cgls

    subroutine linearize(vec, scalar)
        real(DP), intent(inout) :: vec(:)
        real(DP), intent(out) :: scalar
        integer :: i

        scalar = dot_product(vec, vec)
        scalar = dsqrt(scalar)

        do i = 1, size(vec)
            vec(i) = vec(i) / scalar
        end do
    end subroutine linearize

    subroutine scaleadd(a, b, flag, scalar)
        real(DP), intent(inout) :: a(:)
        real(DP), intent(in) :: b(:), scalar
        character(len=*), intent(in) :: flag
        integer :: i

        if (flag(1:1) .eq. '+') then
            do i = 1, size(a)
                a(i) = b(i) + a(i) * scalar
            end do
        else
            do i = 1, size(a)
                a(i) = b(i) - a(i) * scalar
            end do
        end if
    end subroutine scaleadd

    subroutine vecmul(x, y)
        real(DP), intent(in) :: x(:)
        real(DP), intent(out) :: y(:)
        real(DP) :: y0
        integer :: i, j, l, l1, l2

        l2 = 0

        do i = 1, mrow
            y0 = 0.0_DP
            l1 = l2 + 1
            l2 = l2 + (row_ptr(i + 1) - row_ptr(i))
            !$omp parallel do reduction(+:y0)
                do l = l1, l2
                    j = col_ind(l)
                    y0 = y0 + val(l) * x(j)
                end do
            !$omp end parallel do
            y(i) = y0
        end do
    end subroutine vecmul

    subroutine tvecmul(x, y)
        real(DP), intent(in) :: x(:)
        real(DP), intent(out) :: y(:)
        real(DP) :: xi
        integer :: i, j, l, l1, l2

        l2 = 0
        y = 0.0_DP

        do i = 1, mrow
            xi = x(i)
            l1 = l2 + 1
            l2 = l2 + (row_ptr(i + 1) - row_ptr(i))
            !$omp parallel do reduction(+:y)
                do l = l1, l2
                    j = col_ind(l)
                    y(j) = y(j) + val(l) * xi
                end do
            !$omp end parallel do
        end do
    end subroutine tvecmul

    subroutine csm_output_n_rows(nrows)
        ! output_n_rows
        ! output the first n rows of the current matrix
        !     nrows: number of rows to print
        integer, intent(in) :: nrows
        integer :: i, j

        if (nrows > mrow) then
            call throw("invalid row number specified")
        end if

        do i = 1, nrows
            write(*, "(A)", advance="no") "| "
            do j = 1, row_ptr(i + 1) - row_ptr(i)
                write(*, "('[', I0, ']: ', F0.2, ' ')", advance="no") &
                    col_ind(row_ptr(i) + j - 1), val(row_ptr(i) + j - 1)
            end do
            write(*, "(A)") "|"
        end do
    end subroutine csm_output_n_rows

    ! throw an error, exit program
    subroutine throw(emsg)
        character(len=*), intent(in) :: emsg

        write(*, "('Fatal: ', A)") trim(emsg)
        stop 1
    end subroutine throw

end module CompressedSparseMatrix
