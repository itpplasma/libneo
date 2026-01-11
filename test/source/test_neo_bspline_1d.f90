program test_neo_bspline_1d
    use libneo_kinds, only : dp
    use math_constants
    use libneo_util,  only : linspace
    use neo_bspline

    implicit none

    real(dp), parameter :: TOL_L2 = 5.0d-7
    real(dp), parameter :: X_MIN = 1.23d0
    real(dp), parameter :: X_MAX = TWOPI + 1.23d0

    call test_bspline_1d_partition_unity()
    call test_bspline_1d_lsq_cos()
    call test_bspline_1d_lsq_batch()

contains

    subroutine test_bspline_1d_partition_unity()
        type(bspline_1d) :: spl
        integer, parameter :: DEGREE = 3
        integer, parameter :: N_CTRL = 20
        integer, parameter :: N_SAMPLE = 200

        real(dp) :: x_eval(N_SAMPLE)
        real(dp) :: sum_basis, x, temp
        real(dp) :: max_err
        real(dp) :: c(N_CTRL)
        integer :: i, j

        print *, "Testing neo_bspline 1D partition of unity"

        call bspline_1d_init_uniform(spl, DEGREE, N_CTRL, X_MIN, X_MAX)
        call linspace(X_MIN, X_MAX, N_SAMPLE, x_eval)

        max_err = 0.0d0
        do i = 1, N_SAMPLE
            x = x_eval(i)
            sum_basis = 0.0d0
            do j = 1, N_CTRL
                c = 0.0d0
                c(j) = 1.0d0
                call bspline_1d_eval(spl, c, x, temp)
                sum_basis = sum_basis + temp
            end do
            max_err = max(max_err, abs(sum_basis - 1.0d0))
        end do

        print *, "  max |sum_j N_j(x) - 1| =", max_err

        if (max_err > 1.0d-10) then
            error stop "neo_bspline 1D partition of unity violated"
        end if

    end subroutine test_bspline_1d_partition_unity


    subroutine test_bspline_1d_lsq_cos()
        type(bspline_1d) :: spl
        integer, parameter :: DEGREE = 5
        integer, parameter :: N_CTRL = 60
        integer, parameter :: N_DATA = 400
        integer, parameter :: N_PLOT = 201

        real(dp) :: x_data(N_DATA), f_data(N_DATA)
        real(dp) :: coeff(N_CTRL)
        real(dp) :: x_plot(N_PLOT), f_true(N_PLOT), f_fit(N_PLOT)
        real(dp) :: err2
        integer :: i

        print *, "Testing neo_bspline 1D LSQ CGLS (cos)"

        call bspline_1d_init_uniform(spl, DEGREE, N_CTRL, X_MIN, X_MAX)

        ! Scattered data in [X_MIN, X_MAX]
        do i = 1, N_DATA
            x_data(i) = X_MIN + (X_MAX - X_MIN) * real(i, dp) / real(N_DATA + 1, dp)
            x_data(i) = x_data(i) + 0.1d0 * sin(3.0d0*real(i, dp))
            if (x_data(i) < X_MIN) x_data(i) = X_MIN
            if (x_data(i) > X_MAX) x_data(i) = X_MAX
            f_data(i) = cos(2.0d0*x_data(i)) + 0.5d0*sin(3.0d0*x_data(i))
        end do

        coeff = 0.0_dp
        call bspline_1d_lsq_cgls(spl, x_data, f_data, coeff, max_iter=400, tol=1.0d-10)

        call linspace(X_MIN, X_MAX, N_PLOT, x_plot)
        do i = 1, N_PLOT
            f_true(i) = cos(2.0d0*x_plot(i)) + 0.5d0*sin(3.0d0*x_plot(i))
            call bspline_1d_eval(spl, coeff, x_plot(i), f_fit(i))
        end do

        err2 = sum((f_fit - f_true)**2)/real(N_PLOT, dp)
        print *, "  LSQ L2 error =", sqrt(err2)

        if (sqrt(err2) > TOL_L2) then
            error stop "neo_bspline 1D LSQ error too large"
        end if

        call write_plot_data("bspline_1d_lsq.dat", N_PLOT, x_plot, f_true, f_fit)

    end subroutine test_bspline_1d_lsq_cos


    subroutine test_bspline_1d_lsq_batch()
        type(bspline_1d) :: spl
        integer, parameter :: DEGREE = 5
        integer, parameter :: N_CTRL = 60
        integer, parameter :: N_DATA = 400
        integer, parameter :: N_PLOT = 201
        integer, parameter :: N_RHS = 2

        real(dp) :: x_data(N_DATA)
        real(dp) :: f_data(N_DATA)
        real(dp) :: coeff(N_CTRL, N_RHS)
        real(dp) :: x_plot(N_PLOT)
        real(dp) :: f_true(N_PLOT, N_RHS)
        real(dp) :: f_fit(N_PLOT, N_RHS)
        real(dp) :: err2(N_RHS)
        integer :: i, k

        print *, "Testing neo_bspline 1D LSQ CGLS (batched)"

        call bspline_1d_init_uniform(spl, DEGREE, N_CTRL, X_MIN, X_MAX)

        do i = 1, N_DATA
            x_data(i) = X_MIN + (X_MAX - X_MIN)*real(i, dp)/real(N_DATA + 1, dp)
            x_data(i) = x_data(i) + 0.1d0*sin(3.0d0*real(i, dp))
            if (x_data(i) < X_MIN) x_data(i) = X_MIN
            if (x_data(i) > X_MAX) x_data(i) = X_MAX
        end do

        do k = 1, N_RHS
            do i = 1, N_DATA
                if (k == 1) then
                    f_data(i) = cos(2.0d0*x_data(i)) + 0.5d0*sin(3.0d0*x_data(i))
                else
                    f_data(i) = sin(1.5d0*x_data(i)) + 0.3d0*cos(4.0d0*x_data(i))
                end if
            end do
            coeff(:, k) = 0.0_dp
            call bspline_1d_lsq_cgls(spl, x_data, f_data, coeff(:, k), &
                max_iter=400, tol=1.0d-10)
        end do

        call linspace(X_MIN, X_MAX, N_PLOT, x_plot)
        do i = 1, N_PLOT
            f_true(i, 1) = cos(2.0d0*x_plot(i)) + 0.5d0*sin(3.0d0*x_plot(i))
            f_true(i, 2) = sin(1.5d0*x_plot(i)) + 0.3d0*cos(4.0d0*x_plot(i))
            do k = 1, N_RHS
                call bspline_1d_eval(spl, coeff(:, k), x_plot(i), f_fit(i, k))
            end do
        end do

        do k = 1, N_RHS
            err2(k) = sum((f_fit(:, k) - f_true(:, k))**2)/real(N_PLOT, dp)
            print *, "  LSQ L2 error (rhs =", k, ") =", sqrt(err2(k))
            if (sqrt(err2(k)) > TOL_L2) then
                error stop "neo_bspline 1D batched LSQ error too large"
            end if
        end do

    end subroutine test_bspline_1d_lsq_batch


    subroutine write_plot_data(fname, n, x, f_true, f_fit)
        character(*), intent(in) :: fname
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), f_true(n), f_fit(n)
        integer :: iunit, i

        open(newunit=iunit, file=fname, status="replace", action="write", &
            form="formatted")
        do i = 1, n
            write(iunit,'(3es24.16)') x(i), f_true(i), f_fit(i)
        end do
        close(iunit)
    end subroutine write_plot_data

end program test_neo_bspline_1d
