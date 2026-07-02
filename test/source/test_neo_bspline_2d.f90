program test_neo_bspline_2d
    use libneo_kinds, only: dp
    use math_constants
    use libneo_util, only: linspace
    use neo_bspline

    implicit none

    real(dp), parameter :: X_MIN = 1.23d0
    real(dp), parameter :: X_MAX = TWOPI + 1.23d0
    real(dp), parameter :: TOL_L2_FULL = 5.0d-4

    call test_bspline_2d_lsq_full_grid()
    call test_bspline_2d_interp_exact()
    call test_bspline_2d_lsq_full_grid_batch()

contains

    subroutine test_bspline_2d_lsq_full_grid()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [32, 28]
        integer, parameter :: N_GRID1 = 40
        integer, parameter :: N_GRID2 = 36

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: f_grid(N_GRID1, N_GRID2)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2))

        integer :: i1, i2
        real(dp) :: x(2), f_true, f_fit
        real(dp) :: err2

        print *, "Testing neo_bspline 2D LSQ CGLS (full grid)"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                f_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        coeff = 0.0_dp
        call bspline_2d_lsq_cgls(spl, x1, x2, f_grid, coeff, &
            max_iter=800, tol=1.0d-10)

        err2 = 0.0d0
        open(unit=20, file="bspline_2d_lsq_grid.dat", status="replace", &
            action="write")
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                x(1) = x1(i1)
                x(2) = x2(i2)
                f_true = cos(x(1))*cos(2.0d0*x(2))
                call bspline_2d_eval(spl, coeff, x, f_fit)
                err2 = err2 + (f_fit - f_true)**2
                write(20,'(4es24.16)') x(1), x(2), f_true, f_fit
            end do
        end do
        close(20)

        err2 = err2/real(N_GRID1*N_GRID2, dp)
        print *, "  2D LSQ full-grid L2 error =", sqrt(err2)

        if (sqrt(err2) > TOL_L2_FULL) then
            error stop "neo_bspline 2D LSQ full-grid error too large"
        end if

    end subroutine test_bspline_2d_lsq_full_grid


    subroutine test_bspline_2d_interp_exact()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N1 = 18
        integer, parameter :: N2 = 14

        real(dp) :: x1(N1), x2(N2)
        real(dp) :: f_grid(N1, N2)
        real(dp) :: coeff(N1, N2)
        real(dp) :: f_fit, f_true, x(2)
        real(dp) :: err_max
        integer :: i1, i2

        print *, "Testing neo_bspline 2D direct interpolation (collocation)"

        call bspline_2d_init_uniform(spl, DEGREE, [N1, N2], [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N1, x1)
        call linspace(X_MIN, X_MAX, N2, x2)

        do i2 = 1, N2
            do i1 = 1, N1
                f_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        call bspline_2d_interp(spl, x1, x2, f_grid, coeff)

        err_max = 0.0d0
        do i2 = 1, N2
            do i1 = 1, N1
                x = [x1(i1), x2(i2)]
                f_true = cos(x(1))*cos(2.0d0*x(2))
                call bspline_2d_eval(spl, coeff, x, f_fit)
                err_max = max(err_max, abs(f_fit - f_true))
            end do
        end do

        print *, "  2D direct interp max error =", err_max
        if (err_max > 1.0d-10) then
            error stop "neo_bspline 2D direct interpolation error too large"
        end if
    end subroutine test_bspline_2d_interp_exact


    subroutine test_bspline_2d_lsq_full_grid_batch()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [32, 28]
        integer, parameter :: N_GRID1 = 40
        integer, parameter :: N_GRID2 = 36
        integer, parameter :: N_RHS = 2

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: f_grid(N_GRID1, N_GRID2)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2), N_RHS)

        real(dp) :: x(2), f_true(N_RHS), f_fit
        real(dp) :: err2(N_RHS)
        integer :: i1, i2, k

        print *, "Testing neo_bspline 2D LSQ CGLS (full grid, batched)"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        do k = 1, N_RHS
            do i2 = 1, N_GRID2
                do i1 = 1, N_GRID1
                    if (k == 1) then
                        f_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
                    else
                        f_grid(i1, i2) = sin(x1(i1))*cos(3.0d0*x2(i2))
                    end if
                end do
            end do
            coeff(:, :, k) = 0.0_dp
            call bspline_2d_lsq_cgls(spl, x1, x2, f_grid, coeff(:, :, k), &
                max_iter=800, tol=1.0d-10)
        end do

        err2 = 0.0d0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                x(1) = x1(i1)
                x(2) = x2(i2)
                f_true(1) = cos(x(1))*cos(2.0d0*x(2))
                f_true(2) = sin(x(1))*cos(3.0d0*x(2))
                do k = 1, N_RHS
                    call bspline_2d_eval(spl, coeff(:, :, k), x, f_fit)
                    err2(k) = err2(k) + (f_fit - f_true(k))**2
                end do
            end do
        end do

        err2 = err2/real(N_GRID1*N_GRID2, dp)
        do k = 1, N_RHS
            print *, "  2D LSQ full-grid L2 error (rhs =", k, ") =", sqrt(err2(k))
            if (sqrt(err2(k)) > TOL_L2_FULL) then
                error stop "neo_bspline 2D batched LSQ full-grid error too large"
            end if
        end do

    end subroutine test_bspline_2d_lsq_full_grid_batch

end program test_neo_bspline_2d
