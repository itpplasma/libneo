program test_neo_bspline_2d
    use libneo_kinds, only : dp
    use math_constants
    use libneo_util,  only : linspace
    use neo_bspline

    implicit none

    real(dp), parameter :: X_MIN = 1.23d0
    real(dp), parameter :: X_MAX = TWOPI + 1.23d0
    real(dp), parameter :: TOL_L2_FULL = 5.0d-4
    real(dp), parameter :: TOL_MAX_CIRC = 1.0d-3

    call test_bspline_2d_adjoint()
    call test_bspline_2d_lsq_full_grid()
    call test_bspline_2d_lsq_circle()
    call test_bspline_2d_lsq_full_grid_batch()

contains

    subroutine test_bspline_2d_adjoint()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [16, 14]
        integer, parameter :: N_GRID1 = 24
        integer, parameter :: N_GRID2 = 20

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: x_data(N_GRID1*N_GRID2)
        real(dp) :: y_data(N_GRID1*N_GRID2)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2))
        real(dp) :: r(N_GRID1*N_GRID2)
        real(dp) :: f(N_GRID1*N_GRID2)
        real(dp) :: g(N_CTRL(1), N_CTRL(2))
        real(dp) :: lhs, rhs, rel_err, denom
        integer :: i1, i2, idx

        print *, "Testing neo_bspline 2D adjoint consistency"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        idx = 0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                idx = idx + 1
                x_data(idx) = x1(i1)
                y_data(idx) = x2(i2)
                r(idx) = cos(0.37d0*x1(i1))*sin(0.41d0*x2(i2))
            end do
        end do

        ! Populate coeff deterministically
        do i2 = 1, N_CTRL(2)
            do i1 = 1, N_CTRL(1)
                coeff(i1, i2) = sin(0.11d0*real(i1, dp)) + cos(0.07d0*real(i2, dp))
            end do
        end do

        call apply_A2D(spl, x_data, y_data, coeff, f)
        call apply_A2D_T(spl, x_data, y_data, r, g)

        lhs = sum(f*r)
        rhs = sum(coeff*g)
        denom = max(abs(lhs), abs(rhs), 1.0d0)
        rel_err = abs(lhs - rhs)/denom

        print *, "  2D adjoint rel error =", rel_err

        if (rel_err > 1.0d-12) then
            error stop "neo_bspline 2D adjoint inconsistency"
        end if

    end subroutine test_bspline_2d_adjoint


    subroutine test_bspline_2d_lsq_full_grid()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [32, 28]
        integer, parameter :: N_GRID1 = 40
        integer, parameter :: N_GRID2 = 36

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: x_data(N_GRID1*N_GRID2)
        real(dp) :: y_data(N_GRID1*N_GRID2)
        real(dp) :: f_data(N_GRID1*N_GRID2)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2))

        integer :: i1, i2, idx
        integer :: i
        real(dp) :: x(2), f_true, f_fit
        real(dp) :: err2

        print *, "Testing neo_bspline 2D LSQ CGLS (full grid)"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        idx = 0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                idx = idx + 1
                x_data(idx) = x1(i1)
                y_data(idx) = x2(i2)
                f_data(idx) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        coeff = 0.0_dp
        call bspline_2d_lsq_cgls(spl, x_data, y_data, f_data, coeff, &
            max_iter=800, tol=1.0d-10)

        err2 = 0.0d0
        idx = 0
        open(unit=20, file="bspline_2d_lsq_grid.dat", status="replace", action="write")
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                idx = idx + 1
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


    subroutine test_bspline_2d_lsq_circle()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [32, 28]
        integer, parameter :: N_GRID1 = 40
        integer, parameter :: N_GRID2 = 36
        integer, parameter :: MAX_DATA = N_GRID1*N_GRID2

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: x_data(MAX_DATA), y_data(MAX_DATA), f_data(MAX_DATA)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2))

        real(dp) :: xc, yc, r0, dx, dy
        real(dp) :: theta, r_eval
        real(dp) :: x(2), f_true, f_fit
        real(dp) :: max_err
        integer :: i1, i2, ndata, i

        print *, "Testing neo_bspline 2D LSQ CGLS (circular mask)"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        xc = 0.5d0*(X_MIN + X_MAX)
        yc = 0.5d0*(X_MIN + X_MAX)
        r0 = 0.45d0*(X_MAX - X_MIN)

        ndata = 0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                dx = x1(i1) - xc
                dy = x2(i2) - yc
                if (dx*dx + dy*dy <= r0*r0) then
                    ndata = ndata + 1
                    x_data(ndata) = x1(i1)
                    y_data(ndata) = x2(i2)
                    f_data(ndata) = cos(x1(i1))*cos(2.0d0*x2(i2))
                end if
            end do
        end do

        if (ndata < 100) then
            error stop "neo_bspline 2D LSQ circle: too few data points"
        end if

        coeff = 0.0_dp
        call bspline_2d_lsq_cgls(spl, x_data(1:ndata), y_data(1:ndata), &
            f_data(1:ndata), coeff, max_iter=800, tol=1.0d-10)

        max_err = 0.0d0
        r_eval = 0.4d0*(X_MAX - X_MIN)

        open(unit=21, file="bspline_2d_lsq_circle.dat", status="replace", &
            action="write")
        do i = 1, 64
            theta = TWOPI*real(i, dp)/64.0d0
            x(1) = xc + r_eval*cos(theta)
            x(2) = yc + r_eval*sin(theta)
            f_true = cos(x(1))*cos(2.0d0*x(2))
            call bspline_2d_eval(spl, coeff, x, f_fit)
            max_err = max(max_err, abs(f_fit - f_true))
            write(21,'(5es24.16)') theta, x(1), x(2), f_true, f_fit
        end do
        close(21)

        print *, "  2D LSQ circle: max error =", max_err

        if (max_err > TOL_MAX_CIRC) then
            error stop "neo_bspline 2D LSQ circle error too large"
        end if

    end subroutine test_bspline_2d_lsq_circle


    subroutine test_bspline_2d_lsq_full_grid_batch()
        type(bspline_2d) :: spl
        integer, parameter :: DEGREE(2) = [3, 3]
        integer, parameter :: N_CTRL(2) = [32, 28]
        integer, parameter :: N_GRID1 = 40
        integer, parameter :: N_GRID2 = 36
        integer, parameter :: N_RHS = 2

        real(dp) :: x1(N_GRID1), x2(N_GRID2)
        real(dp) :: x_data(N_GRID1*N_GRID2)
        real(dp) :: y_data(N_GRID1*N_GRID2)
        real(dp) :: f_data(N_GRID1*N_GRID2, N_RHS)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2), N_RHS)

        real(dp) :: x(2), f_true(N_RHS), f_fit
        real(dp) :: err2(N_RHS)
        integer :: i1, i2, idx, k

        print *, "Testing neo_bspline 2D LSQ CGLS (full grid, batched)"

        call bspline_2d_init_uniform(spl, DEGREE, N_CTRL, [X_MIN, X_MIN], &
            [X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)

        idx = 0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                idx = idx + 1
                x_data(idx) = x1(i1)
                y_data(idx) = x2(i2)
                f_data(idx, 1) = cos(x1(i1))*cos(2.0d0*x2(i2))
                f_data(idx, 2) = sin(x1(i1))*cos(3.0d0*x2(i2))
            end do
        end do

        coeff = 0.0_dp
        call bspline_2d_lsq_cgls_batch(spl, x_data, y_data, f_data, coeff, &
            max_iter=800, tol=1.0d-10)

        err2 = 0.0d0
        idx = 0
        do i2 = 1, N_GRID2
            do i1 = 1, N_GRID1
                idx = idx + 1
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
