program test_neo_bspline_3d
    use libneo_kinds, only: dp
    use math_constants
    use libneo_util, only: linspace
    use neo_bspline

    implicit none

    real(dp), parameter :: X_MIN = 1.23d0
    real(dp), parameter :: X_MAX = TWOPI + 1.23d0
    real(dp), parameter :: TOL_L2_3D = 2.0d-3

    call test_bspline_3d_lsq_full_grid()
    call test_bspline_3d_lsq_batch()

contains

    subroutine test_bspline_3d_lsq_full_grid()
        type(bspline_3d) :: spl
        integer, parameter :: DEGREE(3) = [3, 3, 3]
        integer, parameter :: N_CTRL(3) = [24, 20, 16]
        integer, parameter :: N_GRID1 = 18
        integer, parameter :: N_GRID2 = 16
        integer, parameter :: N_GRID3 = 14

        real(dp) :: x1(N_GRID1), x2(N_GRID2), x3(N_GRID3)
        real(dp) :: f_grid(N_GRID1, N_GRID2, N_GRID3)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2), N_CTRL(3))

        real(dp) :: x(3), f_true, f_fit
        real(dp) :: err2
        integer :: i1, i2, i3

        print *, "Testing neo_bspline 3D LSQ CGLS (full grid)"

        call bspline_3d_init_uniform(spl, DEGREE, N_CTRL, &
            [X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)
        call linspace(X_MIN, X_MAX, N_GRID3, x3)

        do i3 = 1, N_GRID3
            do i2 = 1, N_GRID2
                do i1 = 1, N_GRID1
                    f_grid(i1, i2, i3) = cos(x1(i1))*cos(2.0d0*x2(i2)) &
                        *cos(3.0d0*x3(i3))
                end do
            end do
        end do

        coeff = 0.0_dp
        call bspline_3d_lsq_cgls(spl, x1, x2, x3, f_grid, coeff, &
            max_iter=800, tol=1.0d-10)

        err2 = 0.0d0
        open(unit=30, file="bspline_3d_lsq_grid.dat", status="replace", &
            action="write")
        do i3 = 1, N_GRID3
            do i2 = 1, N_GRID2
                do i1 = 1, N_GRID1
                    x(1) = x1(i1)
                    x(2) = x2(i2)
                    x(3) = x3(i3)
                    f_true = cos(x(1))*cos(2.0d0*x(2))*cos(3.0d0*x(3))
                    call bspline_3d_eval(spl, coeff, x, f_fit)
                    err2 = err2 + (f_fit - f_true)**2
                    write(30,'(5es24.16)') x(1), x(2), x(3), f_true, f_fit
                end do
            end do
        end do
        close(30)

        err2 = err2/real(N_GRID1*N_GRID2*N_GRID3, dp)
        print *, "  3D LSQ full-grid L2 error =", sqrt(err2)

        if (sqrt(err2) > TOL_L2_3D) then
            error stop "neo_bspline 3D LSQ full-grid error too large"
        end if

    end subroutine test_bspline_3d_lsq_full_grid

    subroutine test_bspline_3d_lsq_batch()
        type(bspline_3d) :: spl
        integer, parameter :: DEGREE(3) = [3, 3, 3]
        integer, parameter :: N_CTRL(3) = [24, 20, 16]
        integer, parameter :: N_GRID1 = 18
        integer, parameter :: N_GRID2 = 16
        integer, parameter :: N_GRID3 = 14
        integer, parameter :: N_RHS = 2

        real(dp) :: x1(N_GRID1), x2(N_GRID2), x3(N_GRID3)
        real(dp) :: f_grid(N_GRID1, N_GRID2, N_GRID3)
        real(dp) :: coeff(N_CTRL(1), N_CTRL(2), N_CTRL(3), N_RHS)

        real(dp) :: x(3), f_true(N_RHS), f_fit
        real(dp) :: err2(N_RHS)
        integer :: i1, i2, i3, k

        print *, "Testing neo_bspline 3D LSQ CGLS (batched)"

        call bspline_3d_init_uniform(spl, DEGREE, N_CTRL, &
            [X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX])

        call linspace(X_MIN, X_MAX, N_GRID1, x1)
        call linspace(X_MIN, X_MAX, N_GRID2, x2)
        call linspace(X_MIN, X_MAX, N_GRID3, x3)

        do k = 1, N_RHS
            do i3 = 1, N_GRID3
                do i2 = 1, N_GRID2
                    do i1 = 1, N_GRID1
                        if (k == 1) then
                            f_grid(i1, i2, i3) = cos(x1(i1))*cos(2.0d0*x2(i2)) &
                                *cos(3.0d0*x3(i3))
                        else
                            f_grid(i1, i2, i3) = sin(x1(i1))*cos(2.0d0*x2(i2)) &
                                *cos(1.5d0*x3(i3))
                        end if
                    end do
                end do
            end do
            coeff(:, :, :, k) = 0.0_dp
            call bspline_3d_lsq_cgls(spl, x1, x2, x3, f_grid, &
                coeff(:, :, :, k), max_iter=800, tol=1.0d-10)
        end do

        err2 = 0.0d0
        do i3 = 1, N_GRID3
            do i2 = 1, N_GRID2
                do i1 = 1, N_GRID1
                    x(1) = x1(i1)
                    x(2) = x2(i2)
                    x(3) = x3(i3)
                    f_true(1) = cos(x(1))*cos(2.0d0*x(2))*cos(3.0d0*x(3))
                    f_true(2) = sin(x(1))*cos(2.0d0*x(2))*cos(1.5d0*x(3))
                    do k = 1, N_RHS
                        call bspline_3d_eval(spl, coeff(:, :, :, k), x, f_fit)
                        err2(k) = err2(k) + (f_fit - f_true(k))**2
                    end do
                end do
            end do
        end do

        err2 = err2/real(N_GRID1*N_GRID2*N_GRID3, dp)
        do k = 1, N_RHS
            print *, "  3D LSQ full-grid L2 error (rhs =", k, ") =", sqrt(err2(k))
            if (sqrt(err2(k)) > TOL_L2_3D) then
                error stop "neo_bspline 3D batched LSQ full-grid error too large"
            end if
        end do

    end subroutine test_bspline_3d_lsq_batch

end program test_neo_bspline_3d
