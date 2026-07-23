program test_spline_cgls_global
    use libneo_kinds, only : dp
    use math_constants
    use libneo_util,  only : linspace
    use interpolate, only : SplineData2D, construct_splines_2d, &
        evaluate_splines_2d, build_design_matrix_2d
    use cgls_dense,  only : cgls_dense_solve

    implicit none

    call test_design_matrix_2d_exact()

contains

    subroutine test_design_matrix_2d_exact()
        integer, parameter :: N_POINTS1 = 16
        integer, parameter :: N_POINTS2 = 14
        integer, parameter :: ORDER1 = 3
        integer, parameter :: ORDER2 = 3

        real(dp), parameter :: TOL = 1.0d-10

        real(dp) :: x_min(2), x_max(2)
        real(dp) :: x1(N_POINTS1), x2(N_POINTS2)
        real(dp) :: y_grid(N_POINTS1, N_POINTS2)
        real(dp) :: x_data(N_POINTS1*N_POINTS2)
        real(dp) :: y_data(N_POINTS1*N_POINTS2)
        real(dp) :: f_data(N_POINTS1*N_POINTS2)
        real(dp) :: phi_mat(N_POINTS1*N_POINTS2, N_POINTS1*N_POINTS2)
        real(dp) :: y_vec(N_POINTS1*N_POINTS2)
        real(dp) :: y_fit(N_POINTS1, N_POINTS2)
        type(SplineData2D) :: spl

        integer :: i1, i2, idx
        real(dp) :: x_eval(2), f_true, f_fit, err, max_err

        print *, "Testing 2D design matrix CGLS exactness (global knots)..."

        x_min(1) = 1.23d0
        x_min(2) = 1.23d0
        x_max(1) = TWOPI + 1.23d0
        x_max(2) = TWOPI + 1.23d0

        call linspace(x_min(1), x_max(1), N_POINTS1, x1)
        call linspace(x_min(2), x_max(2), N_POINTS2, x2)

        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                y_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        call construct_splines_2d(x_min, x_max, y_grid, &
            [ORDER1, ORDER2], [.false., .false.], spl)

        idx = 0
        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                idx = idx + 1
                x_data(idx) = x1(i1)
                y_data(idx) = x2(i2)
                x_eval(1) = x1(i1)
                x_eval(2) = x2(i2)
                call evaluate_splines_2d(spl, x_eval, f_data(idx))
            end do
        end do

        call build_design_matrix_2d(x_min, x_max, [ORDER1, ORDER2], &
            [.false., .false.], [N_POINTS1, N_POINTS2], x_data, y_data, phi_mat)

        y_vec = 0.0d0
        call cgls_dense_solve(phi_mat, f_data, y_vec, max_iter=N_POINTS1*N_POINTS2)

        y_fit = reshape(y_vec, shape=[N_POINTS1, N_POINTS2])

        max_err = 0.0d0
        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                x_eval(1) = x1(i1)
                x_eval(2) = x2(i2)
                f_true = cos(x_eval(1))*cos(2.0d0*x_eval(2))
                call evaluate_splines_2d(spl, x_eval, f_fit)
                err = abs(f_fit - f_true)
                if (err > max_err) max_err = err
            end do
        end do

        if (max_err > TOL) then
            print *, "  2D grid spline max error =", max_err
            error stop "2D grid spline does not match analytic in design-matrix test"
        end if

        max_err = 0.0d0
        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                err = abs(y_fit(i1, i2) - y_grid(i1, i2))
                if (err > max_err) max_err = err
            end do
        end do

        if (max_err > TOL) then
            print *, "  2D CGLS design-matrix max knot error =", max_err
            error stop "2D CGLS failed to reconstruct global knots"
        end if

        print *, "  PASSED: 2D global-knot design matrix CGLS"
    end subroutine test_design_matrix_2d_exact

end program test_spline_cgls_global
