program test_batch_eval_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_util, only: linspace
    use interpolate, only: SplineData1D, SplineData2D, SplineData3D, &
                           construct_splines_1d, construct_splines_2d, &
                           construct_splines_3d, &
                           destroy_splines_1d, destroy_splines_2d, &
                           destroy_splines_3d, &
                           evaluate_splines_1d, evaluate_splines_1d_der, &
                           evaluate_splines_1d_der2, &
                           evaluate_splines_2d, evaluate_splines_2d_der, &
                           evaluate_splines_3d, evaluate_splines_3d_der, &
                           evaluate_splines_3d_der2
    use spline_eval_reference, only: evaluate_splines_1d_ref, &
                                     evaluate_splines_1d_der_ref, &
                                     evaluate_splines_1d_der2_ref, &
                                     evaluate_splines_2d_ref, &
                                     evaluate_splines_2d_der_ref, &
                                     evaluate_splines_3d_ref, &
                                     evaluate_splines_3d_der_ref, &
                                     evaluate_splines_3d_der2_ref

    implicit none

    real(dp), parameter :: REL_TOL = 1.0d-11

    call test_1d_evaluation()
    call test_2d_evaluation()
    call test_3d_evaluation()

    print *, "All facade vs reference oracle tests PASSED!"

contains

    subroutine test_1d_evaluation()
        integer, parameter :: n_points = 50, n_eval = 100
        real(dp), parameter :: x_min = 1.0d0, x_max = 5.0d0
        integer :: order, ie, iper
        logical :: periodic
        real(dp) :: x_grid(n_points), y_grid(n_points)
        real(dp) :: x_eval, period
        real(dp) :: y_facade, dy_facade, d2y_facade
        real(dp) :: y_ref, dy_ref, d2y_ref
        type(SplineData1D) :: spl

        call linspace(x_min, x_max, n_points, x_grid)
        do ie = 1, n_points
            y_grid(ie) = exp(0.2d0*cos(x_grid(ie))) * (2.0d0 + sin(x_grid(ie)))
        end do

        period = x_max - x_min
        do order = 3, 5
            do iper = 0, 1
                periodic = (iper == 1)
                call construct_splines_1d(x_min, x_max, y_grid, order, periodic, spl)

                do ie = 1, n_eval
                    if (periodic) then
                        x_eval = x_min + 3.0d0*period*dble(ie)/dble(n_eval+1) - period
                    else
                        x_eval = x_min + period*dble(ie)/dble(n_eval+1)
                    end if

                    call evaluate_splines_1d(spl, x_eval, y_facade)
                    call evaluate_splines_1d_ref(spl, x_eval, y_ref)
                    call assert_close("1d eval", y_facade, y_ref, REL_TOL)

                    call evaluate_splines_1d_der(spl, x_eval, y_facade, dy_facade)
                    call evaluate_splines_1d_der_ref(spl, x_eval, y_ref, dy_ref)
                    call assert_close("1d der y", y_facade, y_ref, REL_TOL)
                    call assert_close("1d der dy", dy_facade, dy_ref, REL_TOL)

                    call evaluate_splines_1d_der2(spl, x_eval, y_facade, dy_facade, d2y_facade)
                    call evaluate_splines_1d_der2_ref(spl, x_eval, y_ref, dy_ref, d2y_ref)
                    call assert_close("1d der2 y", y_facade, y_ref, REL_TOL)
                    call assert_close("1d der2 dy", dy_facade, dy_ref, REL_TOL)
                    call assert_close("1d der2 d2y", d2y_facade, d2y_ref, REL_TOL)
                end do

                call destroy_splines_1d(spl)
            end do
        end do

        print *, "  PASSED: 1D facade vs reference (value + derivatives)"
    end subroutine test_1d_evaluation

    subroutine test_2d_evaluation()
        integer, parameter :: n1 = 25, n2 = 27, n_eval = 50
        real(dp), parameter :: x_min(2) = [1.0d0, -0.5d0]
        real(dp), parameter :: x_max(2) = [4.0d0, 2.5d0]
        integer :: order(2), ie, i1, i2, icase
        logical :: periodic(2)
        real(dp) :: x1(n1), x2(n2), y_grid(n1, n2)
        real(dp) :: x_eval(2), period(2)
        real(dp) :: y_facade, dy_facade(2)
        real(dp) :: y_ref, dy_ref(2)
        type(SplineData2D) :: spl
        integer :: order_cases(2, 3)
        logical :: periodic_cases(2, 3)

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)

        do i2 = 1, n2
            do i1 = 1, n1
                y_grid(i1, i2) = exp(0.15d0*cos(x1(i1))) * (1.5d0 + 0.5d0*sin(x2(i2)))
            end do
        end do

        period = x_max - x_min
        order_cases = reshape([3, 3, 5, 5, 5, 3], [2, 3])
        periodic_cases = reshape([.false., .false., .true., .true., .true., .false.], &
                                 [2, 3])

        do icase = 1, 3
            order = order_cases(:, icase)
            periodic = periodic_cases(:, icase)
            call construct_splines_2d(x_min, x_max, y_grid, order, periodic, spl)

            do ie = 1, n_eval
                x_eval(1) = x_min(1) + period(1)*dble(ie)/dble(n_eval+1)
                x_eval(2) = x_min(2) + period(2)*dble(ie+3)/dble(n_eval+5)

                call evaluate_splines_2d(spl, x_eval, y_facade)
                call evaluate_splines_2d_ref(spl, x_eval, y_ref)
                call assert_close("2d eval", y_facade, y_ref, REL_TOL)

                call evaluate_splines_2d_der(spl, x_eval, y_facade, dy_facade)
                call evaluate_splines_2d_der_ref(spl, x_eval, y_ref, dy_ref)
                call assert_close("2d der y", y_facade, y_ref, REL_TOL)
                call assert_close("2d der dy1", dy_facade(1), dy_ref(1), REL_TOL)
                call assert_close("2d der dy2", dy_facade(2), dy_ref(2), REL_TOL)
            end do

            call destroy_splines_2d(spl)
        end do

        print *, "  PASSED: 2D facade vs reference (value + derivatives)"
    end subroutine test_2d_evaluation

    subroutine test_3d_evaluation()
        integer, parameter :: n1 = 12, n2 = 13, n3 = 14, n_eval = 30
        real(dp), parameter :: x_min(3) = [0.5d0, -1.0d0, 2.0d0]
        real(dp), parameter :: x_max(3) = [2.0d0, 1.0d0, 4.0d0]
        integer :: order(3), ie, i1, i2, i3, icase
        logical :: periodic(3)
        real(dp) :: x1(n1), x2(n2), x3(n3), y_grid(n1, n2, n3)
        real(dp) :: x_eval(3), period(3)
        real(dp) :: y_facade, dy_facade(3), d2y_facade(6)
        real(dp) :: y_ref, dy_ref(3), d2y_ref(6)
        type(SplineData3D) :: spl
        integer :: order_cases(3, 2)
        logical :: periodic_cases(3, 2)

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)
        call linspace(x_min(3), x_max(3), n3, x3)

        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    y_grid(i1, i2, i3) = exp(0.1d0*cos(x1(i1))) * &
                                         (1.2d0 + 0.3d0*sin(x2(i2))) * &
                                         (1.4d0 + 0.2d0*cos(x3(i3)))
                end do
            end do
        end do

        period = x_max - x_min
        order_cases = reshape([5, 3, 3, 3, 5, 5], [3, 2])
        periodic_cases = reshape([.true., .true., .true., .false., .true., .false.], &
                                 [3, 2])

        do icase = 1, 2
            order = order_cases(:, icase)
            periodic = periodic_cases(:, icase)
            call construct_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

            do ie = 1, n_eval
                x_eval(1) = x_min(1) + period(1)*dble(ie)/dble(n_eval+1)
                x_eval(2) = x_min(2) + period(2)*dble(ie+2)/dble(n_eval+3)
                x_eval(3) = x_min(3) + period(3)*dble(ie+5)/dble(n_eval+7)

                call evaluate_splines_3d(spl, x_eval, y_facade)
                call evaluate_splines_3d_ref(spl, x_eval, y_ref)
                call assert_close("3d eval", y_facade, y_ref, REL_TOL)

                call evaluate_splines_3d_der(spl, x_eval, y_facade, dy_facade)
                call evaluate_splines_3d_der_ref(spl, x_eval, y_ref, dy_ref)
                call assert_close("3d der y", y_facade, y_ref, REL_TOL)
                call assert_close("3d der dy1", dy_facade(1), dy_ref(1), REL_TOL)
                call assert_close("3d der dy2", dy_facade(2), dy_ref(2), REL_TOL)
                call assert_close("3d der dy3", dy_facade(3), dy_ref(3), REL_TOL)

                call evaluate_splines_3d_der2(spl, x_eval, y_facade, dy_facade, d2y_facade)
                call evaluate_splines_3d_der2_ref(spl, x_eval, y_ref, dy_ref, d2y_ref)
                call assert_close("3d der2 y", y_facade, y_ref, REL_TOL)
                call assert_close("3d der2 dy1", dy_facade(1), dy_ref(1), REL_TOL)
                call assert_close("3d der2 dy2", dy_facade(2), dy_ref(2), REL_TOL)
                call assert_close("3d der2 dy3", dy_facade(3), dy_ref(3), REL_TOL)
                call assert_close("3d der2 d2y1", d2y_facade(1), d2y_ref(1), REL_TOL)
                call assert_close("3d der2 d2y2", d2y_facade(2), d2y_ref(2), REL_TOL)
                call assert_close("3d der2 d2y3", d2y_facade(3), d2y_ref(3), REL_TOL)
                call assert_close("3d der2 d2y4", d2y_facade(4), d2y_ref(4), REL_TOL)
                call assert_close("3d der2 d2y5", d2y_facade(5), d2y_ref(5), REL_TOL)
                call assert_close("3d der2 d2y6", d2y_facade(6), d2y_ref(6), REL_TOL)
            end do

            call destroy_splines_3d(spl)
        end do

        print *, "  PASSED: 3D facade vs reference (value + derivatives)"
    end subroutine test_3d_evaluation

    subroutine assert_close(label, a, b, rel_tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: a, b, rel_tol

        real(dp) :: rel_err, scale_val

        scale_val = max(abs(a), abs(b))
        if (scale_val < 1.0d-100) return
        rel_err = abs(a - b) / scale_val
        if (rel_err > rel_tol) then
            print *, trim(label), " FAILED"
            print *, "  facade: ", a, " ref: ", b
            print *, "  rel_err: ", rel_err
            error stop "Facade vs reference oracle test failed"
        end if
    end subroutine assert_close

end program test_batch_eval_oracle
