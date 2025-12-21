program test_batch_interpolate_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_util, only: linspace
    use batch_interpolate_types, only: BatchSplineData1D, BatchSplineData2D, BatchSplineData3D
    use batch_interpolate_1d, only: construct_batch_splines_1d, &
                                    construct_batch_splines_1d_legacy, &
                                    destroy_batch_splines_1d, &
                                    evaluate_batch_splines_1d_many
    use batch_interpolate_2d, only: construct_batch_splines_2d, &
                                    construct_batch_splines_2d_legacy, &
                                    destroy_batch_splines_2d, &
                                    evaluate_batch_splines_2d_many
    use batch_interpolate_3d, only: construct_batch_splines_3d, &
                                    construct_batch_splines_3d_legacy, &
                                    destroy_batch_splines_3d, &
                                    evaluate_batch_splines_3d_many

    implicit none

    call test_oracle_1d()
    call test_oracle_2d()
    call test_oracle_3d()

    print *, "All batch spline oracle tests passed!"

contains

    subroutine test_oracle_1d()
        real(dp), parameter :: tol = 5.0d-12
        real(dp), parameter :: x_min = 1.23d0
        real(dp), parameter :: x_max = 7.89d0
        integer, parameter :: n_points = 64
        integer, parameter :: nq = 5
        integer, parameter :: n_eval = 200

        integer :: order, iq, ip, ie
        integer :: iper
        logical :: periodic
        real(dp) :: x_grid(n_points)
        real(dp) :: y_grid(n_points, nq)
        real(dp) :: x_eval(n_eval)
        real(dp) :: y_new(nq, n_eval)
        real(dp) :: y_old(nq, n_eval)
        real(dp) :: period, r

        type(BatchSplineData1D) :: spl_new, spl_old

        call linspace(x_min, x_max, n_points, x_grid)

        do iq = 1, nq
            do ip = 1, n_points
                y_grid(ip, iq) = cos(x_grid(ip) + 0.5d0*dble(iq - 1)) + &
                                 0.1d0*sin(2.0d0*x_grid(ip))
            end do
        end do

        do order = 3, 5
            do iper = 0, 1
                periodic = (iper == 1)
                call construct_batch_splines_1d(x_min, x_max, y_grid, order, periodic, &
                                                spl_new)
                call construct_batch_splines_1d_legacy(x_min, x_max, y_grid, order, periodic, &
                                                       spl_old)

                period = (x_max - x_min)
                do ie = 1, n_eval
                    call random_number(r)
                    if (periodic) then
                        x_eval(ie) = x_min + (4.0d0*r - 2.0d0)*period
                    else
                        x_eval(ie) = x_min + r*period
                    end if
                end do

                call evaluate_batch_splines_1d_many(spl_new, x_eval, y_new)
                call evaluate_batch_splines_1d_many(spl_old, x_eval, y_old)
                call assert_close_2d("1d many oracle", y_new, y_old, tol)

                call destroy_batch_splines_1d(spl_new)
                call destroy_batch_splines_1d(spl_old)
            end do
        end do

        print *, "  PASSED: 1D new vs legacy (many-point, periodic/nonperiodic, order 3..5)"
    end subroutine test_oracle_1d


    subroutine test_oracle_2d()
        real(dp), parameter :: tol = 1.0d-11
        real(dp), parameter :: x_min(2) = [1.23d0, -0.75d0]
        real(dp), parameter :: x_max(2) = [4.56d0,  2.25d0]
        integer, parameter :: n1 = 32, n2 = 33
        integer, parameter :: nq = 4
        integer, parameter :: n_eval = 250
        integer, parameter :: n_cases = 3

        integer :: iq, i1, i2, ie
        integer :: icase
        integer :: order(2)
        logical :: periodic(2)
        real(dp) :: x1(n1), x2(n2)
        real(dp) :: y_grid(n1, n2, nq)
        real(dp) :: x_eval(2, n_eval)
        real(dp) :: y_new(nq, n_eval)
        real(dp) :: y_old(nq, n_eval)
        real(dp) :: period(2), r(2)

        type(BatchSplineData2D) :: spl_new, spl_old
        integer :: order_cases(2, n_cases)
        logical :: periodic_cases(2, n_cases)

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)

        do iq = 1, nq
            do i2 = 1, n2
                do i1 = 1, n1
                    y_grid(i1, i2, iq) = cos(x1(i1) + 0.25d0*dble(iq - 1)) * &
                                         cos(x2(i2) - 0.10d0*dble(iq - 1))
                end do
            end do
        end do

        period = x_max - x_min

        order_cases = reshape([3, 3, 5, 5, 5, 3], [2, n_cases])
        periodic_cases = reshape([.false., .false., .true., .true., .true., .false.], &
                                 [2, n_cases])

        do icase = 1, n_cases
            order = order_cases(:, icase)
            periodic = periodic_cases(:, icase)

            call construct_batch_splines_2d(x_min, x_max, y_grid, order, periodic, spl_new)
            call construct_batch_splines_2d_legacy(x_min, x_max, y_grid, order, periodic, &
                                                   spl_old)

            do ie = 1, n_eval
                call random_number(r)
                if (periodic(1)) then
                    x_eval(1, ie) = x_min(1) + (4.0d0*r(1) - 2.0d0)*period(1)
                else
                    x_eval(1, ie) = x_min(1) + r(1)*period(1)
                end if
                if (periodic(2)) then
                    x_eval(2, ie) = x_min(2) + (4.0d0*r(2) - 2.0d0)*period(2)
                else
                    x_eval(2, ie) = x_min(2) + r(2)*period(2)
                end if
            end do

            call evaluate_batch_splines_2d_many(spl_new, x_eval, y_new)
            call evaluate_batch_splines_2d_many(spl_old, x_eval, y_old)
            call assert_close_2d("2d many oracle", y_new, y_old, tol)

            call destroy_batch_splines_2d(spl_new)
            call destroy_batch_splines_2d(spl_old)
        end do

        print *, "  PASSED: 2D new vs legacy (many-point, mixed periodic, order sets)"
    end subroutine test_oracle_2d


    subroutine test_oracle_3d()
        real(dp), parameter :: tol = 2.0d-11
        real(dp), parameter :: x_min(3) = [0.10d0, -1.20d0, 2.20d0]
        real(dp), parameter :: x_max(3) = [1.70d0,  0.80d0, 3.70d0]
        integer, parameter :: n1 = 16, n2 = 17, n3 = 18
        integer, parameter :: nq = 3
        integer, parameter :: n_eval = 200
        integer, parameter :: n_cases = 2

        integer :: iq, i1, i2, i3, ie
        integer :: icase
        integer :: order(3)
        logical :: periodic(3)
        real(dp) :: x1(n1), x2(n2), x3(n3)
        real(dp) :: y_grid(n1, n2, n3, nq)
        real(dp) :: x_eval(3, n_eval)
        real(dp) :: y_new(nq, n_eval)
        real(dp) :: y_old(nq, n_eval)
        real(dp) :: period(3), r(3)

        type(BatchSplineData3D) :: spl_new, spl_old
        integer :: order_cases(3, n_cases)
        logical :: periodic_cases(3, n_cases)

        call linspace(x_min(1), x_max(1), n1, x1)
        call linspace(x_min(2), x_max(2), n2, x2)
        call linspace(x_min(3), x_max(3), n3, x3)

        do iq = 1, nq
            do i3 = 1, n3
                do i2 = 1, n2
                    do i1 = 1, n1
                        y_grid(i1, i2, i3, iq) = &
                            cos(x1(i1) + 0.10d0*dble(iq - 1)) * &
                            cos(x2(i2) - 0.15d0*dble(iq - 1)) * &
                            cos(x3(i3) + 0.05d0*dble(iq - 1))
                    end do
                end do
            end do
        end do

        period = x_max - x_min

        order_cases = reshape([5, 3, 3, 5, 5, 3], [3, n_cases])
        periodic_cases = reshape([.true., .true., .true., .true., .false., .true.], &
                                 [3, n_cases])

        do icase = 1, n_cases
            order = order_cases(:, icase)
            periodic = periodic_cases(:, icase)

            call construct_batch_splines_3d(x_min, x_max, y_grid, order, periodic, spl_new)
            call construct_batch_splines_3d_legacy(x_min, x_max, y_grid, order, periodic, &
                                                   spl_old)

            do ie = 1, n_eval
                call random_number(r)
                x_eval(1, ie) = x_min(1) + (4.0d0*r(1) - 2.0d0)*period(1)
                if (periodic(2)) then
                    x_eval(2, ie) = x_min(2) + (4.0d0*r(2) - 2.0d0)*period(2)
                else
                    x_eval(2, ie) = x_min(2) + r(2)*period(2)
                end if
                x_eval(3, ie) = x_min(3) + (4.0d0*r(3) - 2.0d0)*period(3)
            end do

            call evaluate_batch_splines_3d_many(spl_new, x_eval, y_new)
            call evaluate_batch_splines_3d_many(spl_old, x_eval, y_old)
            call assert_close_2d("3d many oracle", y_new, y_old, tol)

            call destroy_batch_splines_3d(spl_new)
            call destroy_batch_splines_3d(spl_old)
        end do

        print *, "  PASSED: 3D new vs legacy (many-point, mixed periodic, order sets)"
    end subroutine test_oracle_3d


    subroutine assert_close_2d(label, a, b, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(in) :: tol

        real(dp) :: err_max

        if (any(shape(a) /= shape(b))) then
            error stop label // ": shape mismatch"
        end if

        err_max = maxval(abs(a - b))
        if (err_max > tol) then
            print *, trim(label), ": max abs error = ", err_max, " tol = ", tol
            error stop trim(label) // ": values differ"
        end if
    end subroutine assert_close_2d

end program test_batch_interpolate_oracle
