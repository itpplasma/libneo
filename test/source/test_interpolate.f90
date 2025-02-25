program test_interpolate
    use libneo_kinds, only : dp
    use math_constants
    use libneo_util

    implicit none

    real(dp), parameter :: TOL_CRUDE = 1.0d-3, TOL = 1.0d-6, TOL_EXACT = 1.0d-11

    real(dp), parameter :: X_MIN = 1.23d0, X_MAX = TWOPI + 1.23d0

    call test_spline_1d(spline_order=3, periodic=.False.)
    call test_spline_1d(spline_order=3, periodic=.True.)
    call test_spline_1d(spline_order=5, periodic=.False.)
    call test_spline_1d(spline_order=5, periodic=.True.)


    call test_spline_2d(spline_order=[3,3], periodic=[.False., .False.])
    call test_spline_2d(spline_order=[5,5], periodic=[.False., .False.])
    call test_spline_2d(spline_order=[3,5], periodic=[.False., .False.])
    call test_spline_2d(spline_order=[5,3], periodic=[.False., .False.])
    call test_spline_2d(spline_order=[5,5], periodic=[.True., .False.])
    call test_spline_2d(spline_order=[5,5], periodic=[.True., .True.])


    call test_spline_3d(spline_order=[5,3,3], periodic=[.False.,.False.,.True.])
    call test_spline_3d(spline_order=[3,5,3], periodic=[.False., .True.,.True.])
    call test_spline_3d(spline_order=[3,3,3], periodic=[.True., .True., .True.])

    call bench_spline_1d

contains

    subroutine bench_spline_1d
        use interpolate
        integer, parameter :: n_iter = 1000
        integer, parameter :: N_POINTS = 100
        integer, parameter :: spline_order = 5
        logical, parameter :: periodic = .False.

        real(dp), dimension(N_POINTS) :: x, y

        real(dp) :: x_eval, expected, actual, t_start, t_end, t_elapsed

        integer :: iter

        type(SplineData1D) :: spl

        print *, "Benchmarking 1D spline with order ", spline_order, &
                 " and periodic = ", periodic

        call linspace(X_MIN, X_MAX, N_POINTS, x)

        y = cos(x)

        call construct_splines_1d(X_MIN, X_MAX, y, spline_order, periodic, spl)

        x_eval = 0.5d0*(x(10) + x(11))

        expected = cos(x_eval)

        call cpu_time(t_start)
        do iter = 1, n_iter
            call evaluate_splines_1d(spl, x_eval, actual)
        end do
        call cpu_time(t_end)
        t_elapsed = t_end - t_start
        print *, "expected, actual: ", expected, actual
        if (abs(expected - actual) > TOL) error stop

        print *, "Benchmark: ", n_iter, " evaluations took ", t_elapsed, " seconds"
        print *, "Average evaluation time: ", (t_elapsed / n_iter)*1.0d9, " nanoseconds per evaluation"

    end subroutine bench_spline_1d

    subroutine test_spline_1d(spline_order, periodic)
        use interpolate

        integer, parameter :: N_POINTS = 100

        integer, intent(in) :: spline_order
        logical, intent(in) :: periodic

        real(dp), dimension(N_POINTS) :: x, y

        real(dp) :: x_eval, expected, d_expected, d2_expected, &
            actual, d_actual, d2_actual

        type(SplineData1D) :: spl

        print *, "Testing 1D spline with order ", spline_order, &
                 " and periodic = ", periodic

        call linspace(X_MIN, X_MAX, N_POINTS, x)

        y = cos(x)

        call construct_splines_1d(X_MIN, X_MAX, y, spline_order, periodic, spl)

        x_eval = 0.5d0*(x(10) + x(11))

        expected = cos(x_eval)
        d_expected = -sin(x_eval)
        d2_expected = -cos(x_eval)

        call evaluate_splines_1d(spl, x_eval, actual)
        print *, "expected, actual: ", expected, actual
        if (abs(expected - actual) > TOL) error stop

        call evaluate_splines_1d_der(spl, x_eval, actual, d_actual)
        print *, "d_expected, d_actual: ", d_expected, d_actual
        if (abs(expected - actual) > TOL) error stop
        if (abs(d_expected - d_actual) > TOL) error stop

        call evaluate_splines_1d_der2(spl, x_eval, actual, d_actual, d2_actual)
        if (abs(expected - actual) > TOL) error stop
        if (abs(d_expected - d_actual) > TOL) error stop
        print *, "d2_expected, d2_actual: ", d2_expected, d2_actual
        if (abs(d2_expected - d2_actual) > 1d-3) error stop

        call destroy_splines_1d(spl)

    end subroutine test_spline_1d


    subroutine test_spline_2d(spline_order, periodic)
        use interpolate

        integer, intent(in) :: spline_order(2)
        logical, intent(in) :: periodic(2)


        integer, parameter :: N_POINTS(2) = [100, 103]

        real(dp), allocatable :: x1(:), x2(:), y(:,:)
        integer :: k1, k2

        real(dp) :: x_eval(2), expected, actual

        type(SplineData2D) :: spl

        allocate(x1(N_POINTS(1)), x2(N_POINTS(2)))
        allocate(y(N_POINTS(1), N_POINTS(2)))

        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)

        do k2 = 1, N_POINTS(2)
            do k1 = 1, N_POINTS(1)
                y(k1, k2) = cos(x1(k1))*cos(x2(k2))
            end do
        end do

        call construct_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], &
            y, spline_order, periodic, spl)

        x_eval(1) = (x1(55) + x1(56))/2.0d0
        x_eval(2) = (x2(28) + x2(29))/2.0d0

        expected = cos(x_eval(1))*cos(x_eval(2))

        call evaluate_splines_2d(spl, x_eval, actual)

        if (abs(expected - actual) > TOL) error stop

        call destroy_splines_2d(spl)

    end subroutine test_spline_2d


    subroutine test_spline_3d(spline_order, periodic)
            use interpolate

            integer, intent(in) :: spline_order(3)
            logical, intent(in) :: periodic(3)


            integer, parameter :: N_POINTS(3) = [82, 93, 87]

            real(dp), allocatable :: x1(:), x2(:), x3(:), y(:,:,:)
            integer :: k1, k2, k3

            real(dp) :: x_eval(3), expected, actual
            real(dp), dimension(3) :: d_expected, d_actual
            real(dp), dimension(6) :: d2_expected, d2_actual

            type(SplineData3D) :: spl

            allocate(x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3)))
            allocate(y(N_POINTS(1), N_POINTS(2), N_POINTS(3)))

            call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
            call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
            call linspace(X_MIN, X_MAX, N_POINTS(3), x3)

            do k3 = 1, N_POINTS(3)
                do k2 = 1, N_POINTS(2)
                    do k1 = 1, N_POINTS(1)
                       y(k1, k2, k3) = cos(x1(k1))*cos(x2(k2))*cos(x3(k3))
                  end do
                end do
            end do

            call construct_splines_3d(&
                [X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                y, spline_order, periodic, spl)

            x_eval(1) = x1(30)
            x_eval(2) = x2(28)
            x_eval(3) = x3(27)

            expected = cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))

            call evaluate_splines_3d(spl, x_eval, actual)

            if (abs(expected - actual) > TOL_EXACT) error stop

            x_eval(1) = (x1(30) + x1(31))/2.0d0
            x_eval(2) = (x2(28) + x2(29))/2.0d0
            x_eval(3) = (x3(27) + x3(28))/2.0d0

            expected = cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))
            d_expected(1) = -sin(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))
            d_expected(2) = -cos(x_eval(1))*sin(x_eval(2))*cos(x_eval(3))
            d_expected(3) = -cos(x_eval(1))*cos(x_eval(2))*sin(x_eval(3))
            d2_expected(1) = -cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3)) !11
            d2_expected(2) = sin(x_eval(1))*sin(x_eval(2))*cos(x_eval(3))  !12
            d2_expected(3) = sin(x_eval(1))*cos(x_eval(2))*sin(x_eval(3))  !13
            d2_expected(4) = -cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3)) !22
            d2_expected(5) = cos(x_eval(1))*sin(x_eval(2))*sin(x_eval(3))  !23
            d2_expected(6) = -cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3)) !33

            call evaluate_splines_3d(spl, x_eval, actual)

            if (abs(expected - actual) > TOL) error stop

            call evaluate_splines_3d_der2( &
                spl, x_eval, actual, d_actual, d2_actual)

            if (abs(expected - actual) > TOL) error stop
            if (any(abs(d_expected - d_actual) > TOL_CRUDE)) error stop
            if (any(abs(d2_expected - d2_actual) > TOL_CRUDE)) error stop

            call destroy_splines_3d(spl)

        end subroutine test_spline_3d

end program test_interpolate
