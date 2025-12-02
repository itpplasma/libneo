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
    call test_spline_2d_masked_ring()
    call test_spline_3d_masked_sphere()

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


    subroutine test_spline_2d_masked_ring()
        use interpolate

        integer, parameter :: N_POINTS(2) = [64, 65]
        real(dp), parameter :: X_MIN_LOC = -1.0d0, X_MAX_LOC = 1.0d0
        real(dp), parameter :: R_IN = 0.25d0, R_OUT = 0.85d0

        real(dp), allocatable :: x1(:), x2(:), y(:,:)
        logical, allocatable :: mask_nodes(:,:)
        real(dp) :: x_eval(2), r_eval, expected, actual, error_val
        logical :: is_valid
        integer :: i1, i2, unit_mask, unit_eval

        type(SplineData2D) :: spl

        allocate(x1(N_POINTS(1)), x2(N_POINTS(2)))
        allocate(y(N_POINTS(1), N_POINTS(2)))
        allocate(mask_nodes(N_POINTS(1), N_POINTS(2)))

        call linspace(X_MIN_LOC, X_MAX_LOC, N_POINTS(1), x1)
        call linspace(X_MIN_LOC, X_MAX_LOC, N_POINTS(2), x2)

        do i2 = 1, N_POINTS(2)
            do i1 = 1, N_POINTS(1)
                y(i1, i2) = x1(i1) + 2.0d0*x2(i2)
                if (sqrt(x1(i1)**2 + x2(i2)**2) >= R_IN .and. &
                    sqrt(x1(i1)**2 + x2(i2)**2) <= R_OUT) then
                    mask_nodes(i1, i2) = .true.
                else
                    mask_nodes(i1, i2) = .false.
                end if
            end do
        end do

        ! Write mask for visualization: one row per grid point
        open(newunit=unit_mask, file="interpolate_mask_2d_ring.dat", &
             status="replace", action="write")
        do i2 = 1, N_POINTS(2)
            do i1 = 1, N_POINTS(1)
                write(unit_mask, *) x1(i1), x2(i2), &
                    merge(1.0d0, 0.0d0, mask_nodes(i1, i2))
            end do
        end do
        close(unit_mask)

        call construct_splines_2d_masked( &
            [X_MIN_LOC, X_MIN_LOC], [X_MAX_LOC, X_MAX_LOC], y, &
            [5, 5], [.false., .false.], mask_nodes, spl)

        ! Point in the hole: should be marked invalid
        x_eval = [0.0d0, 0.0d0]
        call evaluate_splines_2d_masked(spl, x_eval, actual, is_valid)
        if (is_valid) error stop "Masked 2D eval inside hole reported valid"

        ! Point on the ring: should be valid and accurate
        x_eval = [0.5d0, 0.0d0]
        r_eval = sqrt(x_eval(1)**2 + x_eval(2)**2)
        if (r_eval <= R_IN .or. r_eval >= R_OUT) &
            error stop "2D masked test point not in ring"

        expected = x_eval(1) + 2.0d0*x_eval(2)
        call evaluate_splines_2d_masked(spl, x_eval, actual, is_valid)
        if (.not. is_valid) error stop "Masked 2D eval on ring reported invalid"
        if (abs(expected - actual) > TOL) &
            error stop "Masked 2D eval error too large"

        ! Write evaluation results and errors at cell centers
        open(newunit=unit_eval, file="interpolate_mask_2d_ring_eval.dat", &
             status="replace", action="write")
        do i2 = 1, N_POINTS(2) - 1
            do i1 = 1, N_POINTS(1) - 1
                x_eval(1) = 0.5d0*(x1(i1) + x1(i1 + 1))
                x_eval(2) = 0.5d0*(x2(i2) + x2(i2 + 1))
                expected = x_eval(1) + 2.0d0*x_eval(2)
                call evaluate_splines_2d_masked(spl, x_eval, actual, is_valid)
                if (is_valid) then
                    error_val = actual - expected
                else
                    error_val = 0.0d0
                end if
                write(unit_eval, *) x_eval(1), x_eval(2), expected, actual, &
                    error_val, merge(1.0d0, 0.0d0, is_valid)
            end do
        end do
        close(unit_eval)

        call destroy_splines_2d(spl)

    end subroutine test_spline_2d_masked_ring


    subroutine test_spline_3d_masked_sphere()
        use interpolate

        integer, parameter :: N_POINTS(3) = [32, 33, 34]
        real(dp), parameter :: X_MIN_LOC = -1.0d0, X_MAX_LOC = 1.0d0
        real(dp), parameter :: R_IN = 0.25d0, R_OUT = 0.9d0

        real(dp), allocatable :: x1(:), x2(:), x3(:), y(:,:,:)
        logical, allocatable :: mask_nodes(:,:,:)
        real(dp) :: x_eval(3), r_eval, expected, actual, error_val
        logical :: is_valid
        integer :: i1, i2, i3, unit_mask, unit_eval, mid3

        type(SplineData3D) :: spl

        allocate(x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3)))
        allocate(y(N_POINTS(1), N_POINTS(2), N_POINTS(3)))
        allocate(mask_nodes(N_POINTS(1), N_POINTS(2), N_POINTS(3)))

        call linspace(X_MIN_LOC, X_MAX_LOC, N_POINTS(1), x1)
        call linspace(X_MIN_LOC, X_MAX_LOC, N_POINTS(2), x2)
        call linspace(X_MIN_LOC, X_MAX_LOC, N_POINTS(3), x3)

        do i3 = 1, N_POINTS(3)
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y(i1, i2, i3) = x1(i1) + 2.0d0*x2(i2) + 3.0d0*x3(i3)
                    if (sqrt(x1(i1)**2 + x2(i2)**2 + x3(i3)**2) >= R_IN .and. &
                        sqrt(x1(i1)**2 + x2(i2)**2 + x3(i3)**2) <= R_OUT) then
                        mask_nodes(i1, i2, i3) = .true.
                    else
                        mask_nodes(i1, i2, i3) = .false.
                    end if
                end do
            end do
        end do

        mid3 = 1
        do i3 = 1, N_POINTS(3)
            if (abs(x3(i3)) < abs(x3(mid3))) mid3 = i3
        end do

        ! Write mask slice near x3 = 0 for visualization
        open(newunit=unit_mask, file="interpolate_mask_3d_shell_slice.dat", &
             status="replace", action="write")
        do i2 = 1, N_POINTS(2)
            do i1 = 1, N_POINTS(1)
                write(unit_mask, *) x1(i1), x2(i2), &
                    merge(1.0d0, 0.0d0, mask_nodes(i1, i2, mid3))
            end do
        end do
        close(unit_mask)

        call construct_splines_3d_masked( &
            [X_MIN_LOC, X_MIN_LOC, X_MIN_LOC], &
            [X_MAX_LOC, X_MAX_LOC, X_MAX_LOC], y, &
            [3, 3, 3], [.false., .false., .false.], mask_nodes, spl)

        ! Write evaluation results and errors on slice near x3 = 0
        open(newunit=unit_eval, &
             file="interpolate_mask_3d_shell_slice_eval.dat", &
             status="replace", action="write")
        do i2 = 1, N_POINTS(2) - 1
            do i1 = 1, N_POINTS(1) - 1
                x_eval(1) = 0.5d0*(x1(i1) + x1(i1 + 1))
                x_eval(2) = 0.5d0*(x2(i2) + x2(i2 + 1))
                x_eval(3) = x3(mid3)
                expected = x_eval(1) + 2.0d0*x_eval(2) + 3.0d0*x_eval(3)
                call evaluate_splines_3d_masked(spl, x_eval, actual, is_valid)
                if (is_valid) then
                    error_val = actual - expected
                else
                    error_val = 0.0d0
                end if
                write(unit_eval, *) x_eval(1), x_eval(2), expected, actual, &
                    error_val, merge(1.0d0, 0.0d0, is_valid)
            end do
        end do
        close(unit_eval)

        ! Point in the hole: should be marked invalid
        x_eval = [0.0d0, 0.0d0, 0.0d0]
        call evaluate_splines_3d_masked(spl, x_eval, actual, is_valid)
        if (is_valid) error stop "Masked 3D eval inside sphere reported valid"

        ! Point in the spherical shell: should be valid and accurate
        x_eval = [0.4d0, 0.0d0, 0.0d0]
        r_eval = sqrt(sum(x_eval**2))
        if (r_eval <= R_IN .or. r_eval >= R_OUT) &
            error stop "3D masked test point not in shell"

        expected = x_eval(1) + 2.0d0*x_eval(2) + 3.0d0*x_eval(3)
        call evaluate_splines_3d_masked(spl, x_eval, actual, is_valid)
        if (.not. is_valid) error stop "Masked 3D eval on shell reported invalid"
        if (abs(expected - actual) > TOL) &
            error stop "Masked 3D eval error too large"

        call destroy_splines_3d(spl)

    end subroutine test_spline_3d_masked_sphere

end program test_interpolate
