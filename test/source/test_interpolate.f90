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

    call test_spline_1d_lsq_scattered()
    call test_spline_2d_lsq_full_grid()
    call test_spline_2d_lsq_circular_domain()
    call test_spline_3d_lsq_spherical_and_toroidal()

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

    subroutine test_spline_1d_lsq_scattered()
        use interpolate

        integer, parameter :: N_POINTS_GRID = 100
        integer, parameter :: N_POINTS_DATA = N_POINTS_GRID
        integer, parameter :: ORDER = 5
        logical, parameter :: PERIODIC = .false.

        real(dp) :: x_grid(N_POINTS_GRID), y_grid(N_POINTS_GRID)
        real(dp) :: x_data(N_POINTS_DATA), f_data(N_POINTS_DATA)
        real(dp) :: x_eval, f_true, f_grid, f_lsq

        type(SplineData1D) :: spl_grid
        type(SplineData1D) :: spl_lsq

        integer :: i

        print *, "Testing 1D LSQ spline with scattered data"

        call linspace(X_MIN, X_MAX, N_POINTS_GRID, x_grid)

        do i = 1, N_POINTS_GRID
            y_grid(i) = sin(2.0d0*x_grid(i)) + 0.5d0*cos(3.0d0*x_grid(i))
        end do

        call construct_splines_1d(X_MIN, X_MAX, y_grid, ORDER, PERIODIC, spl_grid)

        do i = 1, N_POINTS_DATA
            x_data(i) = x_grid(i)
            f_data(i) = y_grid(i)
        end do

        call construct_splines_1d_lsq(X_MIN, X_MAX, ORDER, PERIODIC, &
            N_POINTS_GRID, x_data, f_data, spl_lsq)

        do i = 1, 50
            x_eval = X_MIN + (X_MAX - X_MIN)*dble(i)/51.0d0
            f_true = sin(2.0d0*x_eval) + 0.5d0*cos(3.0d0*x_eval)

            call evaluate_splines_1d(spl_grid, x_eval, f_grid)
            call evaluate_splines_1d(spl_lsq, x_eval, f_lsq)

            if (abs(f_grid - f_true) > TOL) then
                error stop "1D grid spline does not match analytic function"
            end if
            if (abs(f_lsq - f_true) > TOL) then
                error stop "1D LSQ spline does not match analytic function"
            end if
        end do

        call destroy_splines_1d(spl_grid)
        call destroy_splines_1d(spl_lsq)

    end subroutine test_spline_1d_lsq_scattered

    subroutine test_spline_2d_lsq_full_grid()
        use interpolate

        integer, parameter :: N_POINTS1 = 40
        integer, parameter :: N_POINTS2 = 36
        integer, parameter :: ORDER1 = 5
        integer, parameter :: ORDER2 = 5

        real(dp) :: x_min(2), x_max(2)
        real(dp) :: x1(N_POINTS1), x2(N_POINTS2)
        real(dp) :: y_grid(N_POINTS1, N_POINTS2)
        real(dp) :: x_data(N_POINTS1*N_POINTS2)
        real(dp) :: y_data(N_POINTS1*N_POINTS2)
        real(dp) :: f_data(N_POINTS1*N_POINTS2)

        type(SplineData2D) :: spl_grid
        type(SplineData2D) :: spl_lsq

        real(dp) :: x_eval(2), f_true, f_grid, f_lsq
        real(dp) :: max_err_grid, max_err_lsq
        integer :: i1, i2, idx
        real(dp) :: xmin_loc, xmax_loc

        print *, "Testing 2D LSQ spline on full grid (sanity check)"

        xmin_loc = 1.23d0
        xmax_loc = TWOPI + 1.23d0

        x_min(1) = xmin_loc
        x_min(2) = xmin_loc
        x_max(1) = xmax_loc
        x_max(2) = xmax_loc

        call linspace(xmin_loc, xmax_loc, N_POINTS1, x1)
        call linspace(xmin_loc, xmax_loc, N_POINTS2, x2)

        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                y_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        call construct_splines_2d(x_min, x_max, y_grid, &
            [ORDER1, ORDER2], [.false., .false.], spl_grid)

        idx = 0
        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                idx = idx + 1
                x_data(idx) = x1(i1)
                y_data(idx) = x2(i2)
                f_data(idx) = y_grid(i1, i2)
            end do
        end do

        call construct_splines_2d_lsq(x_min, x_max, &
            [ORDER1, ORDER2], [.false., .false.], &
            [N_POINTS1, N_POINTS2], x_data, y_data, f_data, spl_lsq)

        max_err_grid = 0.0d0
        max_err_lsq = 0.0d0

        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                x_eval(1) = x1(i1)
                x_eval(2) = x2(i2)

                f_true = cos(x_eval(1))*cos(2.0d0*x_eval(2))

                call evaluate_splines_2d(spl_grid, x_eval, f_grid)
                call evaluate_splines_2d(spl_lsq, x_eval, f_lsq)

                max_err_grid = max(max_err_grid, abs(f_grid - f_true))
                max_err_lsq = max(max_err_lsq, abs(f_lsq - f_true))
            end do
        end do

        print *, "  2D full-grid LSQ: max grid error =", max_err_grid
        print *, "  2D full-grid LSQ: max LSQ  error =", max_err_lsq

        if (max_err_grid > TOL) then
            error stop "2D grid spline does not match analytic (full-grid LSQ)"
        end if
        if (max_err_lsq > TOL) then
            error stop "2D LSQ spline does not match analytic (full-grid LSQ)"
        end if

        call destroy_splines_2d(spl_grid)
        call destroy_splines_2d(spl_lsq)

    end subroutine test_spline_2d_lsq_full_grid

    subroutine test_spline_2d_lsq_circular_domain()
        use interpolate

        integer, parameter :: N_POINTS1 = 40
        integer, parameter :: N_POINTS2 = 36
        integer, parameter :: ORDER1 = 5
        integer, parameter :: ORDER2 = 5

        integer, parameter :: MAX_DATA = N_POINTS1*N_POINTS2

        real(dp) :: x_min(2), x_max(2)
        real(dp) :: x1(N_POINTS1), x2(N_POINTS2)
        real(dp) :: y_grid(N_POINTS1, N_POINTS2)

        real(dp) :: x_data(MAX_DATA), y_data(MAX_DATA), f_data(MAX_DATA)

        type(SplineData2D) :: spl_grid
        type(SplineData2D) :: spl_lsq

        real(dp) :: x_eval(2), f_true, f_grid, f_lsq
        real(dp) :: max_err_grid, max_err_lsq
        real(dp) :: xc, yc, r0, dx, dy
        real(dp) :: theta, r
        real(dp) :: xmin_loc, xmax_loc
        real(dp), parameter :: TOL_LSQ_MASK = 5.0d-2

        integer :: i1, i2, ndata, i

        print *, "Testing 2D LSQ spline with circular domain cut-out"

        max_err_grid = 0.0d0
        max_err_lsq = 0.0d0

        xmin_loc = 1.23d0
        xmax_loc = TWOPI + 1.23d0

        x_min(1) = xmin_loc
        x_min(2) = xmin_loc
        x_max(1) = xmax_loc
        x_max(2) = xmax_loc

        call linspace(xmin_loc, xmax_loc, N_POINTS1, x1)
        call linspace(xmin_loc, xmax_loc, N_POINTS2, x2)

        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
                y_grid(i1, i2) = cos(x1(i1))*cos(2.0d0*x2(i2))
            end do
        end do

        call construct_splines_2d(x_min, x_max, y_grid, &
            [ORDER1, ORDER2], [.false., .false.], spl_grid)

        xc = 0.5d0*(xmin_loc + xmax_loc)
        yc = 0.5d0*(xmin_loc + xmax_loc)
        r0 = 0.45d0*(xmax_loc - xmin_loc)

        ndata = 0
        do i2 = 1, N_POINTS2
            do i1 = 1, N_POINTS1
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

        call construct_splines_2d_lsq(x_min, x_max, &
            [ORDER1, ORDER2], [.false., .false.], &
            [N_POINTS1, N_POINTS2], x_data(1:ndata), y_data(1:ndata), &
            f_data(1:ndata), spl_lsq)

        do i = 1, 40
            theta = TWOPI*dble(i)/40.0d0
            r = 0.4d0*(xmax_loc - xmin_loc)
            x_eval(1) = xc + r*cos(theta)
            x_eval(2) = yc + r*sin(theta)

            f_true = cos(x_eval(1))*cos(2.0d0*x_eval(2))

            call evaluate_splines_2d(spl_grid, x_eval, f_grid)
            call evaluate_splines_2d(spl_lsq, x_eval, f_lsq)

            max_err_grid = max(max_err_grid, abs(f_grid - f_true))
            max_err_lsq = max(max_err_lsq, abs(f_lsq - f_true))
        end do

        print *, "  2D LSQ: max grid error  =", max_err_grid
        print *, "  2D LSQ: max LSQ  error  =", max_err_lsq

        if (max_err_grid > TOL) then
            error stop "2D grid spline does not match analytic function"
        end if
        if (max_err_lsq > TOL_LSQ_MASK) then
            error stop "2D LSQ spline does not match analytic function"
        end if

        call destroy_splines_2d(spl_grid)
        call destroy_splines_2d(spl_lsq)

    end subroutine test_spline_2d_lsq_circular_domain

    subroutine test_spline_3d_lsq_spherical_and_toroidal()
        use interpolate

        integer, parameter :: N_POINTS1 = 18
        integer, parameter :: N_POINTS2 = 16
        integer, parameter :: N_POINTS3 = 14

        integer, parameter :: ORDER1 = 3
        integer, parameter :: ORDER2 = 3
        integer, parameter :: ORDER3 = 3

        integer, parameter :: MAX_DATA = N_POINTS1*N_POINTS2*N_POINTS3

        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x1(N_POINTS1), x2(N_POINTS2), x3(N_POINTS3)
        real(dp) :: y_grid(N_POINTS1, N_POINTS2, N_POINTS3)

        real(dp) :: xs_data(MAX_DATA), ys_data(MAX_DATA), zs_data(MAX_DATA)
        real(dp) :: fs_sphere(MAX_DATA), fs_torus(MAX_DATA)

        type(SplineData3D) :: spl_grid
        type(SplineData3D) :: spl_lsq_sphere
        type(SplineData3D) :: spl_lsq_torus

        real(dp) :: x_eval(3), f_true, f_grid, f_sphere, f_torus

        real(dp) :: xc, yc, zc, r_sphere, r_major, r_minor
        real(dp) :: dx, dy, dz, R_xy, dist_torus
        real(dp) :: theta1, theta2, r, theta, phi
        real(dp) :: xmin_loc, xmax_loc
        real(dp) :: max_err_grid_sph, max_err_lsq_sph
        real(dp) :: max_err_grid_tor, max_err_lsq_tor
        real(dp), parameter :: TOL_GRID_3D_LSQ = 1.0d-3
        real(dp), parameter :: TOL_LSQ_3D_MASK = 3.0d-1

        integer :: i1, i2, i3, n_sphere, n_torus, i

        print *, "Testing 3D LSQ spline with spherical and toroidal domains"

        xmin_loc = 1.23d0
        xmax_loc = TWOPI + 1.23d0

        x_min = xmin_loc
        x_max = xmax_loc

        do i1 = 1, N_POINTS1
            x1(i1) = xmin_loc + (xmax_loc - xmin_loc)*dble(i1 - 1)/ &
                dble(N_POINTS1 - 1)
        end do
        do i2 = 1, N_POINTS2
            x2(i2) = xmin_loc + (xmax_loc - xmin_loc)*dble(i2 - 1)/ &
                dble(N_POINTS2 - 1)
        end do
        do i3 = 1, N_POINTS3
            x3(i3) = xmin_loc + (xmax_loc - xmin_loc)*dble(i3 - 1)/ &
                dble(N_POINTS3 - 1)
        end do

        do i3 = 1, N_POINTS3
            do i2 = 1, N_POINTS2
                do i1 = 1, N_POINTS1
                    y_grid(i1, i2, i3) = cos(x1(i1))*cos(x2(i2))*cos(x3(i3))
                end do
            end do
        end do

        call construct_splines_3d(x_min, x_max, y_grid, &
            [ORDER1, ORDER2, ORDER3], &
            [.false., .false., .false.], spl_grid)

        xc = 0.5d0*(xmin_loc + xmax_loc)
        yc = 0.5d0*(xmin_loc + xmax_loc)
        zc = 0.5d0*(xmin_loc + xmax_loc)
        r_sphere = 0.35d0*(xmax_loc - xmin_loc)

        r_major = 0.35d0*(xmax_loc - xmin_loc)
        r_minor = 0.15d0*(xmax_loc - xmin_loc)

        n_sphere = 0
        n_torus = 0
        max_err_grid_sph = 0.0d0
        max_err_lsq_sph = 0.0d0
        max_err_grid_tor = 0.0d0
        max_err_lsq_tor = 0.0d0

        do i3 = 1, N_POINTS3
            do i2 = 1, N_POINTS2
                do i1 = 1, N_POINTS1
                    dx = x1(i1) - xc
                    dy = x2(i2) - yc
                    dz = x3(i3) - zc

                    if (dx*dx + dy*dy + dz*dz <= r_sphere*r_sphere) then
                        n_sphere = n_sphere + 1
                        xs_data(n_sphere) = x1(i1)
                        ys_data(n_sphere) = x2(i2)
                        zs_data(n_sphere) = x3(i3)
                        fs_sphere(n_sphere) = cos(x1(i1))*cos(x2(i2))*cos(x3(i3))
                    end if

                    R_xy = sqrt(dx*dx + dy*dy)
                    dist_torus = sqrt((R_xy - r_major)**2 + dz*dz)
                    if (dist_torus <= r_minor) then
                        n_torus = n_torus + 1
                        xs_data(MAX_DATA - n_torus + 1) = x1(i1)
                        ys_data(MAX_DATA - n_torus + 1) = x2(i2)
                        zs_data(MAX_DATA - n_torus + 1) = x3(i3)
                        fs_torus(MAX_DATA - n_torus + 1) = &
                            cos(x1(i1))*cos(x2(i2))*cos(x3(i3))
                    end if
                end do
            end do
        end do

        if (n_sphere < 10) then
            error stop "Too few points in spherical LSQ sample"
        end if
        if (n_torus < 10) then
            error stop "Too few points in toroidal LSQ sample"
        end if

        call construct_splines_3d_lsq(x_min, x_max, &
            [ORDER1, ORDER2, ORDER3], &
            [.false., .false., .false.], &
            [N_POINTS1, N_POINTS2, N_POINTS3], &
            xs_data(1:n_sphere), ys_data(1:n_sphere), &
            zs_data(1:n_sphere), fs_sphere(1:n_sphere), spl_lsq_sphere)

        call construct_splines_3d_lsq(x_min, x_max, &
            [ORDER1, ORDER2, ORDER3], &
            [.false., .false., .false.], &
            [N_POINTS1, N_POINTS2, N_POINTS3], &
            xs_data(MAX_DATA - n_torus + 1:MAX_DATA), &
            ys_data(MAX_DATA - n_torus + 1:MAX_DATA), &
            zs_data(MAX_DATA - n_torus + 1:MAX_DATA), &
            fs_torus(MAX_DATA - n_torus + 1:MAX_DATA), spl_lsq_torus)

        do i = 1, 40
            theta1 = TWOPI*dble(i)/40.0d0
            theta2 = 0.5d0*TWOPI*dble(i)/40.0d0
            r = 0.3d0*(xmax_loc - xmin_loc)

            x_eval(1) = xc + r*cos(theta1)*cos(theta2)
            x_eval(2) = yc + r*sin(theta1)*cos(theta2)
            x_eval(3) = zc + r*sin(theta2)

            f_true = cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))

            call evaluate_splines_3d(spl_grid, x_eval, f_grid)
            call evaluate_splines_3d(spl_lsq_sphere, x_eval, f_sphere)

            max_err_grid_sph = max(max_err_grid_sph, abs(f_grid - f_true))
            max_err_lsq_sph = max(max_err_lsq_sph, abs(f_sphere - f_true))
        end do

        do i = 1, 40
            theta = TWOPI*dble(i)/40.0d0
            phi = 0.5d0*TWOPI*dble(i)/40.0d0

            x_eval(1) = xc + (r_major + r_minor*cos(phi))*cos(theta)
            x_eval(2) = yc + (r_major + r_minor*cos(phi))*sin(theta)
            x_eval(3) = zc + r_minor*sin(phi)

            f_true = cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))

            call evaluate_splines_3d(spl_grid, x_eval, f_grid)
            call evaluate_splines_3d(spl_lsq_torus, x_eval, f_torus)

            max_err_grid_tor = max(max_err_grid_tor, abs(f_grid - f_true))
            max_err_lsq_tor = max(max_err_lsq_tor, abs(f_torus - f_true))
        end do

        print *, "  3D LSQ sphere: max grid error =", max_err_grid_sph
        print *, "  3D LSQ sphere: max LSQ  error =", max_err_lsq_sph
        print *, "  3D LSQ torus : max grid error =", max_err_grid_tor
        print *, "  3D LSQ torus : max LSQ  error =", max_err_lsq_tor

        if (max_err_grid_sph > TOL_GRID_3D_LSQ) then
            error stop "3D grid spline does not match analytic function (sphere)"
        end if
        if (max_err_lsq_sph > TOL_LSQ_3D_MASK) then
            error stop "3D spherical LSQ spline does not match analytic"
        end if
        if (max_err_grid_tor > TOL_GRID_3D_LSQ) then
            error stop "3D grid spline does not match analytic (torus)"
        end if
        if (max_err_lsq_tor > TOL_LSQ_3D_MASK) then
            error stop "3D toroidal LSQ spline does not match analytic"
        end if

        call destroy_splines_3d(spl_grid)
        call destroy_splines_3d(spl_lsq_sphere)
        call destroy_splines_3d(spl_lsq_torus)

    end subroutine test_spline_3d_lsq_spherical_and_toroidal

end program test_interpolate
