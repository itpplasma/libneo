program test_interpolate
    use math_constants
    use util

    implicit none

    real(8), parameter :: TOL = 1.0d-6

    real(8), parameter :: X_MIN = 0.0d0, X_MAX = TWOPI

    call test_spline_1d(spline_order=3, periodic=.False.)
    call test_spline_1d(spline_order=3, periodic=.True.)
    call test_spline_1d(spline_order=4, periodic=.False.)
    call test_spline_1d(spline_order=4, periodic=.True.)
    call test_spline_1d(spline_order=5, periodic=.False.)
    call test_spline_1d(spline_order=5, periodic=.True.)


    call test_spline_2d(spline_order=(/3,3/), periodic=(/.False., .False./))
    call test_spline_2d(spline_order=(/5,5/), periodic=(/.False., .False./))
    call test_spline_2d(spline_order=(/3,5/), periodic=(/.False., .False./))
    call test_spline_2d(spline_order=(/5,3/), periodic=(/.False., .False./))
    call test_spline_2d(spline_order=(/5,5/), periodic=(/.True., .False./))
    call test_spline_2d(spline_order=(/5,5/), periodic=(/.True., .True./))


    call test_spline_3d(spline_order=(/3,3,3/), periodic=(/.False., .False., .True./))

contains

    subroutine test_spline_1d(spline_order, periodic)
        use interpolate

        integer, parameter :: N_POINTS = 100

        integer, intent(in) :: spline_order
        logical, intent(in) :: periodic

        real(8), dimension(N_POINTS) :: x, y

        real(8) :: x_eval, expected, actual

        type(SplineData1D) :: spl

        call linspace(0.0d0, 2.0d0 * pi, N_POINTS, x)

        y = cos(x)

        call construct_splines_1d(x, y, spline_order, periodic, spl)

        x_eval = (x(30) + x(31))/2.0d0

        expected = cos(x_eval)

        call evaluate_splines_1d(x_eval, spl, actual)

        if (abs(expected - actual) > TOL) error stop

        call destroy_splines_1d(spl)

    end subroutine test_spline_1d


    subroutine test_spline_2d(spline_order, periodic)
        use interpolate

        integer, intent(in) :: spline_order(2)
        logical, intent(in) :: periodic(2)


        integer, parameter :: N_POINTS(2) = (/100, 103/)

        real(8), allocatable :: x1(:), x2(:), y(:,:)
        integer :: k1, k2

        real(8) :: x_eval(2), expected, actual

        type(SplineData2D) :: spl

        allocate(x1(N_POINTS(1)), x2(N_POINTS(2)))
        allocate(y(N_POINTS(1), N_POINTS(2)))

        call linspace(0.0d0, 2.0d0 * pi, N_POINTS(1), x1)
        call linspace(0.0d0, 2.0d0 * pi, N_POINTS(2), x2)

        do k2 = 1, N_POINTS(2)
            do k1 = 1, N_POINTS(1)
                y(k1, k2) = cos(x1(k1))*cos(x2(k2))
            end do
        end do

        call construct_splines_2d(x1, x2, y, spline_order, periodic, spl)

        x_eval(1) = (x1(30) + x1(31))/2.0d0
        x_eval(2) = (x2(28) + x2(29))/2.0d0

        expected = cos(x_eval(1))*cos(x_eval(2))

        call evaluate_splines_2d(x_eval, spl, actual)

        if (abs(expected - actual) > TOL) error stop

        call destroy_splines_2d(spl)

    end subroutine test_spline_2d


    subroutine test_spline_3d(spline_order, periodic)
            use interpolate

            integer, intent(in) :: spline_order(3)
            logical, intent(in) :: periodic(3)


            integer, parameter :: N_POINTS(3) = (/82, 93, 87/)

            real(8), allocatable :: x1(:), x2(:), x3(:), y(:,:,:)
            integer :: k1, k2, k3

            real(8) :: x_eval(3), expected, actual

            type(SplineData3D) :: spl

            allocate(x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3)))
            allocate(y(N_POINTS(1), N_POINTS(2), N_POINTS(3)))

            call linspace(0.0d0, 2.0d0 * pi, N_POINTS(1), x1)
            call linspace(0.0d0, 2.0d0 * pi, N_POINTS(2), x2)

            do k3 = 1, N_POINTS(2)
                do k2 = 1, N_POINTS(2)
                    do k1 = 1, N_POINTS(1)
                       y(k1, k2, k3) = cos(x1(k1))*cos(x2(k2))*cos(x3(k3))
                  end do
                end do
            end do

            call construct_splines_3d(x1, x2, x3, y, spline_order, periodic, spl)

            x_eval(1) = (x1(30) + x1(31))/2.0d0
            x_eval(2) = (x2(28) + x2(29))/2.0d0
            x_eval(2) = (x2(27) + x2(28))/2.0d0

            expected = cos(x_eval(1))*cos(x_eval(2))*cos(x_eval(3))

            call evaluate_splines_3d(x_eval, spl, actual)

            if (abs(expected - actual) > TOL) error stop

            call destroy_splines_3d(spl)

        end subroutine test_spline_3d

end program test_interpolate
