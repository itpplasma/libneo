program test_interpolate
    use math_constants
    use util

    implicit none

    real(8), parameter :: TOL = 1.0d-6

    integer, parameter :: N_POINTS = 100
    real(8), parameter :: X_MIN = 0.0d0, X_MAX = TWOPI

    call test_spl_reg(spline_order=3, periodic=.False.)

contains

    subroutine generate_test_data_1d(x, y, dy, d2y)
        real(8), dimension(:), intent(in)  :: x
        real(8), dimension(:), intent(out) :: y, dy, d2y

        y = cos(x)
        dy = -sin(x)
        d2y = -cos(x)

    end subroutine generate_test_data_1d

    subroutine test_spl_reg(spline_order, periodic)
        use interpolate

        integer, intent(in) :: spline_order
        logical, intent(in) :: periodic

        real(8), dimension(N_POINTS) :: x, y, dy, d2y

        real(8) :: x_eval, expected, actual

        real(kind=real_kind), dimension(0:spline_order, N_POINTS) :: spline_coeff

        type(SplineData1D) :: spl

        call linspace(0.0d0, 2.0d0 * pi, 100, x)
        call generate_test_data_1d(x, y, dy, d2y)

        call construct_splines_1d(x, y, spline_order, periodic, spl)

        x_eval = (x(30) + x(31))/2.0d0

        expected = cos(x_eval)

        call evaluate_splines_1d(x_eval, spl, actual)

        if (abs(expected - actual) > TOL) error stop

        call destroy_splines_1d(spl)

    end subroutine test_spl_reg

end program test_interpolate
