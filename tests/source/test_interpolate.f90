program test_interpolate
    use math_constants
    use util

    implicit none

    integer, parameter :: N_POINTS = 100
    real(8), parameter :: X_MIN = 0.0d0, X_MAX = TWOPI

    call test_spl_reg(spline_order=3)

contains

    subroutine generate_test_data_1d(x, y, dy, d2y)
        real(8), dimension(:), intent(in)  :: x
        real(8), dimension(:), intent(out) :: y, dy, d2y

        y = cos(x)
        dy = -sin(x)
        d2y = -cos(x)

    end subroutine generate_test_data_1d

    subroutine test_spl_reg(spline_order)

        integer, intent(in) :: spline_order

        real(8), dimension(N_POINTS) :: x, y, dy, d2y

        real(8) :: expected, actual

        real(kind=real_kind), dimension(0:spline_order, N_POINTS) :: spline_coeff

        real(8) :: h_step, x_norm
        integer :: i, k, point_index

        call linspace(0.0d0, 2.0d0 * pi, 100, x)
        call generate_test_data_1d(x, y, dy, d2y)

        h_step = x(2) - x(1)

        spline_coeff(0,:) = y
        call spl_reg(spline_order, N_POINTS, h_step, spline_coeff)

        i = 30

        expected = cos(x(i))

        x_norm = x(i)/h_step
        point_index = max(0, min(N_POINTS-1, int(x_norm)))
        x_norm = (x_norm - dble(point_index))*h_step
        point_index = point_index + 1

        actual = spline_coeff(spline_order+1, point_index)

        do k = spline_order, 0, -1
            actual = spline_coeff(k, point_index) + x_norm*actual
        enddo

        print *, 'expected = ', expected
        print *, 'actual = ', actual

    end subroutine test_spl_reg

end program test_interpolate
