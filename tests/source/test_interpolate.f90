program test_interpolate
    use interpolate

    implicit none

    integer, parameter :: N_POINTS = 100
    real(8), parameter :: X_MIN = 0.0d0, X_MAX = 2.0d0 * pi

    real(8), dimension(N_POINTS) :: x, y, dy, d2y

    call linspace(x, 0.0d0, 2.0d0 * pi, 100)

    call test_interpolate_1d

contains

    subroutine generate_test_data_1d(x, y, dy, d2y)
        real, dimension(:), intent(in)  :: x
        real, dimension(:), intent(out) :: y, dy, d2y
        integer :: i

        y = cos(x)
        dy = -sin(x)
        d2y = -cos(x)

    end subroutine generate_test_data_1d

    subroutine test_interpolate_1d

    end subroutine test_interpolate_1d

end program test_interpolate
