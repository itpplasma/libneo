module util
    implicit none

contains

    subroutine arange(a, b, dx, x)
        real(8), intent(in) :: a, b, dx
        real(8), dimension(:), intent(out) :: x

        integer :: i, n

        n = int((b - a) / dx) + 1
        do i = 1, n
            x(i) = a + (i - 1) * dx
        end do
    end subroutine arange

    subroutine linspace(a, b, n, x)
        real(8), intent(in) :: a, b
        integer, intent(in) :: n
        real(8), dimension(:), intent(out) :: x

        real(8) :: dx
        integer :: i

        dx = (b - a) / (n - 1)
        do i = 1, n
            x(i) = a + (i - 1) * dx
        end do
    end subroutine linspace

end module util
