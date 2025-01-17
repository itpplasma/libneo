module libneo_util
    use iso_fortran_env, only: dp => real64

    implicit none

contains

    subroutine linspace(a, b, n, x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp), dimension(:), intent(out) :: x

        real(dp) :: dx
        integer :: i

        dx = (b - a) / (n - 1)
        do i = 1, n
            x(i) = a + (i - 1) * dx
        end do
    end subroutine linspace


    function get_random_numbers(xmin, xmax, n, seed) result(x)
        real(dp), intent(in) :: xmin, xmax
        integer, intent(in) :: n
        integer, dimension(:), intent(in), optional :: seed
        real(dp), dimension(:), allocatable :: x

        if (present(seed)) then
            call random_seed(put=seed)
        end if
        allocate(x(n))
        call random_number(x)
        x = xmin + (xmax - xmin) * x
    end function get_random_numbers

end module libneo_util
