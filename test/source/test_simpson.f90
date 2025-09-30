program test_simpson
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use simpson_integration, only : simpson_nonequi, simpson_prefix

    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: n = 55

    real(dp) :: approx
    real(dp) :: prefix(n)
    real(dp) :: x(n), f(n)
    integer :: i

    do i = 1, n
        x(i) = 0.5_dp * pi * (1.0_dp - cos((i - 1.0_dp) * pi / (n - 1.0_dp)))
        f(i) = sin(x(i))
    end do

    call simpson_prefix(prefix, x, f)
    if (abs(prefix(n) - 2.0_dp) > 1.0e-6_dp) then
        write(*,*) 'Simpson prefix failed. Value:', prefix(n)
        error stop
    end if

    call simpson_nonequi(approx, x, f)
    if (abs(approx - 2.0_dp) > 1.0e-6_dp) then
        write(*,*) 'Simpson total integral failed. Value:', approx
        error stop
    end if

end program test_simpson
