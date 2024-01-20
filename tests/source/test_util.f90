program test_util
    use iso_fortran_env, only: dp => real64

    implicit none

    real(dp), parameter :: TOL = 1.0d-12

    call test_linspace

contains

    subroutine test_linspace
        use util, only : linspace

        real(8), dimension(5) :: expected, actual

        expected = (/ 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)
        call linspace(0.0d0, 4.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > TOL) error stop
    end subroutine test_linspace

end program test_util
