program test_interpolate
    implicit none

    real(8), parameter :: TOL = 1.0d-12

    call test_arange
    call test_linspace

contains

    subroutine test_arange
        use util, only : arange

        real(8), dimension(5) :: expected, actual

        expected = (/ 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)
        call arange(0.0d0, 5.0d0, 1.0d0, actual)

        if (maxval(abs(expected - actual)) > TOL) error stop
    end subroutine test_arange

    subroutine test_linspace
        use util, only : linspace

        real(8), dimension(5) :: expected, actual

        expected = (/ 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)
        call linspace(0.0d0, 4.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > TOL) error stop
    end subroutine test_linspace

end program test_interpolate
