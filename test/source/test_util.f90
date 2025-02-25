program test_util
    use libneo_kinds, only : dp

    implicit none

    real(dp), parameter :: TOL = 1.0d-12

    call test_linspace

contains

    subroutine test_linspace
        use libneo_util, only : linspace

        real(dp), dimension(5) :: expected, actual

        expected = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
        call linspace(0.0d0, 4.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > TOL) error stop
    end subroutine test_linspace

end program test_util
