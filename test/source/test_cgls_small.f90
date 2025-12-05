program test_cgls_small
    use libneo_kinds, only : dp
    use cgls_small, only : cgls_solve_small

    implicit none

    real(dp), parameter :: TOL = 1.0d-12

    call test_spd_2x2()
    call test_spd_3x3()

contains

    subroutine test_spd_2x2()
        integer, parameter :: order = 1
        real(dp) :: M(0:order, 0:order)
        real(dp) :: x_true(0:order)
        real(dp) :: rhs(0:order)
        real(dp) :: x(0:order)
        real(dp) :: err
        integer :: i, j

        M(0,0) = 4.0d0
        M(0,1) = 1.0d0
        M(1,0) = 1.0d0
        M(1,1) = 3.0d0

        x_true(0) = 1.0d0
        x_true(1) = 2.0d0

        rhs = 0.0d0
        do i = 0, order
            do j = 0, order
                rhs(i) = rhs(i) + M(i,j)*x_true(j)
            end do
        end do

        call cgls_solve_small(order, M, rhs, x)

        err = 0.0d0
        do i = 0, order
            err = err + (x(i) - x_true(i))**2
        end do

        if (sqrt(err) > TOL) then
            error stop "cgls_small 2x2 SPD test failed"
        end if

    end subroutine test_spd_2x2


    subroutine test_spd_3x3()
        integer, parameter :: order = 2
        real(dp) :: M(0:order, 0:order)
        real(dp) :: x_true(0:order)
        real(dp) :: rhs(0:order)
        real(dp) :: x(0:order)
        real(dp) :: err
        integer :: i, j

        M(0,0) = 4.0d0
        M(0,1) = 1.0d0
        M(0,2) = 0.0d0
        M(1,0) = 1.0d0
        M(1,1) = 3.0d0
        M(1,2) = 1.0d0
        M(2,0) = 0.0d0
        M(2,1) = 1.0d0
        M(2,2) = 2.0d0

        x_true(0) = 1.0d0
        x_true(1) = -1.0d0
        x_true(2) = 0.5d0

        rhs = 0.0d0
        do i = 0, order
            do j = 0, order
                rhs(i) = rhs(i) + M(i,j)*x_true(j)
            end do
        end do

        call cgls_solve_small(order, M, rhs, x)

        err = 0.0d0
        do i = 0, order
            err = err + (x(i) - x_true(i))**2
        end do

        if (sqrt(err) > TOL) then
            error stop "cgls_small 3x3 SPD test failed"
        end if

    end subroutine test_spd_3x3

end program test_cgls_small

