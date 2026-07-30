program test_cgls_dense
    use libneo_kinds, only : dp
    use cgls_dense,   only : cgls_dense_solve

    implicit none

    real(dp), parameter :: TOL = 1.0d-10

    call test_small_unweighted()
    call test_small_weighted()

contains

    subroutine test_small_unweighted()
        real(dp) :: phi(3, 2)
        real(dp) :: x_true(2)
        real(dp) :: rhs(3)
        real(dp) :: x(2)
        real(dp) :: err
        integer :: i, j

        phi = reshape([ &
            1.0d0, 2.0d0, &
            0.0d0, 1.0d0, &
            1.0d0, 0.0d0 ], shape(phi))

        x_true(1) = 0.5d0
        x_true(2) = -1.5d0

        do i = 1, 3
            rhs(i) = 0.0d0
            do j = 1, 2
                rhs(i) = rhs(i) + phi(i, j)*x_true(j)
            end do
        end do

        call cgls_dense_solve(phi, rhs, x)

        err = 0.0d0
        do i = 1, 2
            err = err + (x(i) - x_true(i))**2
        end do
        err = sqrt(err)

        if (err > TOL) then
            error stop "cgls_dense unweighted test failed"
        end if
    end subroutine test_small_unweighted


    subroutine test_small_weighted()
        real(dp) :: phi(3, 2)
        real(dp) :: x_true(2)
        real(dp) :: rhs(3)
        real(dp) :: x(2)
        real(dp) :: w(3)
        real(dp) :: err
        integer :: i, j

        phi(1, 1) = 1.0d0
        phi(1, 2) = 0.0d0
        phi(2, 1) = 1.0d0
        phi(2, 2) = 1.0d0
        phi(3, 1) = 0.0d0
        phi(3, 2) = 1.0d0

        x_true(1) = -0.3d0
        x_true(2) = 2.1d0

        do i = 1, 3
            rhs(i) = 0.0d0
            do j = 1, 2
                rhs(i) = rhs(i) + phi(i, j)*x_true(j)
            end do
        end do

        w(1) = 1.0d0
        w(2) = 0.5d0
        w(3) = 2.0d0

        call cgls_dense_solve(phi, rhs, x, w)

        err = 0.0d0
        do i = 1, 2
            err = err + (x(i) - x_true(i))**2
        end do
        err = sqrt(err)

        if (err > TOL) then
            print *, "cgls_dense weighted err =", err
            print *, "  x_true =", x_true
            print *, "  x      =", x
            error stop "cgls_dense weighted test failed"
        end if
    end subroutine test_small_weighted

end program test_cgls_dense
