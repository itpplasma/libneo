module cgls_small
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none

contains

    subroutine cgls_solve_small(order, amat, rhs, x)
        integer, intent(in) :: order
        real(dp), intent(in) :: amat(0:order, 0:order)
        real(dp), intent(in) :: rhs(0:order)
        real(dp), intent(out) :: x(0:order)

        integer :: n, i, j, k, iter
        real(dp) :: alpha, beta
        real(dp) :: r_norm2, r_norm2_new, denom
        real(dp) :: tol

        real(dp), allocatable :: r(:), s(:), p(:), q(:)

        n = order + 1
        tol = 1.0d-12

        allocate(r(0:order))
        allocate(s(0:order))
        allocate(p(0:order))
        allocate(q(0:order))

        x = 0.0d0

        r = rhs

        do i = 0, order
            s(i) = 0.0d0
            do j = 0, order
                s(i) = s(i) + amat(j, i)*r(j)
            end do
        end do

        p = s
        r_norm2 = 0.0d0
        do i = 0, order
            r_norm2 = r_norm2 + s(i)*s(i)
        end do

        if (r_norm2 <= tol) then
            deallocate(r, s, p, q)
            return
        end if

        do iter = 1, 64

            do i = 0, order
                q(i) = 0.0d0
                do j = 0, order
                    q(i) = q(i) + amat(i, j)*p(j)
                end do
            end do

            denom = 0.0d0
            do i = 0, order
                denom = denom + q(i)*q(i)
            end do

            if (denom <= tol) exit

            alpha = r_norm2/denom

            do i = 0, order
                x(i) = x(i) + alpha*p(i)
                r(i) = r(i) - alpha*q(i)
            end do

            do i = 0, order
                s(i) = 0.0d0
                do j = 0, order
                    s(i) = s(i) + amat(j, i)*r(j)
                end do
            end do

            r_norm2_new = 0.0d0
            do i = 0, order
                r_norm2_new = r_norm2_new + s(i)*s(i)
            end do

            if (r_norm2_new <= tol*r_norm2) exit

            beta = r_norm2_new/r_norm2
            r_norm2 = r_norm2_new

            do i = 0, order
                p(i) = s(i) + beta*p(i)
            end do
        end do

        deallocate(r, s, p, q)

    end subroutine cgls_solve_small

end module cgls_small
