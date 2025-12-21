module acc_tridiag_pcr
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: tridiag_pcr_solve_inplace
    public :: tridiag_pcr_solve_cyclic

contains

    subroutine tridiag_pcr_solve_inplace(a, b, c, d, a_tmp, b_tmp, c_tmp, d_tmp, x)
        real(dp), intent(inout) :: a(:, :)     ! (n, nsys)
        real(dp), intent(inout) :: b(:, :)
        real(dp), intent(inout) :: c(:, :)
        real(dp), intent(inout) :: d(:, :)
        real(dp), intent(inout) :: a_tmp(:, :)
        real(dp), intent(inout) :: b_tmp(:, :)
        real(dp), intent(inout) :: c_tmp(:, :)
        real(dp), intent(inout) :: d_tmp(:, :)
        real(dp), intent(out) :: x(:, :)

        integer :: n, nsys
        integer :: step, stride
        integer :: i, jsys
        integer :: im, ip
        real(dp) :: alpha, gamma, denom

        n = size(b, 1)
        nsys = size(b, 2)

        if (size(a, 1) /= n .or. size(c, 1) /= n .or. size(d, 1) /= n) then
            error stop "tridiag_pcr_solve_inplace: inconsistent first dimension sizes"
        end if
        if (size(a, 2) /= nsys .or. size(c, 2) /= nsys .or. size(d, 2) /= nsys) then
            error stop "tridiag_pcr_solve_inplace: inconsistent system counts"
        end if
        if (size(a_tmp, 1) /= n .or. size(b_tmp, 1) /= n .or. size(c_tmp, 1) /= n) then
            error stop "tridiag_pcr_solve_inplace: tmp size mismatch"
        end if
        if (size(d_tmp, 1) /= n .or. size(a_tmp, 2) /= nsys .or. size(b_tmp, 2) /= nsys) then
            error stop "tridiag_pcr_solve_inplace: tmp system count mismatch"
        end if
        if (size(x, 1) /= n .or. size(x, 2) /= nsys) then
            error stop "tridiag_pcr_solve_inplace: output size mismatch"
        end if

        stride = 1
        do step = 1, 64
            if (stride >= n) exit

            !$acc parallel loop collapse(2) gang vector present(a, b, c, d, a_tmp, &
            !$acc& b_tmp, c_tmp, d_tmp) private(i, jsys, im, ip, alpha, gamma, denom)
            do jsys = 1, nsys
                do i = 1, n
                    im = i - stride
                    ip = i + stride

                    alpha = 0.0d0
                    if (im >= 1) alpha = a(i, jsys)/b(im, jsys)

                    gamma = 0.0d0
                    if (ip <= n) gamma = c(i, jsys)/b(ip, jsys)

                    denom = b(i, jsys)
                    if (im >= 1) denom = denom - alpha*c(im, jsys)
                    if (ip <= n) denom = denom - gamma*a(ip, jsys)
                    b_tmp(i, jsys) = denom

                    d_tmp(i, jsys) = d(i, jsys)
                    if (im >= 1) d_tmp(i, jsys) = d_tmp(i, jsys) - alpha*d(im, jsys)
                    if (ip <= n) d_tmp(i, jsys) = d_tmp(i, jsys) - gamma*d(ip, jsys)

                    a_tmp(i, jsys) = 0.0d0
                    if (im >= 1) a_tmp(i, jsys) = -alpha*a(im, jsys)

                    c_tmp(i, jsys) = 0.0d0
                    if (ip <= n) c_tmp(i, jsys) = -gamma*c(ip, jsys)
                end do
            end do

            !$acc parallel loop collapse(2) gang vector present(a, b, c, d, a_tmp, &
            !$acc& b_tmp, c_tmp, d_tmp)
            do jsys = 1, nsys
                do i = 1, n
                    a(i, jsys) = a_tmp(i, jsys)
                    b(i, jsys) = b_tmp(i, jsys)
                    c(i, jsys) = c_tmp(i, jsys)
                    d(i, jsys) = d_tmp(i, jsys)
                end do
            end do

            stride = stride*2
        end do

        !$acc parallel loop collapse(2) gang vector present(b, d, x)
        do jsys = 1, nsys
            do i = 1, n
                x(i, jsys) = d(i, jsys)/b(i, jsys)
            end do
        end do
    end subroutine tridiag_pcr_solve_inplace

    subroutine tridiag_pcr_solve_cyclic(bdiag, rhs, a, b, c, d, a_tmp, b_tmp, c_tmp, &
                                        d_tmp, y, z, u, factor, x)
        real(dp), intent(in) :: bdiag(:, :)   ! (n, nsys)
        real(dp), intent(in) :: rhs(:, :)

        real(dp), intent(inout) :: a(:, :)
        real(dp), intent(inout) :: b(:, :)
        real(dp), intent(inout) :: c(:, :)
        real(dp), intent(inout) :: d(:, :)
        real(dp), intent(inout) :: a_tmp(:, :)
        real(dp), intent(inout) :: b_tmp(:, :)
        real(dp), intent(inout) :: c_tmp(:, :)
        real(dp), intent(inout) :: d_tmp(:, :)
        real(dp), intent(inout) :: y(:, :)
        real(dp), intent(inout) :: z(:, :)
        real(dp), intent(inout) :: u(:, :)
        real(dp), intent(inout) :: factor(:)
        real(dp), intent(out) :: x(:, :)

        integer :: n, nsys
        integer :: i, jsys
        real(dp), parameter :: alpha = 1.0d0
        real(dp), parameter :: gamma = 1.0d0
        real(dp) :: beta, denom

        n = size(bdiag, 1)
        nsys = size(bdiag, 2)

        if (size(rhs, 1) /= n .or. size(rhs, 2) /= nsys) then
            error stop "tridiag_pcr_solve_cyclic: rhs size mismatch"
        end if
        if (size(factor) /= nsys) then
            error stop "tridiag_pcr_solve_cyclic: factor size mismatch"
        end if

        !$acc parallel loop gang vector present(u, bdiag)
        do jsys = 1, nsys
            u(1, jsys) = -bdiag(1, jsys)
            u(n, jsys) = alpha
        end do
        !$acc parallel loop collapse(2) gang vector present(u)
        do jsys = 1, nsys
            do i = 2, n - 1
                u(i, jsys) = 0.0d0
            end do
        end do

        !$acc parallel loop collapse(2) gang vector present(a, b, c, d, bdiag, rhs)
        do jsys = 1, nsys
            do i = 1, n
                a(i, jsys) = 1.0d0
                b(i, jsys) = bdiag(i, jsys)
                c(i, jsys) = 1.0d0
                d(i, jsys) = rhs(i, jsys)
            end do
        end do
        !$acc parallel loop gang vector present(a, b, c, d, bdiag)
        do jsys = 1, nsys
            beta = -bdiag(1, jsys)
            a(1, jsys) = 0.0d0
            c(n, jsys) = 0.0d0
            b(1, jsys) = bdiag(1, jsys) - beta
            b(n, jsys) = bdiag(n, jsys) - alpha*gamma/beta
        end do
        call tridiag_pcr_solve_inplace(a, b, c, d, a_tmp, b_tmp, c_tmp, d_tmp, y)

        !$acc parallel loop collapse(2) gang vector present(a, b, c, d, bdiag, u)
        do jsys = 1, nsys
            do i = 1, n
                a(i, jsys) = 1.0d0
                b(i, jsys) = bdiag(i, jsys)
                c(i, jsys) = 1.0d0
                d(i, jsys) = u(i, jsys)
            end do
        end do
        !$acc parallel loop gang vector present(a, b, c, bdiag)
        do jsys = 1, nsys
            beta = -bdiag(1, jsys)
            a(1, jsys) = 0.0d0
            c(n, jsys) = 0.0d0
            b(1, jsys) = bdiag(1, jsys) - beta
            b(n, jsys) = bdiag(n, jsys) - alpha*gamma/beta
        end do
        call tridiag_pcr_solve_inplace(a, b, c, d, a_tmp, b_tmp, c_tmp, d_tmp, z)

        !$acc parallel loop gang vector present(bdiag, y, z, factor) private(beta, denom)
        do jsys = 1, nsys
            beta = -bdiag(1, jsys)
            denom = 1.0d0 + z(1, jsys) + gamma*z(n, jsys)/beta
            factor(jsys) = (y(1, jsys) + gamma*y(n, jsys)/beta) / denom
        end do

        !$acc parallel loop collapse(2) gang vector present(y, z, factor, x)
        do jsys = 1, nsys
            do i = 1, n
                x(i, jsys) = y(i, jsys) - factor(jsys)*z(i, jsys)
            end do
        end do
    end subroutine tridiag_pcr_solve_cyclic

end module acc_tridiag_pcr
