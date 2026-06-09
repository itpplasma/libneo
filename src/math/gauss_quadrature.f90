module neo_gauss_quadrature
    ! Gauss quadrature rules via the Golub-Welsch algorithm:
    ! G. H. Golub, J. H. Welsch, Calculation of Gauss quadrature rules,
    ! Math. Comp. 23 (1969) 221-230. Recurrence coefficients from
    ! DLMF 18.9 (https://dlmf.nist.gov/18.9), zeroth moments from
    ! DLMF Table 18.3.1.
    use libneo_kinds, only: dp

    implicit none
    private

    public :: gauss_legendre, gauss_legendre_ab, gauss_gen_laguerre

contains

    pure subroutine gauss_legendre(n, x, w)
        ! Nodes are roots of P_n, converged by Newton iteration from the
        ! asymptotic estimate cos(pi (i - 1/4)/(n + 1/2)), Abramowitz &
        ! Stegun 22.16.6; weights 2/((1 - x^2) P_n'(x)^2), A&S 25.4.29.
        ! Agrees with the Golub-Welsch eigensolve at machine precision.
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        real(dp), parameter :: pi = 3.14159265358979324_dp
        integer, parameter :: max_iter = 16

        integer :: m, i, it
        real(dp) :: xh((n + 1)/2), p((n + 1)/2), dpdx((n + 1)/2)
        real(dp) :: dx((n + 1)/2)

        if (n < 1) then
            error stop "neo_gauss_quadrature: n must be >= 1"
        end if
        m = (n + 1)/2
        do i = 1, m
            xh(i) = cos(pi*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
        end do
        do it = 1, max_iter
            call legendre_eval(n, m, xh, p, dpdx)
            dx = p/dpdx
            xh = xh - dx
            if (maxval(abs(dx)) <= 4.0_dp*epsilon(1.0_dp)) exit
        end do
        call legendre_eval(n, m, xh, p, dpdx)
        do i = 1, m
            x(n + 1 - i) = xh(i)
            x(i) = -xh(i)
            w(i) = 2.0_dp/((1.0_dp - xh(i)*xh(i))*dpdx(i)*dpdx(i))
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre

    pure subroutine legendre_eval(n, m, x, p, dpdx)
        ! Three-term recurrence DLMF 18.9.8 and derivative DLMF 18.9.17,
        ! evaluated for all m nodes at once to allow vectorization
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(m)
        real(dp), intent(out) :: p(m), dpdx(m)

        real(dp) :: p_prev(m), p_next(m), c_cur, c_prev
        integer :: k

        p_prev = 1.0_dp
        p = x
        do k = 1, n - 1
            c_cur = real(2*k + 1, dp)/real(k + 1, dp)
            c_prev = real(k, dp)/real(k + 1, dp)
            p_next = c_cur*x*p - c_prev*p_prev
            p_prev = p
            p = p_next
        end do
        dpdx = real(n, dp)*(x*p - p_prev)/(x*x - 1.0_dp)
    end subroutine legendre_eval

    subroutine gauss_legendre_ab(n, a, b, x, w)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: x(n), w(n)

        real(dp) :: half_length, midpoint

        call gauss_legendre(n, x, w)
        half_length = 0.5_dp*(b - a)
        midpoint = 0.5_dp*(a + b)
        x = midpoint + half_length*x
        w = half_length*w
    end subroutine gauss_legendre_ab

    subroutine gauss_gen_laguerre(n, alpha, x, w)
        ! Weight function x**alpha * exp(-x) on [0, infinity), alpha > -1
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha
        real(dp), intent(out) :: x(n), w(n)

        real(dp) :: diag(n), offdiag(n)
        integer :: k

        if (alpha <= -1.0_dp) then
            error stop "gauss_gen_laguerre: alpha must be > -1"
        end if
        do k = 1, n
            diag(k) = 2.0_dp*real(k - 1, dp) + alpha + 1.0_dp
        end do
        do k = 1, n - 1
            offdiag(k) = sqrt(real(k, dp)*(real(k, dp) + alpha))
        end do
        call golub_welsch(n, diag, offdiag, gamma(alpha + 1.0_dp), x, w)
    end subroutine gauss_gen_laguerre

    subroutine golub_welsch(n, diag, offdiag, mu0, x, w)
        ! Nodes are eigenvalues of the symmetric tridiagonal Jacobi matrix,
        ! DLMF 3.5.31 (LAPACK dsterf); weights are the Christoffel numbers
        ! w_i = 1 / sum_k q_k(x_i)^2 over the orthonormal polynomials q_k
        ! (Szego, Orthogonal Polynomials, 1975, Thm 3.4.2), equivalent to
        ! the first-eigenvector-component form DLMF 3.5.32
        integer, intent(in) :: n
        real(dp), intent(in) :: diag(n), offdiag(n)
        real(dp), intent(in) :: mu0
        real(dp), intent(out) :: x(n), w(n)

        interface
            subroutine dsterf(n, d, e, info)
                import :: dp
                integer, intent(in) :: n
                real(dp), intent(inout) :: d(*), e(*)
                integer, intent(out) :: info
            end subroutine dsterf
        end interface

        real(dp) :: e(max(1, n - 1)), p_prev, p_cur, p_next, b_prev, s
        integer :: info, i, k

        if (n < 1) then
            error stop "neo_gauss_quadrature: n must be >= 1"
        end if
        x = diag
        e(1:n - 1) = offdiag(1:n - 1)
        call dsterf(n, x, e, info)
        if (info /= 0) then
            error stop "neo_gauss_quadrature: dsterf failed to converge"
        end if
        do i = 1, n
            p_prev = 0.0_dp
            p_cur = 1.0_dp/sqrt(mu0)
            b_prev = 0.0_dp
            s = p_cur*p_cur
            do k = 1, n - 1
                p_next = ((x(i) - diag(k))*p_cur - b_prev*p_prev)/offdiag(k)
                p_prev = p_cur
                p_cur = p_next
                b_prev = offdiag(k)
                s = s + p_cur*p_cur
            end do
            w(i) = 1.0_dp/s
        end do
    end subroutine golub_welsch

end module neo_gauss_quadrature
