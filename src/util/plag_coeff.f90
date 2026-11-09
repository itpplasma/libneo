module plag_coeff_sub
    use libneo_kinds, only: dp
    implicit none

    private
    public :: plag_coeff, plag_coeff_uniform, lagrange_stencil

contains

    !> Lagrange weights at x for the nodes xp(1:npoi).
    !> sum(f(1:npoi)*coef(0, 1:npoi)) is the interpolant at x and, for
    !> nder >= 1, sum(f(1:npoi)*coef(1, 1:npoi)) its first derivative.
    !> Rows above 1 are left untouched.
    pure subroutine plag_coeff(npoi, nder, x, xp, coef)
        integer, intent(in) :: npoi, nder
        real(dp), intent(in) :: x
        real(dp), dimension(npoi), intent(in) :: xp
        real(dp), dimension(0:nder, npoi), intent(out) :: coef

        ! prod_{k /= i} (x - xp(k)) and its x derivative from prefix and
        ! suffix products, so the cost grows like npoi**2 only in the
        ! node-only denominators.
        real(dp) :: pre(npoi), dpre(npoi), suf(npoi), dsuf(npoi), den
        integer :: i, k

        pre(1) = 1.0_dp
        dpre(1) = 0.0_dp
        do i = 2, npoi
            pre(i) = pre(i - 1)*(x - xp(i - 1))
            dpre(i) = dpre(i - 1)*(x - xp(i - 1)) + pre(i - 1)
        end do
        suf(npoi) = 1.0_dp
        dsuf(npoi) = 0.0_dp
        do i = npoi - 1, 1, -1
            suf(i) = suf(i + 1)*(x - xp(i + 1))
            dsuf(i) = dsuf(i + 1)*(x - xp(i + 1)) + suf(i + 1)
        end do
        do i = 1, npoi
            den = 1.0_dp
            do k = 1, npoi
                if (k /= i) den = den*(xp(i) - xp(k))
            end do
            coef(0, i) = pre(i)*suf(i)/den
            if (nder >= 1) coef(1, i) = (dpre(i)*suf(i) + pre(i)*dsuf(i))/den
        end do
    end subroutine plag_coeff

    !> plag_coeff for the uniform nodes x1 + (i - 1)/hm1, i = 1..npoi; hm1 is
    !> the reciprocal mesh step. The node denominators are the integers
    !> prod_{k /= i} (i - k), so no division depends on the mesh.
    pure subroutine plag_coeff_uniform(npoi, nder, x, x1, hm1, coef)
        integer, intent(in) :: npoi, nder
        real(dp), intent(in) :: x, x1, hm1
        real(dp), dimension(0:nder, npoi), intent(out) :: coef

        real(dp) :: t, pre(npoi), dpre(npoi), suf(npoi), dsuf(npoi), rden
        integer :: i, k, den

        t = (x - x1)*hm1
        pre(1) = 1.0_dp
        dpre(1) = 0.0_dp
        do i = 2, npoi
            pre(i) = pre(i - 1)*(t - (i - 2))
            dpre(i) = dpre(i - 1)*(t - (i - 2)) + pre(i - 1)
        end do
        suf(npoi) = 1.0_dp
        dsuf(npoi) = 0.0_dp
        do i = npoi - 1, 1, -1
            suf(i) = suf(i + 1)*(t - i)
            dsuf(i) = dsuf(i + 1)*(t - i) + suf(i + 1)
        end do
        do i = 1, npoi
            den = 1
            do k = 1, npoi
                if (k /= i) den = den*(i - k)
            end do
            rden = 1.0_dp/real(den, dp)
            coef(0, i) = pre(i)*suf(i)*rden
            if (nder >= 1) coef(1, i) = (dpre(i)*suf(i) + pre(i)*dsuf(i))*rden*hm1
        end do
    end subroutine plag_coeff_uniform

    !> Consecutive node indices of a stencil of size(indu) starting at first,
    !> shifted inside 1..nup at the mesh ends.
    pure subroutine lagrange_stencil(first, nup, indu)
        integer, intent(in) :: first, nup
        integer, intent(out) :: indu(:)

        integer :: i, start

        start = max(first, 1)
        if (start + size(indu) - 1 > nup) start = nup - size(indu) + 1
        do i = 1, size(indu)
            indu(i) = start + i - 1
        end do
    end subroutine lagrange_stencil

end module plag_coeff_sub
