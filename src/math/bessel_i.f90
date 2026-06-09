module neo_bessel_i
    ! Modified Bessel functions of the first kind I_n(x) for integer order.
    ! Power series: DLMF 10.25.2. Downward (Miller) recurrence with the
    ! normalization e^x = I_0(x) + 2 sum_{k>=1} I_k(x): DLMF 10.29.1 and
    ! Abramowitz & Stegun 9.6.36 (theta = 0). Large-argument asymptotic
    ! expansion: DLMF 10.40.1.

    use libneo_kinds, only: dp

    implicit none
    private

    public :: bessel_in, bessel_in_array

    ! Below this argument the series converges in few terms; above it the
    ! recurrence factor 2k/x stays small enough that rescaling cannot
    ! overflow within a single step.
    real(dp), parameter :: series_x_max = 2.0_dp
    real(dp), parameter :: rescale_limit = 1.0e250_dp
    real(dp), parameter :: rescale_factor = 1.0e-250_dp
    ! Start the downward recurrence high enough that the trial value at the
    ! starting order is negligible: I_m/I_k0 < start_decay_tol, estimated
    ! with the ratio I_{k+1}/I_k ~ x/(k + 1 + sqrt((k+1)^2 + x^2))
    ! (leading term of the continued fraction DLMF 10.33.1 for I).
    real(dp), parameter :: start_decay_tol = 1.0e-18_dp
    real(dp), parameter :: seed = 1.0e-30_dp
    real(dp), parameter :: two_pi = 6.2831853071795864769_dp
    ! The asymptotic expansion DLMF 10.40.1 reaches double precision for
    ! x >= max(20, 0.7 n^2) (verified against 50-digit mpmath references;
    ! worst relative error 5.2e-16 over n <= 31, x <= 700).
    real(dp), parameter :: asym_x_min = 20.0_dp
    real(dp), parameter :: asym_n_factor = 0.7_dp

contains

    elemental function bessel_in(n, x) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: value

        integer :: na
        real(dp) :: ax

        na = abs(n)
        ax = abs(x)
        if (ax <= series_x_max) then
            value = series_in(na, ax)
        else if (ax >= max(asym_x_min, asym_n_factor*real(na, dp)**2)) then
            value = asym_in(na, ax)
        else
            value = miller_in(na, ax)
        end if
        if (x < 0.0_dp .and. mod(na, 2) == 1) value = -value
    end function bessel_in


    pure subroutine bessel_in_array(nmax, x, values)
        integer, intent(in) :: nmax
        real(dp), intent(in) :: x
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: ax

        ax = abs(x)
        if (ax <= series_x_max) then
            call series_in_array(nmax, ax, values)
        else
            call miller_in_array(nmax, ax, values)
        end if
        if (x < 0.0_dp .and. nmax >= 1) then
            values(1:nmax:2) = -values(1:nmax:2)
        end if
    end subroutine bessel_in_array


    elemental function series_in(n, ax) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: pref, halfx
        integer :: j

        halfx = 0.5_dp*ax
        pref = 1.0_dp
        do j = 1, n
            pref = pref*halfx/real(j, dp)
            if (pref == 0.0_dp) exit
        end do
        if (pref == 0.0_dp) then
            value = 0.0_dp
        else
            value = pref*series_sum(n, 0.25_dp*ax*ax)
        end if
    end function series_in


    elemental function series_sum(n, q) result(total)
        ! Sum over k of q^k / (k! (n+k)! / n!), DLMF 10.25.2 with the
        ! prefactor (x/2)^n / n! split off.
        integer, intent(in) :: n
        real(dp), intent(in) :: q
        real(dp) :: total

        real(dp) :: term
        integer :: k

        term = 1.0_dp
        total = 1.0_dp
        k = 0
        do
            k = k + 1
            term = term*q/(real(k, dp)*real(n + k, dp))
            total = total + term
            if (term <= epsilon(1.0_dp)*total) exit
        end do
    end function series_sum


    elemental function miller_in(n, ax) result(value)
        integer, intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: p_hi, p_cur, p_lo, s, res
        integer :: m, k

        m = start_order(max(n, int(ax) + 1), ax)
        p_hi = 0.0_dp
        p_cur = seed
        s = 0.0_dp
        res = 0.0_dp
        do k = m, 1, -1
            p_lo = p_hi + (2.0_dp*real(k, dp)/ax)*p_cur
            s = s + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_limit) then
                p_lo = p_lo*rescale_factor
                p_cur = p_cur*rescale_factor
                s = s*rescale_factor
                res = res*rescale_factor
            end if
            p_hi = p_cur
            p_cur = p_lo
            if (k - 1 == n) res = p_cur
        end do
        s = s + p_cur
        value = (res/s)*exp(ax)
    end function miller_in


    elemental function asym_in(n, ax) result(value)
        ! DLMF 10.40.1: I_n(x) ~ e^x/sqrt(2 pi x) sum_k (-1)^k a_k(n)/x^k
        ! with a_k(n) = (4n^2 - 1)(4n^2 - 9)...(4n^2 - (2k-1)^2)/(k! 8^k).
        integer, intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: mu, term, total, prev, r8x
        integer :: k

        mu = 4.0_dp*real(n, dp)**2
        r8x = 1.0_dp/(8.0_dp*ax)
        term = 1.0_dp
        total = 1.0_dp
        prev = 1.0_dp
        do k = 1, 40
            term = -term*(mu - real(2*k - 1, dp)**2)*r8x/real(k, dp)
            if (abs(term) >= prev) exit
            total = total + term
            prev = abs(term)
            if (prev <= epsilon(1.0_dp)*abs(total)) exit
        end do
        value = exp(ax)/sqrt(two_pi*ax)*total
    end function asym_in


    pure subroutine series_in_array(nmax, ax, values)
        integer, intent(in) :: nmax
        real(dp), intent(in) :: ax
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: pref, halfx, q
        integer :: j

        halfx = 0.5_dp*ax
        q = 0.25_dp*ax*ax
        pref = 1.0_dp
        do j = 0, nmax
            if (j > 0) pref = pref*halfx/real(j, dp)
            if (pref == 0.0_dp) then
                values(j:nmax) = 0.0_dp
                return
            end if
            values(j) = pref*series_sum(j, q)
        end do
    end subroutine series_in_array


    pure subroutine miller_in_array(nmax, ax, values)
        integer, intent(in) :: nmax
        real(dp), intent(in) :: ax
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: p_hi, p_cur, p_lo, s
        integer :: m, k

        m = start_order(max(nmax, int(ax) + 1), ax)
        p_hi = 0.0_dp
        p_cur = seed
        s = 0.0_dp
        do k = m, 1, -1
            p_lo = p_hi + (2.0_dp*real(k, dp)/ax)*p_cur
            s = s + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_limit) then
                p_lo = p_lo*rescale_factor
                p_cur = p_cur*rescale_factor
                s = s*rescale_factor
                if (k <= nmax) values(k:nmax) = values(k:nmax)*rescale_factor
            end if
            p_hi = p_cur
            p_cur = p_lo
            if (k - 1 <= nmax) values(k - 1) = p_cur
        end do
        s = s + p_cur
        values = values*(exp(ax)/s)
    end subroutine miller_in_array


    elemental function start_order(k0, ax) result(m)
        integer, intent(in) :: k0
        real(dp), intent(in) :: ax
        integer :: m

        real(dp) :: decay

        m = k0
        decay = 1.0_dp
        do while (decay > start_decay_tol)
            m = m + 1
            decay = decay*ax/(real(m, dp) + sqrt(real(m, dp)**2 + ax*ax))
        end do
        m = m + 2
    end function start_order

end module neo_bessel_i
