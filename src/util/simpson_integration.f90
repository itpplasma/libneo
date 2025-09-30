module simpson_integration

    use, intrinsic :: iso_fortran_env, only : dp => real64

    implicit none

contains

    subroutine simpson_nonequi(integral, x, f)
        real(dp), intent(out) :: integral
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: f(:)

        real(dp) :: prefix(size(x))

        call simpson_prefix(prefix, x, f)
        integral = prefix(size(prefix))
    end subroutine simpson_nonequi

    subroutine simpson_prefix(prefix, x, f)
        real(dp), intent(out) :: prefix(:)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: f(:)

        integer :: npts, nmax, i
        real(dp) :: xi, xi1, xi2
        real(dp) :: w0_total, w1_total, w2_total
        real(dp) :: w0_first, w1_first, w2_first
        real(dp) :: w0_second, w1_second, w2_second
        real(dp) :: hlast
        logical :: has_extra

        if (size(prefix) /= size(x) .or. size(prefix) /= size(f)) then
            error stop 'simpson_prefix: inconsistent array sizes'
        end if

        npts = size(x)
        prefix = 0.0_dp

        if (npts < 2) then
            error stop 'simpson_prefix: need at least two points'
        end if
        if (size(f) /= npts) then
            error stop 'simpson_prefix: size(x) /= size(f)'
        end if

        do i = 1, npts - 1
            if (x(i + 1) <= x(i)) then
                error stop 'simpson_prefix: x must be strictly increasing'
            end if
        end do

        if (mod(npts - 1, 2) == 0) then
            nmax = npts
            has_extra = .false.
        else
            nmax = npts - 1
            has_extra = .true.
        end if

        do i = 1, nmax - 2, 2
            xi = x(i)
            xi1 = x(i + 1)
            xi2 = x(i + 2)

            call weights_interval(xi, xi1, xi2, xi, xi2, w0_total, w1_total, w2_total)
            call weights_interval(xi, xi1, xi2, xi, xi1, w0_first, w1_first, w2_first)

            w0_second = w0_total - w0_first
            w1_second = w1_total - w1_first
            w2_second = w2_total - w2_first

            prefix(i + 1) = prefix(i) + w0_first * f(i) + w1_first * f(i + 1) + w2_first * f(i + 2)
            prefix(i + 2) = prefix(i + 1) + w0_second * f(i) + w1_second * f(i + 1) + w2_second * f(i + 2)
        end do

        if (has_extra) then
            hlast = x(npts) - x(npts - 1)
            prefix(npts) = prefix(nmax) + 0.5_dp * (f(npts - 1) + f(npts)) * hlast
        else
            prefix(npts) = prefix(nmax)
        end if
    end subroutine simpson_prefix

    subroutine weights_interval(x0, x1, x2, a, b, w0, w1, w2)
        real(dp), intent(in) :: x0, x1, x2
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: w0, w1, w2

        real(dp) :: d0, d1, d2

        d0 = (x0 - x1) * (x0 - x2)
        d1 = (x1 - x0) * (x1 - x2)
        d2 = (x2 - x0) * (x2 - x1)

        w0 = F(b, x1, x2) / d0 - F(a, x1, x2) / d0
        w1 = F(b, x0, x2) / d1 - F(a, x0, x2) / d1
        w2 = F(b, x0, x1) / d2 - F(a, x0, x1) / d2
    contains
        pure real(dp) function F(x, p, q)
            real(dp), intent(in) :: x, p, q
            F = x**3 / 3.0_dp - (p + q) * x**2 / 2.0_dp + p * q * x
        end function F
    end subroutine weights_interval

end module simpson_integration
