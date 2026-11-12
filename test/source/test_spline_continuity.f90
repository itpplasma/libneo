!> Every equidistant spline builder (orders 3 to 5, regular and periodic, and
!> the complex fifth-order one) must return a piecewise polynomial that
!> interpolates the data and is continuous up to derivative order ns - 1 at
!> each interior knot, and across the seam for periodic splines.
program test_spline_continuity
    use libneo_kinds, only: dp, cdp
    use spl_three_to_five_sub, only: spl_reg, spl_per
    implicit none

    interface
        subroutine cspl_five_reg(n, h, a, b, c, d, e, f)
            import :: dp, cdp
            integer, intent(in) :: n
            real(dp), intent(in) :: h
            complex(cdp), dimension(n), intent(in) :: a
            complex(cdp), dimension(n), intent(out) :: b, c, d, e, f
        end subroutine cspl_five_reg
    end interface

    integer, parameter :: n = 41
    real(dp), parameter :: tol = 1.0e-7_dp
    logical :: ok
    integer :: ns

    ok = .true.
    do ns = 3, 5
        call check(ns, .false., ok)
        call check(ns, .true., ok)
    end do
    call check_complex(ok)
    if (.not. ok) error stop "spline builder is not smooth"

contains

    subroutine check(ns, periodic, ok)
        integer, intent(in) :: ns
        logical, intent(in) :: periodic
        logical, intent(inout) :: ok

        real(dp) :: coef(0:ns, n), data(n), h, x, worst
        integer :: i

        h = 2.0_dp*acos(-1.0_dp)/real(n - 1, dp)
        do i = 1, n
            x = real(i - 1, dp)*h
            data(i) = cos(x) + 0.25_dp*sin(3.0_dp*x)
            if (.not. periodic) data(i) = data(i) + 0.1_dp*x**2
        end do
        coef = 0.0_dp
        coef(0, :) = data
        if (periodic) then
            call spl_per(ns, n, h, coef)
        else
            call spl_reg(ns, n, h, coef)
        end if
        worst = max(maxval(abs(coef(0, :) - data)), max_jump(ns, h, coef, periodic))
        if (worst > tol) then
            write (*, '(a, i2, l2, es12.3)') "FAIL order, periodic, jump", ns, &
                periodic, worst
            ok = .false.
        end if
    end subroutine check

    subroutine check_complex(ok)
        logical, intent(inout) :: ok

        complex(cdp) :: a(n), b(n), c(n), d(n), e(n), f(n)
        real(dp) :: coef(0:5, n), h, x, worst
        integer :: i

        h = 2.0_dp*acos(-1.0_dp)/real(n - 1, dp)
        do i = 1, n
            x = real(i - 1, dp)*h
            a(i) = cmplx(cos(x) + 0.1_dp*x**2, sin(2.0_dp*x) - x, cdp)
        end do
        call cspl_five_reg(n, h, a, b, c, d, e, f)
        coef(0, :) = real(a, dp)
        coef(1, :) = real(b, dp)
        coef(2, :) = real(c, dp)
        coef(3, :) = real(d, dp)
        coef(4, :) = real(e, dp)
        coef(5, :) = real(f, dp)
        worst = max_jump(5, h, coef, .false.)
        coef(0, :) = aimag(a)
        coef(1, :) = aimag(b)
        coef(2, :) = aimag(c)
        coef(3, :) = aimag(d)
        coef(4, :) = aimag(e)
        coef(5, :) = aimag(f)
        worst = max(worst, max_jump(5, h, coef, .false.))
        if (worst > tol) then
            write (*, '(a, es12.3)') "FAIL cspl_five_reg jump", worst
            ok = .false.
        end if
    end subroutine check_complex

    !> Largest jump of derivative k = 0..ns-1 between the end of one interval
    !> and the start of the next, relative to the size of that derivative.
    function max_jump(ns, h, coef, periodic) result(worst)
        integer, intent(in) :: ns
        real(dp), intent(in) :: h, coef(0:, :)
        logical, intent(in) :: periodic
        real(dp) :: worst

        real(dp) :: left, scale
        integer :: i, k, m, last, next

        worst = 0.0_dp
        last = n - 2
        if (periodic) last = n - 1
        do k = 0, ns - 1
            scale = max(1.0_dp, maxval(abs(coef(k, :)))*factorial(k))
            do i = 1, last
                next = i + 1
                if (i == n - 1) next = 1
                left = 0.0_dp
                do m = k, ns
                    left = left + coef(m, i)*factorial(m)/factorial(m - k)*h**(m - k)
                end do
                worst = max(worst, abs(left - coef(k, next)*factorial(k))/scale)
            end do
        end do
    end function max_jump

    pure real(dp) function factorial(k)
        integer, intent(in) :: k
        integer :: j

        factorial = 1.0_dp
        do j = 2, k
            factorial = factorial*real(j, dp)
        end do
    end function factorial

end program test_spline_continuity
