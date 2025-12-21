module spline_build_lines
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: spl_build_line_inplace
    public :: spl_three_reg_line
    public :: spl_three_per_line
    public :: spl_four_reg_line
    public :: spl_four_per_line
    public :: spl_five_reg_line
    public :: spl_five_per_line

contains

    subroutine spl_build_line_inplace(order, periodic, n, h, work, line)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: order
        logical, intent(in) :: periodic
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(inout), contiguous :: work(:, :, 0:)
        integer, intent(in) :: line

        select case (order)
        case (3)
            if (periodic) then
                call spl_three_per_line(n, h, work(:, line, 0), work(:, line, 1), &
                                        work(:, line, 2), work(:, line, 3))
            else
                call spl_three_reg_line(n, h, work(:, line, 0), work(:, line, 1), &
                                        work(:, line, 2), work(:, line, 3))
            end if
        case (4)
            if (periodic) then
                call spl_four_per_line(n, h, work(:, line, 0), work(:, line, 1), &
                                       work(:, line, 2), work(:, line, 3), &
                                       work(:, line, 4))
            else
                call spl_four_reg_line(n, h, work(:, line, 0), work(:, line, 1), &
                                       work(:, line, 2), work(:, line, 3), &
                                       work(:, line, 4))
            end if
        case (5)
            if (periodic) then
                call spl_five_per_line(n, h, work(:, line, 0), work(:, line, 1), &
                                       work(:, line, 2), work(:, line, 3), &
                                       work(:, line, 4), work(:, line, 5))
            else
                call spl_five_reg_line(n, h, work(:, line, 0), work(:, line, 1), &
                                       work(:, line, 2), work(:, line, 3), &
                                       work(:, line, 4), work(:, line, 5))
            end if
        end select
    end subroutine spl_build_line_inplace

    subroutine spl_three_reg_line(n, h, a, b, c, d)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:)

        integer :: i, k, i5
        real(dp) :: c0, e, c1

        b(1) = 0.0d0
        d(1) = 0.0d0
        c0 = -4.0d0*h

        do i = 1, n - 2
            e = -3.0d0*((a(i + 2) - a(i + 1)) - (a(i + 1) - a(i)))/h
            c1 = c0 - b(i)*h
            b(i + 1) = h/c1
            d(i + 1) = (h*d(i) + e)/c1
        end do

        c(n) = 0.0d0
        k = n - 1
        do i = 1, k
            i5 = n - i
            c(i5) = b(i5)*c(i5 + 1) + d(i5)
        end do

        do i = 1, n - 1
            b(i) = (a(i + 1) - a(i))/h - h*(c(i + 1) + 2.0d0*c(i))/3.0d0
            d(i) = (c(i + 1) - c(i))/h/3.0d0
        end do
        b(n) = 0.0d0
        d(n) = 0.0d0
    end subroutine spl_three_reg_line

    subroutine tridiag_thomas(n, bdiag, rhs, x)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in), contiguous :: bdiag(:)
        real(dp), intent(inout), contiguous :: rhs(:)
        real(dp), intent(out), contiguous :: x(:)

        integer :: i
        real(dp) :: bet
        real(dp) :: gam(n)

        bet = bdiag(1)
        x(1) = rhs(1)/bet
        do i = 2, n
            gam(i) = 1.0d0/bet
            bet = bdiag(i) - gam(i)
            x(i) = (rhs(i) - x(i - 1))/bet
        end do
        do i = n - 1, 1, -1
            x(i) = x(i) - gam(i + 1)*x(i + 1)
        end do
    end subroutine tridiag_thomas

    subroutine spl_three_per_line(n, h, a, b, c, d)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:)

        integer :: nmx, i
        real(dp) :: psi, inv_h, beta, denom, fact
        real(dp) :: rhs(n - 1), y(n - 1), z(n - 1), u(n - 1)
        real(dp) :: bb(n - 1)

        nmx = n - 1
        inv_h = 1.0d0/h
        psi = 3.0d0*inv_h*inv_h

        rhs(1) = (a(2) - a(1) - a(n) + a(nmx))*psi
        do i = 3, nmx
            rhs(i - 1) = (a(i) - 2.0d0*a(i - 1) + a(i - 2))*psi
        end do
        rhs(nmx) = (a(n) - 2.0d0*a(nmx) + a(nmx - 1))*psi

        bb = 4.0d0
        beta = -bb(1)
        bb(1) = bb(1) - beta
        bb(nmx) = bb(nmx) - 1.0d0/beta

        u = 0.0d0
        u(1) = beta
        u(nmx) = 1.0d0

        call tridiag_thomas(nmx, bb, rhs, y)
        call tridiag_thomas(nmx, bb, u, z)

        denom = 1.0d0 + z(1) + z(nmx)/beta
        fact = (y(1) + y(nmx)/beta)/denom
        do i = 1, nmx
            c(i) = y(i) - fact*z(i)
        end do
        c(n) = c(1)

        do i = 1, nmx - 1
            b(i) = (a(i + 1) - a(i))*inv_h - h*(c(i + 1) + 2.0d0*c(i))/3.0d0
            d(i) = (c(i + 1) - c(i))*inv_h/3.0d0
        end do
        b(nmx) = (a(n) - a(nmx))*inv_h - h*(c(1) + 2.0d0*c(nmx))/3.0d0
        d(nmx) = (c(1) - c(nmx))*inv_h/3.0d0

        b(n) = b(1)
        d(n) = d(1)
    end subroutine spl_three_per_line

    subroutine spl_four_reg_line(n, h, a, b, c, d, e)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:), e(:)

        integer :: i, ip1
        real(dp) :: fac, fpl31, fpl40, fmn31, fmn40
        real(dp) :: alp(n), bet(n), gam(n)

        fpl31 = 0.5d0*(a(2) + a(4)) - a(3)
        fpl40 = 0.5d0*(a(1) + a(5)) - a(3)
        fmn31 = 0.5d0*(a(4) - a(2))
        fmn40 = 0.5d0*(a(5) - a(1))
        d(3) = (fmn40 - 2.0d0*fmn31)/6.0d0
        e(3) = (fpl40 - 4.0d0*fpl31)/12.0d0
        d(2) = d(3) - 4.0d0*e(3)
        d(1) = d(3) - 8.0d0*e(3)

        alp(1) = 0.0d0
        bet(1) = d(1) + d(2)

        do i = 1, n - 3
            ip1 = i + 1
            alp(ip1) = -1.0d0/(10.0d0 + alp(i))
            bet(ip1) = alp(ip1)*(bet(i) - 4.0d0*(a(i + 3) - 3.0d0*(a(i + 2) - a(ip1)) - a(i)))
        end do

        fpl31 = 0.5d0*(a(n - 3) + a(n - 1)) - a(n - 2)
        fpl40 = 0.5d0*(a(n - 4) + a(n)) - a(n - 2)
        fmn31 = 0.5d0*(a(n - 1) - a(n - 3))
        fmn40 = 0.5d0*(a(n) - a(n - 4))
        d(n - 2) = (fmn40 - 2.0d0*fmn31)/6.0d0
        e(n - 2) = (fpl40 - 4.0d0*fpl31)/12.0d0
        d(n - 1) = d(n - 2) + 4.0d0*e(n - 2)
        d(n) = d(n - 2) + 8.0d0*e(n - 2)

        gam(n - 1) = d(n) + d(n - 1)

        do i = n - 2, 1, -1
            gam(i) = gam(i + 1)*alp(i) + bet(i)
            d(i) = gam(i) - d(i + 1)
            e(i) = (d(i + 1) - d(i))/4.0d0
            c(i) = 0.5d0*(a(i + 2) + a(i)) - a(i + 1) - 0.125d0*(d(i + 2) + 12.0d0*d(i + 1) + &
                   11.0d0*d(i))
            b(i) = a(i + 1) - a(i) - c(i) - (3.0d0*d(i) + d(i + 1))/4.0d0
        end do

        b(n - 1) = b(n - 2) + 2.0d0*c(n - 2) + 3.0d0*d(n - 2) + 4.0d0*e(n - 2)
        c(n - 1) = c(n - 2) + 3.0d0*d(n - 2) + 6.0d0*e(n - 2)
        e(n - 1) = a(n) - a(n - 1) - b(n - 1) - c(n - 1) - d(n - 1)
        b(n) = b(n - 1) + 2.0d0*c(n - 1) + 3.0d0*d(n - 1) + 4.0d0*e(n - 1)
        c(n) = c(n - 1) + 3.0d0*d(n - 1) + 6.0d0*e(n - 1)
        e(n) = e(n - 1)

        fac = 1.0d0/h
        b = b*fac
        fac = fac/h
        c = c*fac
        fac = fac/h
        d = d*fac
        fac = fac/h
        e = e*fac
    end subroutine spl_four_reg_line

    subroutine spl_four_per_line(n, h, a, b, c, d, e)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:), e(:)

        integer :: i, ip1
        real(dp) :: fac, base1, base2, phi1, phi2, phi
        real(dp) :: alp(n), bet(n), gam(n)

        base1 = -5.0d0 + 2.0d0*sqrt(6.0d0)
        base2 = -5.0d0 - 2.0d0*sqrt(6.0d0)

        alp(1) = 0.0d0
        bet(1) = 0.0d0

        do i = 1, n - 3
            ip1 = i + 1
            alp(ip1) = -1.0d0/(10.0d0 + alp(i))
            bet(ip1) = alp(ip1)*(bet(i) - 4.0d0*(a(i + 3) - 3.0d0*(a(i + 2) - a(ip1)) - a(i)))
        end do
        alp(n - 1) = -1.0d0/(10.0d0 + alp(n - 2))
        bet(n - 1) = alp(n - 1)*(bet(n - 2) - 4.0d0*(a(2) - 3.0d0*(a(n) - a(n - 1)) - a(n - 2)))
        alp(n) = -1.0d0/(10.0d0 + alp(n - 1))
        bet(n) = alp(n)*(bet(n - 1) - 4.0d0*(a(3) - 3.0d0*(a(2) - a(n)) - a(n - 1)))

        gam(n) = bet(n)
        do i = n - 1, 1, -1
            gam(i) = gam(i + 1)*alp(i) + bet(i)
        end do

        phi1 = (gam(n)*base2 + gam(2))/(base2 - base1)/(1.0d0 - base1**(n - 1))
        phi2 = (gam(n)*base1 + gam(2))/(base2 - base1)/(1.0d0 - (1.0d0/base2)**(n - 1))

        do i = n, 1, -1
            gam(i) = gam(i) + phi2
            phi2 = phi2/base2
        end do

        do i = 1, n
            gam(i) = gam(i) + phi1
            phi1 = phi1*base1
        end do

        d(n) = 0.0d0
        do i = n - 1, 1, -1
            d(i) = gam(i) - d(i + 1)
        end do

        phi = -0.5d0*d(1)
        do i = 1, n
            d(i) = d(i) + phi
            phi = -phi
        end do

        e(n) = (d(2) - d(n))/4.0d0
        c(n) = 0.5d0*(a(3) + a(n)) - a(2) - 0.125d0*(d(3) + 12.0d0*d(2) + 11.0d0*d(n))
        b(n) = a(2) - a(n) - c(n) - (3.0d0*d(n) + d(2))/4.0d0
        e(n - 1) = (d(1) - d(n - 1))/4.0d0
        c(n - 1) = 0.5d0*(a(2) + a(n - 1)) - a(1) - 0.125d0*(d(2) + 12.0d0*d(1) + &
                   11.0d0*d(n - 1))
        b(n - 1) = a(1) - a(n - 1) - c(n - 1) - (3.0d0*d(n - 1) + d(1))/4.0d0

        do i = n - 2, 1, -1
            e(i) = (d(i + 1) - d(i))/4.0d0
            c(i) = 0.5d0*(a(i + 2) + a(i)) - a(i + 1) - 0.125d0*(d(i + 2) + 12.0d0*d(i + 1) + &
                   11.0d0*d(i))
            b(i) = a(i + 1) - a(i) - c(i) - (3.0d0*d(i) + d(i + 1))/4.0d0
        end do

        fac = 1.0d0/h
        b = b*fac
        fac = fac/h
        c = c*fac
        fac = fac/h
        d = d*fac
        fac = fac/h
        e = e*fac
    end subroutine spl_four_per_line

    subroutine spl_five_reg_line(n, h, a, b, c, d, e, f)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:), e(:), f(:)

        integer :: i, ip1
        real(dp) :: rhop, rhom, fac
        real(dp) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, det
        real(dp) :: b1, b2, b3
        real(dp) :: abeg, bbeg, cbeg, dbeg, ebeg, fbeg
        real(dp) :: aend, bend, cend, dend, eend, fend

        rhop = 13.0d0 + sqrt(105.0d0)
        rhom = 13.0d0 - sqrt(105.0d0)

        a11 = 1.0d0
        a12 = 1.0d0/4.0d0
        a13 = 1.0d0/16.0d0
        a21 = 3.0d0
        a22 = 27.0d0/4.0d0
        a23 = 9.0d0*27.0d0/16.0d0
        a31 = 5.0d0
        a32 = 125.0d0/4.0d0
        a33 = 5.0d0**5/16.0d0
        det = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a12*a21*a33 - &
              a13*a22*a31 - a11*a23*a32
        b1 = a(4) - a(3)
        b2 = a(5) - a(2)
        b3 = a(6) - a(1)
        bbeg = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a12*b2*a33 - &
               a13*a22*b3 - b1*a23*a32
        bbeg = bbeg/det
        dbeg = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - b1*a21*a33 - &
               a13*b2*a31 - a11*a23*b3
        dbeg = dbeg/det
        fbeg = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - a12*a21*b3 - &
               b1*a22*a31 - a11*b2*a32
        fbeg = fbeg/det
        b1 = a(n - 2) - a(n - 3)
        b2 = a(n - 1) - a(n - 4)
        b3 = a(n) - a(n - 5)
        bend = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a12*b2*a33 - &
               a13*a22*b3 - b1*a23*a32
        bend = bend/det
        dend = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - b1*a21*a33 - &
               a13*b2*a31 - a11*a23*b3
        dend = dend/det
        fend = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - a12*a21*b3 - &
               b1*a22*a31 - a11*b2*a32
        fend = fend/det

        a11 = 2.0d0
        a12 = 1.0d0/2.0d0
        a13 = 1.0d0/8.0d0
        a21 = 2.0d0
        a22 = 9.0d0/2.0d0
        a23 = 81.0d0/8.0d0
        a31 = 2.0d0
        a32 = 25.0d0/2.0d0
        a33 = 625.0d0/8.0d0
        det = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a12*a21*a33 - &
              a13*a22*a31 - a11*a23*a32
        b1 = a(4) + a(3)
        b2 = a(5) + a(2)
        b3 = a(6) + a(1)
        abeg = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a12*b2*a33 - &
               a13*a22*b3 - b1*a23*a32
        abeg = abeg/det
        cbeg = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - b1*a21*a33 - &
               a13*b2*a31 - a11*a23*b3
        cbeg = cbeg/det
        ebeg = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - a12*a21*b3 - &
               b1*a22*a31 - a11*b2*a32
        ebeg = ebeg/det
        b1 = a(n - 2) + a(n - 3)
        b2 = a(n - 1) + a(n - 4)
        b3 = a(n) + a(n - 5)
        aend = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a12*b2*a33 - &
               a13*a22*b3 - b1*a23*a32
        aend = aend/det
        cend = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - b1*a21*a33 - &
               a13*b2*a31 - a11*a23*b3
        cend = cend/det
        eend = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - a12*a21*b3 - &
               b1*a22*a31 - a11*b2*a32
        eend = eend/det

        b(1) = 0.0d0
        c(1) = ebeg*(2.0d0 + rhom) - 5.0d0*fbeg*(3.0d0 + 1.5d0*rhom)

        do i = 1, n - 4
            ip1 = i + 1
            b(ip1) = -1.0d0/(rhop + b(i))
            c(ip1) = b(ip1)*(c(i) - 5.0d0*(a(i + 4) - 4.0d0*a(i + 3) + &
                                          6.0d0*a(i + 2) - 4.0d0*a(ip1) + a(i)))
        end do

        e(n - 2) = eend*(2.0d0 + rhom) + 5.0d0*fend*(3.0d0 + 1.5d0*rhom)
        do i = n - 3, 1, -1
            e(i) = e(i + 1)*b(i) + c(i)
        end do

        b(1) = 0.0d0
        c(1) = ebeg - 2.5d0*5.0d0*fbeg

        do i = 1, n - 2
            ip1 = i + 1
            b(ip1) = -1.0d0/(rhom + b(i))
            c(ip1) = b(ip1)*(c(i) - e(i))
        end do

        e(n) = eend + 2.5d0*5.0d0*fend
        e(n - 1) = e(n)*b(n - 1) + c(n - 1)
        f(n - 1) = (e(n) - e(n - 1))/5.0d0
        e(n - 2) = e(n - 1)*b(n - 2) + c(n - 2)
        f(n - 2) = (e(n - 1) - e(n - 2))/5.0d0
        d(n - 2) = dend + 1.5d0*4.0d0*eend + 1.5d0**2*10.0d0*fend

        do i = n - 3, 1, -1
            e(i) = e(i + 1)*b(i) + c(i)
            f(i) = (e(i + 1) - e(i))/5.0d0
            d(i) = (a(i + 3) - 3.0d0*a(i + 2) + 3.0d0*a(i + 1) - a(i))/6.0d0 - &
                   (e(i + 3) + 27.0d0*e(i + 2) + 93.0d0*e(i + 1) + &
                    59.0d0*e(i))/30.0d0
            c(i) = 0.5d0*(a(i + 2) + a(i)) - a(i + 1) - 0.5d0*d(i + 1) - &
                   2.5d0*d(i) - 0.1d0*(e(i + 2) + 18.0d0*e(i + 1) + &
                                       31.0d0*e(i))
            b(i) = a(i + 1) - a(i) - c(i) - d(i) - 0.2d0*(4.0d0*e(i) + e(i + 1))
        end do

        do i = n - 3, n
            b(i) = b(i - 1) + 2.0d0*c(i - 1) + 3.0d0*d(i - 1) + 4.0d0*e(i - 1) + &
                   5.0d0*f(i - 1)
            c(i) = c(i - 1) + 3.0d0*d(i - 1) + 6.0d0*e(i - 1) + 10.0d0*f(i - 1)
            d(i) = d(i - 1) + 4.0d0*e(i - 1) + 10.0d0*f(i - 1)
            if (i /= n) f(i) = a(i + 1) - a(i) - b(i) - c(i) - d(i) - e(i)
        end do
        f(n) = f(n - 1)

        fac = 1.0d0/h
        b = b*fac
        fac = fac/h
        c = c*fac
        fac = fac/h
        d = d*fac
        fac = fac/h
        e = e*fac
        fac = fac/h
        f = f*fac
    end subroutine spl_five_reg_line

    subroutine spl_five_per_line(n, h, a, b, c, d, e, f)
#ifdef _OPENACC
        !$acc routine seq
#endif
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in), contiguous :: a(:)
        real(dp), intent(out), contiguous :: b(:), c(:), d(:), e(:), f(:)

        integer :: i, ip1
        real(dp) :: rhop, rhom, fac, xplu, xmin, gammao_p, gammao_m_redef
        real(dp) :: dummy

        rhop = 13.0d0 + sqrt(105.0d0)
        rhom = 13.0d0 - sqrt(105.0d0)

        ! Stage 1: solve for gamma (stored in e) using rhop.
        b(1) = 0.0d0
        c(1) = 0.0d0

        do i = 1, n - 4
            ip1 = i + 1
            b(ip1) = -1.0d0/(rhop + b(i))
            c(ip1) = b(ip1)*(c(i) - 5.0d0*(a(i + 4) - 4.0d0*a(i + 3) + &
                                          6.0d0*a(i + 2) - 4.0d0*a(ip1) + a(i)))
        end do
        b(n - 2) = -1.0d0/(rhop + b(n - 3))
        c(n - 2) = b(n - 2)*(c(n - 3) - 5.0d0*(a(2) - 4.0d0*a(1) + &
                                              6.0d0*a(n - 1) - 4.0d0*a(n - 2) + &
                                              a(n - 3)))
        b(n - 1) = -1.0d0/(rhop + b(n - 2))
        c(n - 1) = b(n - 1)*(c(n - 2) - 5.0d0*(a(3) - 4.0d0*a(2) + &
                                              6.0d0*a(1) - 4.0d0*a(n - 1) + &
                                              a(n - 2)))
        b(n) = -1.0d0/(rhop + b(n - 1))
        c(n) = b(n)*(c(n - 1) - 5.0d0*(a(4) - 4.0d0*a(3) + 6.0d0*a(2) - &
                                       4.0d0*a(1) + a(n - 1)))

        e(n) = c(n)
        do i = n - 1, 1, -1
            e(i) = e(i + 1)*b(i) + c(i)
        end do

        xplu = sqrt(0.25d0*rhop**2 - 1.0d0) - 0.5d0*rhop
        xmin = -sqrt(0.25d0*rhop**2 - 1.0d0) - 0.5d0*rhop
        dummy = (1.0d0/xmin)**(n - 1)
        gammao_m_redef = (e(2) + xplu*e(n))/(1.0d0 - dummy)/(xmin - xplu)
        gammao_p = (e(2) + xmin*e(n))/(xplu**(n - 1) - 1.0d0)/(xplu - xmin)
        e(1) = e(1) + gammao_m_redef*dummy + gammao_p
        do i = 2, n
            e(i) = e(i) + gammao_m_redef*(1.0d0/xmin)**(n - i) + &
                   gammao_p*xplu**(i - 1)
        end do

        ! Stage 2: solve for e using rhom.
        b(1) = 0.0d0
        c(1) = 0.0d0

        do i = 1, n - 1
            ip1 = i + 1
            b(ip1) = -1.0d0/(rhom + b(i))
            c(ip1) = b(ip1)*(c(i) - e(i))
        end do

        e(n) = c(n)
        do i = n - 1, 1, -1
            e(i) = e(i + 1)*b(i) + c(i)
        end do

        xplu = sqrt(0.25d0*rhom**2 - 1.0d0) - 0.5d0*rhom
        xmin = -sqrt(0.25d0*rhom**2 - 1.0d0) - 0.5d0*rhom
        dummy = (1.0d0/xmin)**(n - 1)
        gammao_m_redef = (e(2) + xplu*e(n))/(1.0d0 - dummy)/(xmin - xplu)
        gammao_p = (e(2) + xmin*e(n))/(xplu**(n - 1) - 1.0d0)/(xplu - xmin)
        e(1) = e(1) + gammao_m_redef*dummy + gammao_p
        do i = 2, n
            e(i) = e(i) + gammao_m_redef*(1.0d0/xmin)**(n - i) + &
                   gammao_p*xplu**(i - 1)
        end do

        do i = n - 1, 1, -1
            f(i) = (e(i + 1) - e(i))/5.0d0
        end do
        f(n) = f(1)

        d(n - 1) = (a(3) - 3.0d0*a(2) + 3.0d0*a(1) - a(n - 1))/6.0d0 - &
                   (e(3) + 27.0d0*e(2) + 93.0d0*e(1) + 59.0d0*e(n - 1))/30.0d0
        d(n - 2) = (a(2) - 3.0d0*a(1) + 3.0d0*a(n - 1) - a(n - 2))/6.0d0 - &
                   (e(2) + 27.0d0*e(1) + 93.0d0*e(n - 1) + 59.0d0*e(n - 2))/30.0d0
        do i = n - 3, 1, -1
            d(i) = (a(i + 3) - 3.0d0*a(i + 2) + 3.0d0*a(i + 1) - a(i))/6.0d0 - &
                   (e(i + 3) + 27.0d0*e(i + 2) + 93.0d0*e(i + 1) + &
                    59.0d0*e(i))/30.0d0
        end do
        d(n) = d(1)

        c(n - 1) = 0.5d0*(a(2) + a(n - 1)) - a(1) - 0.5d0*d(1) - 2.5d0*d(n - 1) - &
                   0.1d0*(e(2) + 18.0d0*e(1) + 31.0d0*e(n - 1))
        b(n - 1) = a(1) - a(n - 1) - c(n - 1) - d(n - 1) - 0.2d0*(4.0d0*e(n - 1) + &
                                                                 e(1))
        do i = n - 2, 1, -1
            c(i) = 0.5d0*(a(i + 2) + a(i)) - a(i + 1) - 0.5d0*d(i + 1) - &
                   2.5d0*d(i) - 0.1d0*(e(i + 2) + 18.0d0*e(i + 1) + 31.0d0*e(i))
            b(i) = a(i + 1) - a(i) - c(i) - d(i) - 0.2d0*(4.0d0*e(i) + e(i + 1))
        end do
        b(n) = b(1)
        c(n) = c(1)

        fac = 1.0d0/h
        b = b*fac
        fac = fac/h
        c = c*fac
        fac = fac/h
        d = d*fac
        fac = fac/h
        e = e*fac
        fac = fac/h
        f = f*fac
    end subroutine spl_five_per_line

end module spline_build_lines
