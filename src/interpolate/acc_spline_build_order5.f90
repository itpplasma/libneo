module acc_spline_build_order5
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: spl_five_reg_line
    public :: spl_five_per_line

contains

    subroutine spl_five_reg_line(n, h, a, b, c, d, e, f)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in) :: a(:)
        real(dp), intent(out) :: b(:), c(:), d(:), e(:), f(:)

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
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(in) :: a(:)
        real(dp), intent(out) :: b(:), c(:), d(:), e(:), f(:)

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

end module acc_spline_build_order5
