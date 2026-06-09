module neo_dawson

    ! Dawson integral F(x) = exp(-x^2) * integral_0^x exp(t^2) dt.
    ! Maclaurin series (DLMF 7.6.4) for |x| < 1, the sampling-theorem method of
    ! Rybicki (1989), Computers in Physics 3, 85 for 1 <= |x| < 10, and the
    ! asymptotic expansion (DLMF 7.12.2) for |x| >= 10.

    use libneo_kinds, only: dp

    implicit none
    private

    public :: dawson

    real(dp), parameter :: inv_sqrt_pi = 0.564189583547756287d0

    ! Maclaurin coefficients (-2)^n / (2n+1)!!, n = 0..20
    integer, parameter :: nmac = 20
    real(dp), parameter :: mac(0:nmac) = [ &
        1.0d0, &
        -0.666666666666666667d0, &
        0.266666666666666667d0, &
        -0.0761904761904761905d0, &
        0.0169312169312169312d0, &
        -0.0030784030784030784d0, &
        0.0004736004736004736d0, &
        -0.0000631467298133964801d0, &
        7.42902703687017413d-6, &
        -7.82002845986334118d-7, &
        7.44764615225080113d-8, &
        -6.47621404543547924d-9, &
        5.18097123634838339d-10, &
        -3.83775647136917288d-11, &
        2.64672860094425716d-12, &
        -1.70756683931887559d-13, &
        1.03488899352659127d-14, &
        -5.91365139158052152d-16, &
        3.19656831977325487d-17, &
        -1.63926580501192558d-18, &
        7.9964185610337833d-20]

    ! Rybicki sampling step h and weights exp(-((2k-1)*h)^2), k = 1..14;
    ! sampling error ~ exp(-(pi/(2h))^2) ~ 7e-18 for h = 1/4
    real(dp), parameter :: h_ryb = 0.25d0
    integer, parameter :: nryb = 14
    real(dp), parameter :: cryb(nryb) = [ &
        0.939413062813475786d0, &
        0.56978282473092301d0, &
        0.209611387151097823d0, &
        0.0467706223839589837d0, &
        0.00632971542748574658d0, &
        0.000519574682154838482d0, &
        0.0000258681002226541213d0, &
        7.8114894083044908d-7, &
        1.43072419185676883d-8, &
        1.58939100945163665d-10, &
        1.07092323825080765d-12, &
        4.37661850287084989d-15, &
        1.0848552640429378d-17, &
        1.63101392267018568d-20]

    ! asymptotic coefficients (2m-1)!!, m = 0..14 (exact in double)
    integer, parameter :: nasy = 14
    real(dp), parameter :: asy(0:nasy) = [ &
        1.0d0, 1.0d0, 3.0d0, 15.0d0, 105.0d0, 945.0d0, 10395.0d0, &
        135135.0d0, 2027025.0d0, 34459425.0d0, 654729075.0d0, &
        13749310575.0d0, 316234143225.0d0, 7905853580625.0d0, &
        213458046676875.0d0]

contains

    elemental function dawson(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f

        real(dp) :: a

        a = abs(x)
        if (a < 1.0d0) then
            f = dawson_series(x)
        else if (a < 10.0d0) then
            f = sign(dawson_rybicki(a), x)
        else
            f = sign(dawson_asymptotic(a), x)
        end if
    end function dawson

    pure function dawson_series(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f

        real(dp) :: t, s
        integer :: n

        t = x*x
        s = mac(nmac)
        do n = nmac - 1, 0, -1
            s = mac(n) + t*s
        end do
        f = x*s
    end function dawson_series

    pure function dawson_rybicki(a) result(f)
        real(dp), intent(in) :: a
        real(dp) :: f

        integer :: n0, k
        real(dp) :: xp, gauss, e1, e2, einv2, ep, em, s

        n0 = 2*nint(0.5d0*a/h_ryb)
        xp = a - real(n0, dp)*h_ryb
        gauss = exp(-xp*xp)
        e1 = exp(2.0d0*h_ryb*xp)
        e2 = e1*e1
        einv2 = 1.0d0/e2
        ep = e1
        em = 1.0d0/e1
        s = 0.0d0
        do k = 1, nryb
            s = s + cryb(k)*(ep/real(n0 + 2*k - 1, dp) + em/real(n0 - 2*k + 1, dp))
            ep = ep*e2
            em = em*einv2
        end do
        f = s*gauss*inv_sqrt_pi
    end function dawson_rybicki

    pure function dawson_asymptotic(a) result(f)
        real(dp), intent(in) :: a
        real(dp) :: f

        real(dp) :: t, s
        integer :: m

        ! t underflows to 0 for a > ~1e154, leaving the exact 1/(2a) tail
        t = 0.5d0/(a*a)
        s = asy(nasy)
        do m = nasy - 1, 0, -1
            s = asy(m) + t*s
        end do
        f = 0.5d0/a*s
    end function dawson_asymptotic

end module neo_dawson
