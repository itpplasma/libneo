module neo_bessel_k
    ! Modified Bessel functions of the second kind K_n(x), integer n >= 0.
    ! K_0, K_1 for x <= 3: power series DLMF 10.31.1/10.31.2 with fixed
    ! polynomial coefficients in q = x^2/4 (generated with mpmath at 50
    ! digits). For x > 3: Chebyshev fits of sqrt(2x/pi) e^x K_nu(x) in
    ! s = 1/x on (1/16, 1/3) and [0, 1/16] (mpmath, 50 digits), matching
    ! the asymptotic form DLMF 10.40.2. K_n: upward recurrence DLMF 10.29.1,
    ! K_{n+1} = K_{n-1} + (2n/x) K_n, forward-stable for K.
    use libneo_kinds, only: dp

    implicit none
    private

    public :: bessel_kn

    real(dp), parameter :: euler_gamma = 0.57721566490153286060651209008240243d0
    real(dp), parameter :: half_pi = 1.5707963267948966192313216916397514d0
    real(dp), parameter :: series_x_max = 3.0d0
    real(dp), parameter :: cheb_hi_x_min = 16.0d0
    real(dp), parameter :: cheb_mid_s_lo = 0.0625d0
    real(dp), parameter :: cheb_mid_s_hi = 0.33333333333333333333333333333333333d0

    ! I_0(x) = sum_k q^k/(k!)^2: coefficients 1/(k!)^2
    real(dp), parameter :: poly_i0(18) = [ &
        1.0d0, 1.0d0, &
        0.25d0, 0.027777777777777777778d0, &
        0.0017361111111111111111d0, 0.000069444444444444444444d0, &
        1.9290123456790123457d-6, 3.9367598891408415218d-8, &
        6.1511873267825648778d-10, 7.5940584281266233059d-12, &
        7.5940584281266233059d-14, 6.2760813455591928148d-16, &
        4.3583898233049950103d-18, 2.5789288895295828463d-20, &
        1.3157800456783585951d-22, 5.8479113141260382003d-25, &
        2.284340357080483672d-27, 7.9042918930120542283d-30]
    ! coefficients h_k/(k!)^2 with h_k = sum_{j=1..k} 1/j
    real(dp), parameter :: poly_s0(18) = [ &
        0.0d0, 1.0d0, &
        0.375d0, 0.050925925925925925926d0, &
        0.0036168981481481481481d0, 0.00015856481481481481481d0, &
        4.7260802469135802469d-6, 1.0207455998272324803d-7, &
        1.6718048413148328114d-9, 2.1483350211950276805d-11, &
        2.2242756054762939135d-13, 1.8952995870061529211d-15, &
        1.3525001839484811536d-17, 8.2013388136826374592d-20, &
        4.2783408265702079911d-22, 1.9404708872364882507d-24, &
        7.7227356755850624588d-27, 2.7187227120298503047d-29]
    ! I_1(x) = (x/2) sum_k q^k/(k!(k+1)!): coefficients 1/(k!(k+1)!)
    real(dp), parameter :: poly_i1(18) = [ &
        1.0d0, 0.5d0, &
        0.083333333333333333333d0, 0.0069444444444444444444d0, &
        0.00034722222222222222222d0, 0.000011574074074074074074d0, &
        2.7557319223985890653d-7, 4.9209498614260519022d-9, &
        6.8346525853139609753d-11, 7.5940584281266233059d-13, &
        6.9036894801151120963d-15, 5.2300677879659940123d-17, &
        3.3526075563884577002d-19, 1.8420920639497020331d-21, &
        8.7718669711890573004d-24, 3.6549445713287738752d-26, &
        1.3437296218120492188d-28, 4.3912732738955856824d-31]
    ! coefficients (h_k + h_{k+1} - 2 gamma)/(k!(k+1)!)
    real(dp), parameter :: poly_s1(18) = [ &
        -0.15443132980306572121d0, 0.67278433509846713939d0, &
        0.18157516696085563434d0, 0.019182189839330562121d0, &
        0.0011153594919665281061d0, 0.000041422476892711430696d0, &
        1.0715459140911808686d-6, 2.0452860035938779413d-8, &
        3.0020487465891878859d-10, 3.4959287296928819208d-12, &
        3.309914735250272068d-14, 2.5986411321011287351d-16, &
        1.7195232826992565241d-18, 9.7212075188236180165d-21, &
        4.7502817433276669896d-23, 2.0264937604328579082d-25, &
        7.6133707277671159281d-28, 2.5382566314880592899d-30]

    real(dp), parameter :: cheb_mid_k0(18) = [ &
        1.9559376047895837408d0, -0.014059053153789941523d0, &
        0.00039802254448432713646d0, -0.000018657916479681936657d0, &
        1.1670921100603393265d-6, -8.8637049995049946834d-8, &
        7.7618435419763736638d-9, -7.5903717486808993347d-10, &
        8.1130444797671898315d-11, -9.3346884671836175433d-12, &
        1.1431343028648823935d-12, -1.4771049122256338483d-13, &
        2.0002339096669601218d-14, -2.8230963387131460904d-15, &
        4.1342702911603966579d-16, -6.2586795525484043848d-17, &
        9.7637329920200614695d-18, -1.5654516590192795234d-18]
    real(dp), parameter :: cheb_mid_k1(18) = [ &
        2.1390906750838713041d0, 0.045796323598856923076d0, &
        -0.00071639260360370820452d0, 0.000028086590879043600368d0, &
        -1.6076666325349922967d-6, 1.1571456080375915070d-7, &
        -9.7722547348374030052d-9, 9.3086831666491458465d-10, &
        -9.7528215763417439059d-11, 1.1045807616711904736d-11, &
        -1.3354794780368619141d-12, 1.7074202453340075633d-13, &
        -2.2914753368530012268d-14, 3.2093797833233986927d-15, &
        -4.6686970748489713955d-16, 7.0264558327001338014d-17, &
        -1.0904794030447221775d-17, 1.7403276745873914004d-18]

    real(dp), parameter :: cheb_hi_k0(12) = [ &
        1.9923831602004083696d0, -0.0037766317031612104235d0, &
        0.000031310052233667130168d0, -4.6784426011890146d-7, &
        1.0013945577406131153d-8, -2.7653023813480403757d-10, &
        9.2831416580801534697d-12, -3.6467062746741891765d-13, &
        1.6325975548816327494d-14, -8.1701412478988243465d-16, &
        4.5035599694971361505d-17, -2.7029444926470642288d-18]
    real(dp), parameter :: cheb_hi_k1(12) = [ &
        2.0231087348982674842d0, 0.011500735467006774991d0, &
        -0.000052954102844829022231d0, 6.6446695239227473635d-7, &
        -1.3058062241211819979d-8, 3.4269978500370871047d-10, &
        -1.1121538753031974548d-11, 4.2645352569241138993d-13, &
        -1.8748591483665248044d-14, 9.2508042350153445531d-16, &
        -5.0417867748703588125d-17, 2.9979881728938231118d-18]

    contains

    pure function bessel_kn(n, x) result(kn)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: kn

        real(dp) :: k_prev, k_curr, k_next
        integer :: n_abs, m

        n_abs = abs(n)
        call bessel_k01(x, k_prev, k_curr)
        if (n_abs == 0) then
            kn = k_prev
        else if (n_abs == 1) then
            kn = k_curr
        else
            do m = 1, n_abs - 1
                k_next = k_prev + (2.0d0*m/x)*k_curr
                k_prev = k_curr
                k_curr = k_next
            end do
            kn = k_curr
        end if
    end function bessel_kn


    pure subroutine bessel_k01(x, k0, k1)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: k0, k1

        real(dp) :: u, pref

        if (x <= series_x_max) then
            call k01_series(x, k0, k1)
        else
            pref = sqrt(half_pi/x)*exp(-x)
            if (x < cheb_hi_x_min) then
                u = (2.0d0/x - cheb_mid_s_lo - cheb_mid_s_hi) &
                    /(cheb_mid_s_hi - cheb_mid_s_lo)
                k0 = pref*cheb_eval(cheb_mid_k0, u)
                k1 = pref*cheb_eval(cheb_mid_k1, u)
            else
                u = 32.0d0/x - 1.0d0
                k0 = pref*cheb_eval(cheb_hi_k0, u)
                k1 = pref*cheb_eval(cheb_hi_k1, u)
            end if
        end if
    end subroutine bessel_k01


    ! DLMF 10.31.2: K_0 = -(ln(x/2) + gamma) I_0 + sum_{k>=1} h_k q^k/(k!)^2
    ! DLMF 10.31.1 (n=1): K_1 = 1/x + ln(x/2) I_1
    !                           - (x/4) sum_{k>=0} (h_k + h_{k+1} - 2 gamma) t_k,
    ! with q = x^2/4, h_k = sum_{j=1..k} 1/j, t_k = q^k/(k! (k+1)!).
    pure subroutine k01_series(x, k0, k1)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: k0, k1

        real(dp) :: q, log_half_x, i0, s0, i1s, s1

        q = 0.25d0*x*x
        log_half_x = log(0.5d0*x)
        i0 = poly_eval(poly_i0, q)
        s0 = poly_eval(poly_s0, q)
        i1s = poly_eval(poly_i1, q)
        s1 = poly_eval(poly_s1, q)
        k0 = -(log_half_x + euler_gamma)*i0 + s0
        k1 = 1.0d0/x + 0.5d0*x*(log_half_x*i1s - 0.5d0*s1)
    end subroutine k01_series


    pure function poly_eval(c, q) result(p)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: q
        real(dp) :: p

        integer :: k

        p = c(size(c))
        do k = size(c) - 1, 1, -1
            p = p*q + c(k)
        end do
    end function poly_eval


    ! Clenshaw recurrence; c(1) is the halved T_0 coefficient convention.
    pure function cheb_eval(c, u) result(f)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: u
        real(dp) :: f

        real(dp) :: b0, b1, b2, two_u
        integer :: j

        two_u = 2.0d0*u
        b1 = 0.0d0
        b2 = 0.0d0
        do j = size(c), 2, -1
            b0 = two_u*b1 - b2 + c(j)
            b2 = b1
            b1 = b0
        end do
        f = u*b1 - b2 + 0.5d0*c(1)
    end function cheb_eval

end module neo_bessel_k
