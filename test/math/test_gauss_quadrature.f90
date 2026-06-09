program test_gauss_quadrature

    ! Reference nodes and weights from mpmath at 40 significant digits:
    ! Golub-Welsch eigensolve of the Jacobi matrix with mp.eigsy,
    ! cross-checked against roots of mpmath legendre/laguerre polynomials.

    use libneo_kinds, only: dp
    use neo_gauss_quadrature, only: gauss_legendre, gauss_legendre_ab, &
                                    gauss_gen_laguerre
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: alpha_list(4) = [0.0_dp, 0.5_dp, 2.5_dp, 3.5_dp]

    call test_legendre_full_small
    call test_legendre_sampled_large
    call test_legendre_polynomial_exactness
    call test_legendre_ab
    call test_gen_laguerre_full_small
    call test_gen_laguerre_sampled_large
    call test_gen_laguerre_moments

contains

    subroutine test_legendre_full_small
        real(dp), parameter :: x_ref_4(4) = [ &
            -0.86113631159405258d0, &
            -0.33998104358485626d0, &
            0.33998104358485626d0, &
            0.86113631159405258d0]
        real(dp), parameter :: w_ref_4(4) = [ &
            0.34785484513745386d0, &
            0.65214515486254614d0, &
            0.65214515486254614d0, &
            0.34785484513745386d0]
        real(dp), parameter :: x_ref_8(8) = [ &
            -0.96028985649753623d0, &
            -0.79666647741362674d0, &
            -0.52553240991632899d0, &
            -0.18343464249564980d0, &
            0.18343464249564980d0, &
            0.52553240991632899d0, &
            0.79666647741362674d0, &
            0.96028985649753623d0]
        real(dp), parameter :: w_ref_8(8) = [ &
            0.10122853629037626d0, &
            0.22238103445337447d0, &
            0.31370664587788729d0, &
            0.36268378337836198d0, &
            0.36268378337836198d0, &
            0.31370664587788729d0, &
            0.22238103445337447d0, &
            0.10122853629037626d0]
        real(dp), parameter :: abs_tol_x = 1.0d-14, rel_tol_w = 1.0d-13

        real(dp) :: x4(4), w4(4), x8(8), w8(8)
        integer :: i

        call print_test("test_legendre_full_small")
        call gauss_legendre(4, x4, w4)
        call gauss_legendre(8, x8, w8)
        do i = 1, 4
            if (abs(x4(i) - x_ref_4(i)) > abs_tol_x .or. &
                abs(w4(i)/w_ref_4(i) - 1.0_dp) > rel_tol_w) then
                call print_fail
                print *, "n=4 i=", i, x4(i), w4(i)
                stop 1
            end if
        end do
        do i = 1, 8
            if (abs(x8(i) - x_ref_8(i)) > abs_tol_x .or. &
                abs(w8(i)/w_ref_8(i) - 1.0_dp) > rel_tol_w) then
                call print_fail
                print *, "n=8 i=", i, x8(i), w8(i)
                stop 1
            end if
        end do
        call print_ok
    end subroutine test_legendre_full_small

    subroutine test_legendre_sampled_large
        ! Indices 1, n/2, n for n = 16, 32, 64
        integer, parameter :: n_list(3) = [16, 32, 64]
        real(dp), parameter :: x_ref(9) = [ &
            -0.98940093499164993d0, &
            -0.095012509837637440d0, &
            0.98940093499164993d0, &
            -0.99726386184948156d0, &
            -0.048307665687738316d0, &
            0.99726386184948156d0, &
            -0.99930504173577214d0, &
            -0.024350292663424433d0, &
            0.99930504173577214d0]
        real(dp), parameter :: w_ref(9) = [ &
            0.027152459411754095d0, &
            0.18945061045506850d0, &
            0.027152459411754095d0, &
            0.0070186100094700966d0, &
            0.096540088514727801d0, &
            0.0070186100094700966d0, &
            0.0017832807216964329d0, &
            0.048690957009139720d0, &
            0.0017832807216964329d0]
        real(dp), parameter :: abs_tol_x = 1.0d-14, rel_tol_w = 1.0d-13

        real(dp) :: x(64), w(64)
        integer :: j, k, m, n, idx(3)

        call print_test("test_legendre_sampled_large")
        m = 0
        do j = 1, size(n_list)
            n = n_list(j)
            idx = [1, n/2, n]
            call gauss_legendre(n, x(1:n), w(1:n))
            do k = 1, 3
                m = m + 1
                if (abs(x(idx(k)) - x_ref(m)) > abs_tol_x .or. &
                    abs(w(idx(k))/w_ref(m) - 1.0_dp) > rel_tol_w) then
                    call print_fail
                    print *, "n=", n, " i=", idx(k), x(idx(k)), w(idx(k))
                    stop 1
                end if
            end do
        end do
        call print_ok
    end subroutine test_legendre_sampled_large

    subroutine test_legendre_polynomial_exactness
        ! n-point rule must integrate x**k over [-1,1] exactly for k <= 2n-1
        integer, parameter :: n_list(6) = [1, 4, 8, 16, 32, 64]
        real(dp), parameter :: abs_tol = 1.0d-13

        real(dp) :: x(64), w(64), s, ref
        integer :: j, k, n

        call print_test("test_legendre_polynomial_exactness")
        do j = 1, size(n_list)
            n = n_list(j)
            call gauss_legendre(n, x(1:n), w(1:n))
            do k = 0, 2*n - 1
                s = sum(w(1:n)*x(1:n)**k)
                ref = 0.0_dp
                if (mod(k, 2) == 0) ref = 2.0_dp/real(k + 1, dp)
                if (abs(s - ref) > abs_tol) then
                    call print_fail
                    print *, "n=", n, " k=", k, " moment=", s, " ref=", ref
                    stop 1
                end if
            end do
        end do
        call print_ok
    end subroutine test_legendre_polynomial_exactness

    subroutine test_legendre_ab
        integer, parameter :: n = 32
        real(dp), parameter :: pi = 3.14159265358979324_dp
        real(dp), parameter :: rel_tol = 1.0d-14

        real(dp) :: x(n), w(n), s

        call print_test("test_legendre_ab")
        call gauss_legendre_ab(n, 0.0_dp, pi, x, w)
        if (abs(sum(w)/pi - 1.0_dp) > rel_tol) then
            call print_fail
            print *, "weight sum on [0,pi]:", sum(w)
            stop 1
        end if
        s = sum(w*sin(x))
        if (abs(s/2.0_dp - 1.0_dp) > rel_tol) then
            call print_fail
            print *, "integral of sin on [0,pi]:", s
            stop 1
        end if
        call gauss_legendre_ab(n, -2.0_dp, 3.0_dp, x, w)
        s = sum(w*exp(x))
        if (abs(s/(exp(3.0_dp) - exp(-2.0_dp)) - 1.0_dp) > rel_tol) then
            call print_fail
            print *, "integral of exp on [-2,3]:", s
            stop 1
        end if
        call print_ok
    end subroutine test_legendre_ab

    subroutine test_gen_laguerre_full_small
        ! n=8, alpha = 0, 0.5, 2.5, 3.5 flattened in order
        real(dp), parameter :: x_ref(32) = [ &
            0.17027963230510100d0, &
            0.90370177679937991d0, &
            2.2510866298661307d0, &
            4.2667001702876588d0, &
            7.0459054023934657d0, &
            10.758516010180995d0, &
            15.740678641278005d0, &
            22.863131736889264d0, &
            0.28263364811659913d0, &
            1.1398738015816137d0, &
            2.6015248434060294d0, &
            4.7241145375277907d0, &
            7.6052562992316137d0, &
            11.417182076545830d0, &
            16.499410797655816d0, &
            23.730003995934708d0, &
            0.86004162376038878d0, &
            2.1662435289922482d0, &
            4.0437426773491275d0, &
            6.5589225384279037d0, &
            9.8172555052257052d0, &
            13.999479750715056d0, &
            19.457447431814625d0, &
            27.096866943714945d0, &
            1.2083204369489536d0, &
            2.7193482009664990d0, &
            4.7857319321828601d0, &
            7.4801580091993414d0, &
            10.911908001846749d0, &
            15.265599663119231d0, &
            20.898773393105829d0, &
            28.730160362630536d0]
        real(dp), parameter :: w_ref(32) = [ &
            0.36918858934163753d0, &
            0.41878678081434296d0, &
            0.17579498663717181d0, &
            0.033343492261215652d0, &
            0.0027945362352256725d0, &
            9.0765087733582131d-5, &
            8.4857467162725315d-7, &
            1.0480011748715104d-9, &
            0.22713936195247155d0, &
            0.39359454280361523d0, &
            0.21290897086722816d0, &
            0.047877483203138173d0, &
            0.0045425174747626382d0, &
            0.00016240460018532570d0, &
            1.6423774138060987d-6, &
            2.1739431266309180d-9, &
            0.30118672531843548d0, &
            1.2529637803437606d0, &
            1.2583457960796291d0, &
            0.44734423636842068d0, &
            0.060582550279059224d0, &
            0.0028904817699157067d0, &
            3.7338528665686646d-5, &
            6.1759956053533482d-8, &
            0.72703666811450358d0, &
            3.8870095950840757d0, &
            4.7364075115745778d0, &
            1.9634439889879067d0, &
            0.30159594888610425d0, &
            0.016007214602785849d0, &
            0.00022705927903432287d0, &
            4.1003846070191387d-7]
        real(dp), parameter :: rel_tol_x = 1.0d-12, rel_tol_w = 1.0d-12

        real(dp) :: x(8), w(8)
        integer :: ia, i, m

        call print_test("test_gen_laguerre_full_small")
        m = 0
        do ia = 1, size(alpha_list)
            call gauss_gen_laguerre(8, alpha_list(ia), x, w)
            do i = 1, 8
                m = m + 1
                if (abs(x(i)/x_ref(m) - 1.0_dp) > rel_tol_x .or. &
                    abs(w(i)/w_ref(m) - 1.0_dp) > rel_tol_w) then
                    call print_fail
                    print *, "alpha=", alpha_list(ia), " i=", i, x(i), w(i)
                    stop 1
                end if
            end do
        end do
        call print_ok
    end subroutine test_gen_laguerre_full_small

    subroutine test_gen_laguerre_sampled_large
        ! Indices 1, n/2, n for n = 16, 32, 64 and alpha = 0, 0.5, 2.5, 3.5;
        ! the smallest weights reach 2e-101 and must stay relatively accurate
        integer, parameter :: n_list(3) = [16, 32, 64]
        real(dp), parameter :: x_ref(36) = [ &
            0.087649410478927840d0, &
            9.4383143363919388d0, &
            51.701160339543318d0, &
            0.14739918461631114d0, &
            9.9212191360724295d0, &
            52.618366255753246d0, &
            0.46921149778855161d0, &
            11.850358073150744d0, &
            56.225744200404799d0, &
            0.67193070292070396d0, &
            12.814186222858152d0, &
            57.996108137129266d0, &
            0.044489365833267018d0, &
            19.855860940336055d0, &
            111.75139809793770d0, &
            0.075352743443543200d0, &
            20.352117419061982d0, &
            112.69994395176783d0, &
            0.24625231492422103d0, &
            22.334331508962957d0, &
            116.45992000951651d0, &
            0.35688470093352815d0, &
            23.324051960829249d0, &
            118.32043834546568d0, &
            0.022415874146705280d0, &
            40.730835444458626d0, &
            234.80957917132616d0, &
            0.038108165657215796d0, &
            41.233926467273473d0, &
            235.77755051198645d0, &
            0.12632875645146994d0, &
            43.244443275963830d0, &
            239.63130475601646d0, &
            0.18433108596198790d0, &
            44.248686994140381d0, &
            241.54759237297520d0]
        real(dp), parameter :: w_ref(36) = [ &
            0.20615171495780099d0, &
            0.00020427191530827846d0, &
            4.1614623703728552d-22, &
            0.097740989137131689d0, &
            0.00040525569008016565d0, &
            1.2137123039229541d-21, &
            0.053120314720276348d0, &
            0.0097424848196979495d0, &
            1.1011483894446578d-19, &
            0.088028076953009481d0, &
            0.059899040621830508d0, &
            1.1881425642969675d-18, &
            0.10921834195238497d0, &
            6.3506022266258067d-9, &
            4.5105361938989742d-48, &
            0.038373959065755049d0, &
            1.7619463660345850d-8, &
            1.8599121523307042d-47, &
            0.0069412003796483614d0, &
            1.3183518729881739d-6, &
            6.0385992432172000d-45, &
            0.0069806358418043373d0, &
            1.2972579168609275d-5, &
            1.1630100138612864d-43, &
            0.056252842339029846d0, &
            5.5648811374540254d-18, &
            2.0890635084369528d-101, &
            0.014322808905506972d0, &
            2.1715475594596204d-17, &
            1.2201742646844473d-100, &
            0.00075639476867535033d0, &
            5.6795188916839963d-15, &
            1.5057590361430175d-97, &
            0.00042400877792188891d0, &
            9.8431188844873609d-14, &
            5.4746171917443813d-96]
        real(dp), parameter :: rel_tol_x = 1.0d-12, rel_tol_w = 1.0d-12

        real(dp) :: x(64), w(64)
        integer :: j, ia, k, m, n, idx(3)

        call print_test("test_gen_laguerre_sampled_large")
        m = 0
        do j = 1, size(n_list)
            n = n_list(j)
            idx = [1, n/2, n]
            do ia = 1, size(alpha_list)
                call gauss_gen_laguerre(n, alpha_list(ia), x(1:n), w(1:n))
                do k = 1, 3
                    m = m + 1
                    if (abs(x(idx(k))/x_ref(m) - 1.0_dp) > rel_tol_x .or. &
                        abs(w(idx(k))/w_ref(m) - 1.0_dp) > rel_tol_w) then
                        call print_fail
                        print *, "n=", n, " alpha=", alpha_list(ia), &
                            " i=", idx(k), x(idx(k)), w(idx(k))
                        stop 1
                    end if
                end do
            end do
        end do
        call print_ok
    end subroutine test_gen_laguerre_sampled_large

    subroutine test_gen_laguerre_moments
        ! Moments of x**alpha * exp(-x): integral of x**k equals
        ! Gamma(k + alpha + 1); exact for k <= 2n-1, capped at k = 100
        ! to keep x**k below double-precision overflow
        integer, parameter :: n_list(6) = [1, 4, 8, 16, 32, 64]
        real(dp), parameter :: rel_tol = 5.0d-12

        real(dp) :: x(64), w(64), s, ref, alpha
        integer :: j, ia, k, n

        call print_test("test_gen_laguerre_moments")
        do j = 1, size(n_list)
            n = n_list(j)
            do ia = 1, size(alpha_list)
                alpha = alpha_list(ia)
                call gauss_gen_laguerre(n, alpha, x(1:n), w(1:n))
                do k = 0, min(2*n - 1, 100)
                    s = sum(w(1:n)*x(1:n)**k)
                    ref = gamma(real(k, dp) + alpha + 1.0_dp)
                    if (abs(s/ref - 1.0_dp) > rel_tol) then
                        call print_fail
                        print *, "n=", n, " alpha=", alpha, " k=", k, &
                            " moment=", s, " ref=", ref
                        stop 1
                    end if
                end do
            end do
        end do
        call print_ok
    end subroutine test_gen_laguerre_moments

end program test_gauss_quadrature
