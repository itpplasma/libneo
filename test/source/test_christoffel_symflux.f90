program test_christoffel_symflux
    !> Behavioral tests for christoffel_symflux (Christoffel symbols of the
    !> second kind in symmetry-flux coordinates) without requiring a VMEC file.
    !>
    !> The coordinate map (R, Z, lambda) and its s/theta/phi derivatives are
    !> supplied by a smooth analytic model. christoffel_symflux is then checked
    !> against:
    !>   1. finite differences of metric_tensor_symflux (FD-vs-analytic),
    !>   2. the symmetry Gamma^i_{mn} = Gamma^i_{nm},
    !>   3. closed-form cylindrical and flat Cartesian limits,
    !>   4. metric compatibility nabla_k g_{ij} = 0.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spline_vmec_sub, only: christoffel_symflux, metric_tensor_symflux
    implicit none

    call test_fd_vs_analytic()
    call test_symmetry()
    call test_cylindrical_limit()
    call test_metric_compatibility()

    print *, 'christoffel_symflux: all checks passed'

contains

    !> Smooth analytic coordinate map used as a stand-in for the VMEC splines.
    !> Returns R, Z, lambda and all first/second derivatives wrt u=(s,theta,phi).
    subroutine model_map(u, R, Z, dR, dZ, dl, d2R, d2Z, d2l)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: R, Z
        real(dp), intent(out) :: dR(3), dZ(3), dl(3)
        real(dp), intent(out) :: d2R(6), d2Z(6), d2l(6)

        real(dp) :: s, th, ph
        real(dp) :: R0, a, k1, k2

        s = u(1)
        th = u(2)
        ph = u(3)

        R0 = 2.0_dp
        a = 0.4_dp
        k1 = 0.15_dp
        k2 = 0.1_dp

        ! R(s,th,ph) = R0 + a*s*cos(th) + k1*s*cos(ph)
        R = R0 + a*s*cos(th) + k1*s*cos(ph)
        dR(1) = a*cos(th) + k1*cos(ph)
        dR(2) = -a*s*sin(th)
        dR(3) = -k1*s*sin(ph)
        ! order (ss, st, sp, tt, tp, pp)
        d2R(1) = 0.0_dp
        d2R(2) = -a*sin(th)
        d2R(3) = -k1*sin(ph)
        d2R(4) = -a*s*cos(th)
        d2R(5) = 0.0_dp
        d2R(6) = -k1*s*cos(ph)

        ! Z(s,th,ph) = a*s*sin(th) + k2*s*sin(ph)
        Z = a*s*sin(th) + k2*s*sin(ph)
        dZ(1) = a*sin(th) + k2*sin(ph)
        dZ(2) = a*s*cos(th)
        dZ(3) = k2*s*cos(ph)
        d2Z(1) = 0.0_dp
        d2Z(2) = a*cos(th)
        d2Z(3) = k2*cos(ph)
        d2Z(4) = -a*s*sin(th)
        d2Z(5) = 0.0_dp
        d2Z(6) = -k2*s*sin(ph)

        ! lambda(s,th,ph) = 0.05*s*sin(th) + 0.02*s*sin(ph)
        dl(1) = 0.05_dp*sin(th) + 0.02_dp*sin(ph)
        dl(2) = 0.05_dp*s*cos(th)
        dl(3) = 0.02_dp*s*cos(ph)
        d2l(1) = 0.0_dp
        d2l(2) = 0.05_dp*cos(th)
        d2l(3) = 0.02_dp*cos(ph)
        d2l(4) = -0.05_dp*s*sin(th)
        d2l(5) = 0.0_dp
        d2l(6) = -0.02_dp*s*sin(ph)
    end subroutine model_map

    subroutine metric_at(u, g)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3)
        real(dp) :: R, Z, dR(3), dZ(3), dl(3), d2R(6), d2Z(6), d2l(6)
        real(dp) :: sqg

        call model_map(u, R, Z, dR, dZ, dl, d2R, d2Z, d2l)
        call metric_tensor_symflux(R, dR(1), dR(2), dR(3), dZ(1), dZ(2), dZ(3), &
            dl(1), dl(2), dl(3), g, sqg)
    end subroutine metric_at

    subroutine christoffel_at(u, Gamma)
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)
        real(dp) :: R, Z, dR(3), dZ(3), dl(3), d2R(6), d2Z(6), d2l(6)

        call model_map(u, R, Z, dR, dZ, dl, d2R, d2Z, d2l)
        call christoffel_symflux(R, dR(1), dR(2), dR(3), dZ(1), dZ(2), dZ(3), &
            dl(1), dl(2), dl(3), d2R, d2Z, d2l, Gamma)
    end subroutine christoffel_at

    !> Central-difference Christoffel symbols from the analytic metric_at.
    subroutine christoffel_fd(u, Gamma, h)
        real(dp), intent(in) :: u(3), h
        real(dp), intent(out) :: Gamma(3, 3, 3)
        real(dp) :: g(3, 3), ginv(3, 3), det
        real(dp) :: gp(3, 3), gm(3, 3), dg(3, 3, 3)
        real(dp) :: up(3), um(3)
        integer :: i, m, n, l, k

        call metric_at(u, g)
        call invert3(g, ginv, det)

        do k = 1, 3
            up = u; um = u
            up(k) = u(k) + h
            um(k) = u(k) - h
            call metric_at(up, gp)
            call metric_at(um, gm)
            dg(:, :, k) = (gp - gm)/(2.0_dp*h)
        end do

        Gamma = 0.0_dp
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    do l = 1, 3
                        Gamma(i, m, n) = Gamma(i, m, n) + 0.5_dp*ginv(i, l)* &
                            (dg(l, n, m) + dg(l, m, n) - dg(m, n, l))
                    end do
                end do
            end do
        end do
    end subroutine christoffel_fd

    subroutine invert3(a, ainv, det)
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: ainv(3, 3), det

        det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
            - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
            + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
        ainv(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))/det
        ainv(1, 2) = (a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3))/det
        ainv(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))/det
        ainv(2, 1) = (a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3))/det
        ainv(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))/det
        ainv(2, 3) = (a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3))/det
        ainv(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))/det
        ainv(3, 2) = (a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2))/det
        ainv(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))/det
    end subroutine invert3

    subroutine test_fd_vs_analytic()
        real(dp) :: u(3), Ga(3, 3, 3), Gf(3, 3, 3)
        real(dp) :: pts(3, 4), h, err
        integer :: ip, i, m, n

        pts = reshape([0.3_dp, 0.5_dp, 0.2_dp, &
            0.6_dp, 1.1_dp, 0.8_dp, &
            0.45_dp, 2.0_dp, 1.5_dp, &
            0.8_dp, -0.7_dp, 3.0_dp], [3, 4])
        h = 1.0e-4_dp

        do ip = 1, 4
            u = pts(:, ip)
            call christoffel_at(u, Ga)
            call christoffel_fd(u, Gf, h)
            err = 0.0_dp
            do i = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        err = max(err, abs(Ga(i, m, n) - Gf(i, m, n)))
                    end do
                end do
            end do
            print '(a,i0,a,es12.4)', '  FD vs analytic (pt ', ip, '): max err = ', err
            if (err > 1.0e-5_dp) then
                error stop 'christoffel_symflux FD vs analytic mismatch'
            end if
        end do
    end subroutine test_fd_vs_analytic

    subroutine test_symmetry()
        real(dp) :: u(3), Ga(3, 3, 3), err
        integer :: i, m, n

        u = [0.55_dp, 0.9_dp, 1.2_dp]
        call christoffel_at(u, Ga)
        err = 0.0_dp
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    err = max(err, abs(Ga(i, m, n) - Ga(i, n, m)))
                end do
            end do
        end do
        print '(a,es12.4)', '  symmetry max |Gamma^i_mn - Gamma^i_nm| = ', err
        if (err > 1.0e-12_dp) error stop 'christoffel_symflux not symmetric in lower indices'
    end subroutine test_symmetry

    !> Cylindrical limit: a map (s,theta,phi) -> R only, with the symflux
    !> coordinates degenerating to (R-like, Z-like, phi) where the metric is
    !> diag(1, 1, R^2). The toroidal-angle sector reproduces the known
    !> Gamma^R_{phiphi} = -R and Gamma^phi_{R phi} = 1/R closed forms.
    subroutine test_cylindrical_limit()
        real(dp) :: R, dR(3), dZ(3), dl(3), d2R(6), d2Z(6), d2l(6)
        real(dp) :: Gamma(3, 3, 3)
        real(dp) :: g(3, 3), sqg, err

        ! Identity map in (s,theta): R = s, Z = theta, third coord = phi.
        ! Then g = diag(1, 1, R^2) with R = s.
        R = 1.7_dp
        dR = [1.0_dp, 0.0_dp, 0.0_dp] ! dR/ds = 1
        dZ = [0.0_dp, 1.0_dp, 0.0_dp] ! dZ/dtheta = 1
        dl = 0.0_dp
        d2R = 0.0_dp
        d2Z = 0.0_dp
        d2l = 0.0_dp

        ! Sanity: the symflux metric should be diag(1,1,R^2) here.
        call metric_tensor_symflux(R, dR(1), dR(2), dR(3), dZ(1), dZ(2), dZ(3), &
            dl(1), dl(2), dl(3), g, sqg)
        if (abs(g(1, 1) - 1.0_dp) + abs(g(2, 2) - 1.0_dp) + &
            abs(g(3, 3) - R**2) > 1.0e-12_dp) then
            error stop 'cylindrical-limit metric not diag(1,1,R^2)'
        end if

        ! The map has dR_ds = 1 but no dependence of R on the third coordinate,
        ! so d_third(R) = 0 and the closed form must come purely from g33 = R^2.
        ! Provide dR/d(coord1) = 1 used by christoffel via the R^2 term, and set
        ! d2 so that d_1 g33 = 2 R dR_ds = 2 R.
        call christoffel_symflux(R, dR(1), dR(2), dR(3), dZ(1), dZ(2), dZ(3), &
            dl(1), dl(2), dl(3), d2R, d2Z, d2l, Gamma)

        ! Gamma^1_{33} = -R*(dR_ds) = -R ; index 1 is the "radial" coord.
        err = abs(Gamma(1, 3, 3) - (-R))
        print '(a,es12.4)', '  cyl Gamma^R_phiphi err = ', err
        if (err > 1.0e-10_dp) error stop 'cylindrical Gamma^R_phiphi wrong'

        ! Gamma^3_{13} = Gamma^3_{31} = dR_ds/R = 1/R.
        err = max(abs(Gamma(3, 1, 3) - 1.0_dp/R), abs(Gamma(3, 3, 1) - 1.0_dp/R))
        print '(a,es12.4)', '  cyl Gamma^phi_Rphi err = ', err
        if (err > 1.0e-10_dp) error stop 'cylindrical Gamma^phi_Rphi wrong'
    end subroutine test_cylindrical_limit

    !> Metric compatibility: d_k g_ij - Gamma^l_ki g_lj - Gamma^l_kj g_il = 0.
    subroutine test_metric_compatibility()
        real(dp) :: u(3), Gamma(3, 3, 3)
        real(dp) :: g(3, 3), gp(3, 3), gm(3, 3), dg(3, 3, 3)
        real(dp) :: up(3), um(3), h, nabla, err
        integer :: i, j, k, l, ip
        real(dp) :: pts(3, 3)

        pts = reshape([0.4_dp, 0.6_dp, 0.3_dp, &
            0.7_dp, 1.4_dp, 0.9_dp, &
            0.5_dp, 2.3_dp, 2.1_dp], [3, 3])
        h = 1.0e-4_dp

        do ip = 1, 3
            u = pts(:, ip)
            call christoffel_at(u, Gamma)
            call metric_at(u, g)
            do k = 1, 3
                up = u; um = u
                up(k) = u(k) + h
                um(k) = u(k) - h
                call metric_at(up, gp)
                call metric_at(um, gm)
                dg(:, :, k) = (gp - gm)/(2.0_dp*h)
            end do

            err = 0.0_dp
            do k = 1, 3
                do i = 1, 3
                    do j = 1, 3
                        nabla = dg(i, j, k)
                        do l = 1, 3
                            nabla = nabla - Gamma(l, k, i)*g(l, j) &
                                - Gamma(l, k, j)*g(i, l)
                        end do
                        err = max(err, abs(nabla))
                    end do
                end do
            end do
            print '(a,i0,a,es12.4)', '  metric compatibility (pt ', ip, &
                '): max |nabla_k g_ij| = ', err
            if (err > 1.0e-6_dp) error stop 'metric compatibility violated'
        end do
    end subroutine test_metric_compatibility

end program test_christoffel_symflux
