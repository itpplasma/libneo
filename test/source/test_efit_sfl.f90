!> Tests of the EFIT straight-field-line construction (src/efit_to_boozer/efit_sfl.f90)
!> on a synthetic Solov'ev equilibrium written as G-EQDSK by this program:
!>     psi = a R^2 Z^2 + b (R^2 - R0^2)^2 + c Z^2 + d R^2 Z    (SI, Wb/rad)
!> solves Delta* psi = -mu0 R^2 p' - F F' exactly with mu0 p' = -(8 b + 2 a) and
!> F F' = -2 c. The qpsi column comes from an independent polar-contour quadrature
!> of the analytic psi (bisection, no spline).
program test_efit_sfl
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use efit_sfl
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp, mu0 = 4.0e-7_dp*pi
    real(dp), parameter :: r0s = 1.7_dp, bb = 0.079_dp, dd = 0.01_dp, b0 = 2.0_dp, psib = 0.2_dp
    real(dp) :: aa = 0.092_dp, cc = 0.2_dp   ! a = -4 b, c = 2 for the p' = 0 file
    integer, parameter :: mp = 12
    integer :: nfail = 0
    real(dp) :: rax, zax, psiax

    call analytic_axis(rax, zax, psiax)
    call write_solovev('solovev_ok.geqdsk', 65, 'ok')
    call write_solovev('solovev_pprime.geqdsk', 65, 'pprime')
    call write_solovev('solovev_qpsi.geqdsk', 65, 'qpsi')
    call write_solovev('solovev_cocos11.geqdsk', 65, 'cocos11')
    call write_solovev('solovev_flip.geqdsk', 65, 'flip')
    aa = -4*bb
    cc = 2.0_dp
    call analytic_axis(rax, zax, psiax)
    call write_solovev('solovev_nopressure.geqdsk', 65, 'nopressure')
    aa = 0.092_dp
    cc = 0.2_dp
    call analytic_axis(rax, zax, psiax)

    call test_checks()
    call test_errors()
    call test_routes_and_boozer_property()
    call test_conventions()
    call test_resolution()
    call test_straight_current()
    call test_derivatives()

    if (nfail > 0) then
        print '(i0, a)', nfail, ' check(s) failed'
        error stop 1
    end if
    print '(a)', 'test_efit_sfl: all checks passed'

contains

    subroutine check(ok, label, value, tol)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tol
        if (ok) then
            print '(a, a, es10.2, a, es10.2)', 'PASS ', label//': ', value, ' tol ', tol
        else
            print '(a, a, es10.2, a, es10.2)', 'FAIL ', label//': ', value, ' tol ', tol
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine init(e, file, route, n_theta, allow)
        type(efit_sfl_t), intent(inout) :: e
        character(*), intent(in) :: file
        integer, intent(in) :: route, n_theta
        logical, intent(in) :: allow
        type(efit_sfl_config_t) :: cfg
        type(efit_sfl_status_t) :: st
        cfg%route = route
        cfg%n_theta = n_theta
        cfg%n_q = 64
        cfg%allow_poor_equilibrium = allow
        cfg%override_reason = 'test'
        call efit_sfl_init(e, file, cfg, st)
        if (st%code /= EFIT_SFL_OK .and. .not. st%overridden) call efit_sfl_raise(st)
    end subroutine init

    subroutine test_checks()
        type(efit_sfl_t) :: e
        type(contour_t) :: ct
        real(dp) :: err, xs(3) = [0.05_dp, 0.5_dp, 0.95_dp], psit
        integer :: k

        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_COVARIANT, 256, .false.)
        call check(e%diag%status%code == EFIT_SFL_OK, 'consistent file passes', 0.0_dp, 0.0_dp)
        call check(e%diag%eps_gs_window < 2e-6_dp, 'eps_gs (exact GS solution)', &
                   e%diag%eps_gs_window, 2e-6_dp)
        call check(e%diag%eps_q_window < 1e-8_dp, 'eps_q (analytic contour qpsi)', &
                   e%diag%eps_q_window, 1e-8_dp)
        call check(e%diag%eps_fb_median < 2e-6_dp, 'eps_fb (force balance, file p'')', &
                   e%diag%eps_fb_median, 2e-6_dp)
        call check(e%diag%exp_bpol == 0 .and. abs(e%diag%gs_factor - 1) < 1e-5_dp, &
                   'GS factor', abs(e%diag%gs_factor - 1), 1e-5_dp)
        err = 0
        do k = 1, 3
            call efit_sfl_contour_at(e, xs(k), ct)
            psit = psiax + xs(k)*(psib - psiax)
            err = max(err, abs(sum(ct%out(1, :))/size(ct%rho)/q_analytic(psit, 1024) - 1))
        end do
        ! limited by the 65 x 65 quintic spline of psi (eps_gs above)
        call check(err < 1e-9_dp, 'contour q vs analytic psi', err, 1e-9_dp)
    end subroutine test_checks

    subroutine test_errors()
        type(efit_sfl_t) :: e
        type(efit_sfl_config_t) :: cfg
        type(efit_sfl_status_t) :: st

        cfg%n_theta = 128
        cfg%n_q = 32
        call efit_sfl_init(e, 'solovev_pprime.geqdsk', cfg, st)
        call check(st%code == EFIT_SFL_ERR_GS .and. .not. st%overridden, &
                   'pprime x 0.7 raises SFL_ERR_GS', st%value, cfg%tol_gs)
        call check(index(st%message, 'Grad-Shafranov') > 0 .and. &
                   index(st%message, 'allow_poor_equilibrium') > 0, 'message names cause, override', &
                   0.0_dp, 0.0_dp)
        call efit_sfl_init(e, 'solovev_qpsi.geqdsk', cfg, st)
        call check(st%code == EFIT_SFL_ERR_Q, 'qpsi x 1.05 raises SFL_ERR_Q', st%value, cfg%tol_q)
        cfg%allow_poor_equilibrium = .true.
        cfg%override_reason = 'test override'
        call efit_sfl_init(e, 'solovev_pprime.geqdsk', cfg, st)
        call check(st%overridden .and. e%diag%override, 'override is kept and logged', &
                   st%value, cfg%tol_gs)
    end subroutine test_errors

    !> Both Boozer routes agree to round-off in axisymmetry; on the Boozer
    !> surfaces psi is constant and B_theta = B . dr/dtheta_B at fixed phi_B is
    !> a flux function (analytic field, independent of the spline).
    subroutine test_routes_and_boozer_property()
        type(efit_sfl_t) :: e
        type(efit_sfl_surface_t) :: o(3), p(3), h(3)
        real(dp) :: ss(3) = [0.1_dp, 0.5_dp, 0.9_dp], dpsi, dbth, db
        integer :: k

        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_COVARIANT, 512, .false.)
        call efit_sfl_transform(e, ss, 40, o)
        call check(maxval(o%d_family) < 1e-12_dp, 'covariant = jacobian_b2 (d_B)', &
                   maxval(o%d_family), 1e-12_dp)
        dpsi = 0
        dbth = 0
        db = 0
        do k = 1, 3
            call boozer_property(o(k), dpsi, dbth, db)
        end do
        call check(dpsi < 1e-9_dp, 'psi constant on Boozer surfaces', dpsi, 1e-9_dp)
        call check(dbth < 1e-7_dp, 'B_theta_B flux function (mpol 40)', dbth, 1e-7_dp)
        call check(db < 1e-9_dp, '|B| series vs analytic field', db, 1e-9_dp)
        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_PEST, 512, .false.)
        call efit_sfl_transform(e, ss, 40, p)
        call check(maxval(abs([p(1)%vmnc, p(2)%vmnc, p(3)%vmns])) == 0.0_dp, 'pest: nu = 0', &
                   maxval(abs([p(1)%vmnc, p(2)%vmnc])), 0.0_dp)
        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_JACOBIAN_CONST, 512, .false.)
        call efit_sfl_transform(e, ss, 40, h)
        call check(abs(h(2)%iota - o(2)%iota) < 1e-14_dp .and. &
                   maxval(abs(h(2)%bmnc(1:) - o(2)%bmnc(1:))) > 1e-4_dp*o(2)%bmnc(0), &
                   'hamada: same iota, different b_mn', abs(h(2)%iota - o(2)%iota), 1e-14_dp)
    end subroutine test_routes_and_boozer_property

    subroutine boozer_property(o, dpsi, dbth, db)
        type(efit_sfl_surface_t), intent(in) :: o
        real(dp), intent(inout) :: dpsi, dbth, db
        integer, parameter :: n = 256
        real(dp) :: th, r, z, v, rt, zt, vt, b, ps(n), bt(n), pr, pz, f
        integer :: j, m

        f = sqrt((r0s*b0)**2 - 4*cc*(o%psi*1e-8_dp - psib))
        do j = 1, n
            th = 2*pi*(j - 1)/n
            r = 0; z = 0; v = 0; rt = 0; zt = 0; vt = 0; b = 0
            do m = 0, ubound(o%rmnc, 1)
                r = r + o%rmnc(m)*cos(m*th) + o%rmns(m)*sin(m*th)
                z = z + o%zmnc(m)*cos(m*th) + o%zmns(m)*sin(m*th)
                v = v + o%vmnc(m)*cos(m*th) + o%vmns(m)*sin(m*th)
                b = b + o%bmnc(m)*cos(m*th) + o%bmns(m)*sin(m*th)
                rt = rt + m*(o%rmns(m)*cos(m*th) - o%rmnc(m)*sin(m*th))
                zt = zt + m*(o%zmns(m)*cos(m*th) - o%zmnc(m)*sin(m*th))
                vt = vt + m*(o%vmns(m)*cos(m*th) - o%vmnc(m)*sin(m*th))
            end do
            r = r*1e-2_dp; z = z*1e-2_dp; rt = rt*1e-2_dp; zt = zt*1e-2_dp
            ps(j) = psi_a(r, z)
            call grad_a(r, z, pr, pz)
            bt(j) = -pz/r*rt + pr/r*zt + f*2*pi*vt
            db = max(db, abs(sqrt(pr**2 + pz**2 + f**2)/r - b*1e-4_dp)/(o%bmnc(0)*1e-4_dp))
        end do
        dpsi = max(dpsi, (maxval(ps) - minval(ps))/(psib - psiax))
        dbth = max(dbth, (maxval(bt) - minval(bt))/abs(sum(bt)/n))
    end subroutine boozer_property

    !> psi in Wb (COCOS 11) is detected and gives the same output; flipping the
    !> sign of psi reverses theta (sine terms change sign, q unchanged).
    subroutine test_conventions()
        type(efit_sfl_t) :: e
        type(efit_sfl_surface_t) :: o(2), c(2), f(2)
        real(dp) :: ss(2) = [0.3_dp, 0.8_dp], d

        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_COVARIANT, 256, .false.)
        call efit_sfl_transform(e, ss, mp, o)
        call init(e, 'solovev_cocos11.geqdsk', EFIT_SFL_COVARIANT, 256, .false.)
        call check(e%diag%exp_bpol == 1, 'psi in Wb detected (exp_bpol = 1)', 0.0_dp, 0.0_dp)
        call efit_sfl_transform(e, ss, mp, c)
        d = maxval(abs([c(1)%bmnc - o(1)%bmnc, c(2)%rmns - o(2)%rmns]))/o(1)%bmnc(0) &
            + abs(c(2)%iota/o(2)%iota - 1)
        call check(d < 1e-9_dp, 'COCOS 11 file gives the same spectra', d, 1e-9_dp)
        call init(e, 'solovev_flip.geqdsk', EFIT_SFL_COVARIANT, 256, .false.)
        call efit_sfl_transform(e, ss, mp, f)
        d = maxval(abs([f(1)%rmnc - o(1)%rmnc, f(2)%rmns + o(2)%rmns, f(2)%bmns + o(2)%bmns])) &
            /o(1)%rmnc(0) + abs(f(2)%iota/o(2)%iota - 1) + abs(f(2)%curr_i/o(2)%curr_i - 1)
        call check(d < 1e-9_dp, 'flipped psi: theta reversed, q and I kept', d, 1e-9_dp)
    end subroutine test_conventions

    subroutine test_resolution()
        type(efit_sfl_t) :: e
        type(efit_sfl_surface_t) :: a(3), b(3)
        real(dp) :: ss(3) = [0.1_dp, 0.5_dp, 0.9_dp], d
        integer :: k

        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_COVARIANT, 256, .false.)
        call efit_sfl_transform(e, ss, 40, a)
        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_COVARIANT, 512, .false.)
        call efit_sfl_transform(e, ss, 40, b)
        d = 0
        do k = 1, 3
            d = max(d, maxval(abs([a(k)%bmnc - b(k)%bmnc, a(k)%bmns - b(k)%bmns]))/b(k)%bmnc(0), &
                    maxval(abs([a(k)%rmnc - b(k)%rmnc, a(k)%zmns - b(k)%zmns]))/b(k)%rmnc(0))
        end do
        call check(d < 1e-12_dp, 'n_theta 256 -> 512', d, 1e-12_dp)
    end subroutine test_resolution

    !> Hamada II (straight current lines) coincides with Hamada I in
    !> equilibrium (d_H, design section 9 C1-C6); without p' it is undefined.
    subroutine test_straight_current()
        type(efit_sfl_t) :: e
        type(efit_sfl_surface_t) :: o(3), h(3)
        type(efit_sfl_status_t) :: st
        real(dp) :: ss(3) = [0.1_dp, 0.5_dp, 0.9_dp], d

        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_STRAIGHT_CURRENT, 512, .false.)
        call efit_sfl_transform(e, ss, 40, o)
        call init(e, 'solovev_ok.geqdsk', EFIT_SFL_JACOBIAN_CONST, 512, .false.)
        call efit_sfl_transform(e, ss, 40, h)
        d = maxval(abs([o(2)%rmnc - h(2)%rmnc, o(2)%zmns - h(2)%zmns]))/h(2)%rmnc(0)
        call check(maxval(o%d_family) < 1e-6_dp .and. d < 1e-6_dp, &
                   'straight_current = jacobian_const in equilibrium (d_H)', &
                   max(maxval(o%d_family), d), 1e-6_dp)
        call init(e, 'solovev_nopressure.geqdsk', EFIT_SFL_STRAIGHT_CURRENT, 128, .false.)
        call efit_sfl_transform(e, ss, mp, o, st)
        call check(st%code == EFIT_SFL_ERR_NO_PRESSURE .and. index(st%message, 'p''') > 0, &
                   'p'' = 0 raises SFL_ERR_NO_PRESSURE', st%value, 1e-8_dp)
    end subroutine test_straight_current

    subroutine test_derivatives()
        integer :: route
        do route = 1, 5
            call derivatives(route)
        end do
    end subroutine test_derivatives

    subroutine derivatives(route)
        integer, intent(in) :: route
        type(efit_sfl_t) :: e, ep, em
        type(efit_sfl_surface_t) :: od(3), ob(3), op(3), om(3)
        type(efit_sfl_tangent_t) :: t, tb
        real(dp) :: ss(3) = [0.1_dp, 0.5_dp, 0.9_dp], lhs, rhs, fd, h, best
        integer :: i, j, k
        character(len=32) :: lab

        call init(e, 'solovev_ok.geqdsk', route, 128, .false.)
        allocate (t%psirz(e%nw, e%nh), t%fpol(e%nw))
        do j = 1, e%nh
            do i = 1, e%nw
                t%psirz(i, j) = (e%psi_edge - e%psi_ax)*sin(0.07_dp*i + 0.05_dp*j) &
                                *cos(0.03_dp*i*j/e%nw)
            end do
        end do
        t%fpol = 0.01_dp*e%g%fpol*cos(0.1_dp*[(i, i=1, e%nw)])
        call efit_sfl_jvp(e, ss, mp, t, od)
        do i = 1, 3
            call seed(ob(i), i)
        end do
        call efit_sfl_vjp(e, ss, mp, ob, tb)
        lhs = 0
        do i = 1, 3
            lhs = lhs + dotp(od(i), ob(i))
        end do
        rhs = sum(t%psirz*tb%psirz) + sum(t%fpol*tb%fpol)
        write (lab, '(a, i0)') 'adjoint duality, route ', route
        call check(abs(lhs - rhs) < 1e-12_dp*abs(lhs), trim(lab), abs(lhs - rhs)/abs(lhs), 1e-12_dp)
        best = huge(1.0_dp)
        do k = 3, 8
            h = 10.0_dp**(-k)
            call efit_sfl_perturb(e, t, h, ep)
            call efit_sfl_perturb(e, t, -h, em)
            call efit_sfl_transform(ep, ss, mp, op)
            call efit_sfl_transform(em, ss, mp, om)
            fd = 0
            do i = 1, 3
                fd = fd + (dotp(op(i), ob(i)) - dotp(om(i), ob(i)))/(2*h)
            end do
            best = min(best, abs(fd - lhs)/abs(lhs))
        end do
        write (lab, '(a, i0)') 'FD V-curve minimum, route ', route
        call check(best < 1e-8_dp, trim(lab), best, 1e-8_dp)
        ! Gram form <J t, J t> = <t, J^T J t>, free of the cancellation of the
        ! fixed seed (3e3 for route 4, whose divisor is proportional to p')
        call efit_sfl_vjp(e, ss, mp, od, tb)
        lhs = 0
        do i = 1, 3
            lhs = lhs + dotp(od(i), od(i))
        end do
        rhs = sum(t%psirz*tb%psirz) + sum(t%fpol*tb%fpol)
        write (lab, '(a, i0)') 'adjoint duality (Gram), route ', route
        call check(abs(lhs - rhs) < 1e-13_dp*lhs, trim(lab), abs(lhs - rhs)/lhs, 1e-13_dp)
    end subroutine derivatives

    subroutine seed(b, i)
        type(efit_sfl_surface_t), intent(inout) :: b
        integer, intent(in) :: i
        integer :: m
        allocate (b%rmnc(0:mp), b%rmns(0:mp), b%zmnc(0:mp), b%zmns(0:mp), b%vmnc(0:mp), &
                  b%vmns(0:mp), b%bmnc(0:mp), b%bmns(0:mp))
        do m = 0, mp
            b%rmnc(m) = cos(1.0_dp*m + i); b%rmns(m) = sin(2.0_dp*m + i)
            b%zmnc(m) = cos(3.0_dp*m - i); b%zmns(m) = sin(0.5_dp*m*i)
            b%vmnc(m) = 100*cos(0.7_dp*m); b%vmns(m) = 100*sin(0.3_dp*m + i)
            b%bmnc(m) = 1e-3_dp*cos(1.3_dp*m + i); b%bmns(m) = 1e-3_dp*sin(0.2_dp*m - i)
        end do
        b%iota = 1.3_dp*i; b%curr_i = 1e-4_dp*i; b%curr_g = 2e-6_dp
        b%x = 0.7_dp; b%psi = 1e-7_dp*i
    end subroutine seed

    real(dp) function dotp(a, b)
        type(efit_sfl_surface_t), intent(in) :: a, b
        dotp = sum(a%rmnc*b%rmnc) + sum(a%rmns*b%rmns) + sum(a%zmnc*b%zmnc) &
               + sum(a%zmns*b%zmns) + sum(a%vmnc*b%vmnc) + sum(a%vmns*b%vmns) &
               + sum(a%bmnc*b%bmnc) + sum(a%bmns*b%bmns) + a%iota*b%iota &
               + a%curr_i*b%curr_i + a%curr_g*b%curr_g + a%x*b%x + a%psi*b%psi
    end function dotp

    ! ------------------------------------------------ analytic Solov'ev field

    pure real(dp) function psi_a(r, z)
        real(dp), intent(in) :: r, z
        psi_a = aa*r**2*z**2 + bb*(r**2 - r0s**2)**2 + cc*z**2 + dd*r**2*z
    end function psi_a

    pure subroutine grad_a(r, z, pr, pz)
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: pr, pz
        pr = 2*aa*r*z**2 + 4*bb*r*(r**2 - r0s**2) + 2*dd*r*z
        pz = 2*aa*r**2*z + 2*cc*z + dd*r**2
    end subroutine grad_a

    subroutine analytic_axis(r, z, p)
        real(dp), intent(out) :: r, z, p
        real(dp) :: pr, pz, hrr, hrz, hzz, det
        integer :: it
        r = r0s
        z = 0
        do it = 1, 50
            call grad_a(r, z, pr, pz)
            hrr = 2*aa*z**2 + 12*bb*r**2 - 4*bb*r0s**2 + 2*dd*z
            hrz = 4*aa*r*z + 2*dd*r
            hzz = 2*aa*r**2 + 2*cc
            det = hrr*hzz - hrz**2
            r = r - (hzz*pr - hrz*pz)/det
            z = z - (-hrz*pr + hrr*pz)/det
        end do
        p = psi_a(r, z)
    end subroutine analytic_axis

    !> q on psi = pt from rays about the axis, bisection on the analytic psi.
    real(dp) function q_analytic(pt, n)
        real(dp), intent(in) :: pt
        integer, intent(in) :: n
        real(dp) :: th, lo, hi, mid, r, z, pr, pz, f, acc
        integer :: j, it
        f = sqrt((r0s*b0)**2 - 4*cc*(pt - psib))
        acc = 0
        do j = 1, n
            th = 2*pi*(j - 1)/n
            lo = 0
            hi = 1.5_dp
            do it = 1, 200
                mid = 0.5_dp*(lo + hi)
                if (psi_a(rax + mid*cos(th), zax + mid*sin(th)) > pt) then
                    hi = mid
                else
                    lo = mid
                end if
            end do
            r = rax + lo*cos(th)
            z = zax + lo*sin(th)
            call grad_a(r, z, pr, pz)
            acc = acc + lo/(r*(pr*cos(th) + pz*sin(th)))
        end do
        q_analytic = f*acc/n
    end function q_analytic

    subroutine write_solovev(path, nw, variant)
        character(*), intent(in) :: path, variant
        integer, intent(in) :: nw
        real(dp) :: rg(nw), zg(nw), ps(nw), psirz(nw, nw), fpol(nw), pres(nw), ffp(nw), &
                    ppr(nw), q(nw), pp, sgn, scale
        integer :: i, j, u

        do i = 1, nw
            rg(i) = 0.9_dp + 1.5_dp*(i - 1)/(nw - 1)
            zg(i) = -0.9_dp + 1.8_dp*(i - 1)/(nw - 1)
            ps(i) = psiax + (psib - psiax)*(i - 1)/(nw - 1)
        end do
        do j = 1, nw
            do i = 1, nw
                psirz(i, j) = psi_a(rg(i), zg(j))
            end do
        end do
        pp = -(8*bb + 2*aa)/mu0   ! 0 for 'nopressure'
        fpol = sqrt((r0s*b0)**2 - 4*cc*(ps - psib))
        pres = pp*(ps - psib)
        ffp = -2*cc
        ppr = pp
        do i = 1, nw
            q(i) = q_analytic(psiax + max(real(i - 1, dp)/(nw - 1), 1e-6_dp)*(psib - psiax), 256)
        end do
        sgn = 1
        scale = 1
        select case (variant)
        case ('pprime')
            ppr = 0.7_dp*ppr
            pres = 0.7_dp*pres
        case ('qpsi')
            q = 1.05_dp*q
        case ('cocos11')
            scale = 2*pi
        case ('flip')
            sgn = -1
        end select
        open (newunit=u, file=path, status='replace')
        write (u, '(a48, 3i4)') 'solovev '//variant, 0, nw, nw
        write (u, '(5es16.9)') rg(nw) - rg(1), zg(nw) - zg(1), r0s, rg(1), 0.0_dp
        write (u, '(5es16.9)') rax, zax, sgn*scale*psiax, sgn*scale*psib, b0
        write (u, '(5es16.9)') 1e6_dp, sgn*scale*psiax, 0.0_dp, rax, 0.0_dp
        write (u, '(5es16.9)') zax, 0.0_dp, sgn*scale*psib, 0.0_dp, 0.0_dp
        write (u, '(5es16.9)') fpol
        write (u, '(5es16.9)') pres
        write (u, '(5es16.9)') sgn*ffp/scale
        write (u, '(5es16.9)') sgn*ppr/scale
        write (u, '(5es16.9)') sgn*scale*psirz
        write (u, '(5es16.9)') q
        if (variant == 'flip') then
            ! no limiter record (as FreeGS writes it)
            write (u, '(2i5)') 1, 0
            write (u, '(5es16.9)') rax, zax
        else
            write (u, '(2i5)') 1, 1
            write (u, '(5es16.9)') rax, zax
            write (u, '(5es16.9)') rg(1), zg(1)
        end if
        close (u)
    end subroutine write_solovev

end program test_efit_sfl
