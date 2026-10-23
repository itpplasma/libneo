!> Straight-field-line routes of the hybrid construction (boozer_hybrid and
!> its submodule boozer_hybrid_routes).
!>
!> Checks:
!>  - the FortSym route kernels: sfl_route_weight and sfl_beta against direct
!>    formulas, sfl_radial_field against the radial derivative of the field
!>    kernel boozer_hybrid_field along a path x(s), Jacobians against central
!>    differences;
!>  - each route enforces its property: with the nu truncation raised
!>    (mnu = 6 (mpol - 1), nnu = 6 ntor) the relative variation of J* B^2
!>    (sfl_jacobian_b2) and of J* (sfl_jacobian_const) is below 1e-6, nu = 0
!>    for sfl_pest; on li383 (p' /= 0) the analytic nu_s, I' and G' of the
!>    Hamada II route match central differences of the covariant route;
!>  - sfl_straight_current on li383 has <eta>_B = 0 and the analytic
!>    G' + iota I' satisfies the averaged radial force balance to the level of
!>    the equilibrium;
!>  - doubling the angle grid changes the spectrum of every route by less than
!>    1e-10 of B_00 at s <= 0.5 and 1e-8 at s = 0.85;
!>  - the JVP of every route matches central differences of the primal
!>    (V-curve minimum below 1e-8) and the VJP is its adjoint to 1e-13;
!> on QA (p = 0), li383 (p' /= 0, net current), the circular tokamak
!> (axisymmetric) and the up-down asymmetric tokamak (lasym). Agreement with the
!> independent JAX implementation and with GVEC is checked in boozer-bench.
program test_sfl_routes
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, boozer_ad_primal, &
                              boozer_ad_jvp, boozer_ad_vjp, wout_tangent_alloc, &
                              spectrum_alloc, boozer_ad_perturb
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, boozer_hybrid_init, &
                             hybrid_surface, hybrid_diag_t, hybrid_surface_diag, &
                             hybrid_field_kernel, sfl_route_name, route_weight_kernel, &
                             route_radial_kernel, route_beta_kernel, SFL_COVARIANT, &
                             SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST, SFL_STRAIGHT_CURRENT, SFL_PEST
    use boozer_hybrid_quality, only: sfl_quality_config_t
    implicit none

    real(dp), parameter :: s_list(3) = [0.15_dp, 0.5_dp, 0.85_dp]
    real(dp), parameter :: MU0 = 16.0e-7_dp*atan(1.0_dp)
    !> The test wouts are the coarse public bench files (ns <= 50), which fail
    !> the equilibrium-quality check; the transform itself is tested here.
    type(sfl_quality_config_t) :: coarse
    logical :: failed
    integer(int64) :: seed_state

    failed = .false.
    coarse%allow_poor_equilibrium = .true.
    coarse%override_reason = 'coarse public test wout (ns <= 50)'
    seed_state = 4242
    call check_kernels()
    call check_properties('wout.nc', 'QA', [2, 3, 5])
    call check_properties('wout_li383.nc', 'li383', [2, 3, 4, 5])
    call check_properties('wout_circular_tokamak.nc', 'circular tokamak', [2, 3, 5])
    call check_properties('wout_updown_tokamak.nc', 'up-down tokamak', [2, 3, 5])
    call check_hamada2('wout_li383.nc')
    call check_ad('wout.nc', 8, 6, 'QA', [2, 3, 5])
    call check_ad('wout_li383.nc', 8, 6, 'li383', [2, 3, 4, 5])
    call check_ad('wout_circular_tokamak.nc', 8, 0, 'circular tokamak', [2, 3, 5])
    call check_ad('wout_updown_tokamak.nc', 8, 0, 'up-down tokamak', [2, 3])
    if (failed) error stop 'test_sfl_routes FAILED'
    print *, 'test_sfl_routes PASSED'

contains

    subroutine report(name, err, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol

        if (err <= tol) then
            print '(a, es10.3, a, es9.2)', '  ok   '//name//': ', err, ' <= ', tol
        else
            print '(a, es10.3, a, es9.2)', '  FAIL '//name//': ', err, ' > ', tol
            failed = .true.
        end if
    end subroutine report

    real(dp) function urand()
        seed_state = modulo(1103515245_int64*seed_state + 12345_int64, 2147483647_int64)
        urand = 2.0_dp*real(seed_state, dp)/2147483647.0_dp - 1.0_dp
    end function urand

    !> The route kernels of src/generated/sfl_route_kernels_leaves.inc.
    subroutine check_kernels()
        real(dp) :: x(10), f(4), jac(4, 10), fp(4), fm(4)
        real(dp) :: rin(19), rp(19), rf(3), rj(3, 19), rfp(3), rfm(3), rjt(3, 19)
        real(dp) :: bin(12), bp(12), bf, bj(12), bfp, bfm, bjt(12), wf(2), wj(2, 6), wfp(2)
        real(dp) :: wfm(2), wjt(2, 6), win(6), wp(6)
        real(dp) :: h, e_val, e_fd, e_rad, e_beta, ds, a, b, c, e, ts, ths, phs
        integer :: k, ip

        e_val = 0.0_dp
        e_fd = 0.0_dp
        e_rad = 0.0_dp
        e_beta = 0.0_dp
        do ip = 1, 4
            x = [1.5_dp, 0.3_dp, -0.2_dp, 0.05_dp, 0.1_dp, 0.4_dp, 0.03_dp, 0.1_dp, &
                 -0.05_dp, 0.4_dp]
            do k = 1, 10
                x(k) = x(k)*(1.0_dp + 0.3_dp*urand())
            end do
            ! weight kernel against the direct formula
            call hybrid_field_kernel(x, f, jac)
            win = [x(1), x(2), x(3), x(5), x(6), f(1)]
            call weight(win, wf, wj)
            e_val = max(e_val, abs(wf(1) - x(1)*(x(3)*x(5) - x(2)*x(6)))/abs(wf(1)), &
                        abs(wf(2) - wf(1)*f(1)**2)/abs(wf(2)))
            ! radial kernel: d/ds of the field kernel along the path x(s) whose
            ! derivative is (rs, rss, rst, rsp, zss, zst, zsp, lst, lsp, iota')
            rin(1:10) = x
            do k = 11, 19
                rin(k) = 0.2_dp*urand()
            end do
            call radial(rin, rf, rj)
            ds = 1.0e-5_dp
            call field_at(rin, ds, fp)
            call field_at(rin, -ds, fm)
            e_rad = max(e_rad, abs((fp(2) - fm(2))/(2*ds) - rf(2))/abs(rf(2)), &
                        abs((fp(3) - fm(3))/(2*ds) - rf(3))/abs(rf(3)))
            ! beta kernel against the 2x2 solve
            do k = 1, 12
                bin(k) = 0.3_dp*urand()
            end do
            bin(2) = 1.0_dp + bin(2)
            bin(3) = 3.0_dp + bin(3)
            call betak(bin, bf, bj)
            a = 1 + bin(4) + bin(11)*bin(8)
            b = bin(5) + bin(11)*bin(9)
            c = bin(8)
            e = 1 + bin(9)
            ts = bin(6) + bin(12)*bin(7) + bin(11)*bin(10)
            ths = -(e*ts - b*bin(10))/(a*e - b*c)
            phs = -(a*bin(10) - c*ts)/(a*e - b*c)
            e_beta = max(e_beta, abs(bf - (bin(1) + bin(2)*ths + bin(3)*phs))/abs(bf))
            ! Jacobians against central differences
            h = 1.0e-6_dp
            do k = 1, 6
                wp = win
                wp(k) = win(k) + h
                call weight(wp, wfp, wjt)
                wp(k) = win(k) - h
                call weight(wp, wfm, wjt)
                e_fd = max(e_fd, maxval(abs((wfp - wfm)/(2*h) - wj(:, k)))/maxval(abs(wj)))
            end do
            do k = 1, 19
                rp = rin
                rp(k) = rin(k) + h
                call radial(rp, rfp, rjt)
                rp(k) = rin(k) - h
                call radial(rp, rfm, rjt)
                e_fd = max(e_fd, maxval(abs((rfp - rfm)/(2*h) - rj(:, k)))/maxval(abs(rj)))
            end do
            do k = 1, 12
                bp = bin
                bp(k) = bin(k) + h
                call betak(bp, bfp, bjt)
                bp(k) = bin(k) - h
                call betak(bp, bfm, bjt)
                e_fd = max(e_fd, abs((bfp - bfm)/(2*h) - bj(k))/maxval(abs(bj)))
            end do
        end do
        call report('route weights = sqrt(g)_V (B^2)', e_val, 1.0e-14_dp)
        call report('radial field kernel = d/ds of the field kernel (FD)', e_rad, 1.0e-8_dp)
        call report('beta kernel = 2x2 solve', e_beta, 1.0e-14_dp)
        call report('route kernel Jacobians vs central differences', e_fd, 1.0e-8_dp)
    end subroutine check_kernels

    !> Field kernel at x(s) = x + s xs + s^2/2 xss for the inputs of the radial
    !> kernel (base x(1:10), radial derivatives in rin).
    subroutine field_at(rin, s, f)
        real(dp), intent(in) :: rin(19), s
        real(dp), intent(out) :: f(4)
        real(dp) :: x(10), xs(10), xss(10), jac(4, 10)

        x = rin(1:10)
        ! d/ds of (r, rs, rt, rp, zs, zt, zp, lt, lp, iota)
        xs = [rin(2), rin(11), rin(12), rin(13), rin(14), rin(15), rin(16), rin(17), &
              rin(18), rin(19)]
        ! r enters with its second derivative rs' = rss; the rest linearly
        xss = 0.0_dp
        xss(1) = rin(11)
        call hybrid_field_kernel(x + s*xs + 0.5_dp*s**2*xss, f, jac)
    end subroutine field_at

    subroutine weight(x, f, jac)
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: f(2), jac(2, 6)

        call route_weight_kernel(x, f, jac)
    end subroutine weight

    subroutine radial(x, f, jac)
        real(dp), intent(in) :: x(19)
        real(dp), intent(out) :: f(3), jac(3, 19)

        call route_radial_kernel(x, f, jac)
    end subroutine radial

    subroutine betak(x, f, jac)
        real(dp), intent(in) :: x(12)
        real(dp), intent(out) :: f, jac(12)

        call route_beta_kernel(x, f, jac)
    end subroutine betak

    subroutine check_properties(wout, label, routes)
        character(*), intent(in) :: wout, label
        integer, intent(in) :: routes(:)
        type(boozer_hybrid_t) :: h
        type(hybrid_diag_t) :: d
        type(hybrid_surface_t) :: sv
        integer :: ir, j, mmax, nmax, nt, np
        real(dp) :: prop
        character(len=120) :: tag

        print '(a)', ' properties: '//label//' ('//wout//')'
        call boozer_hybrid_init(h, wout, 8, 6)
        mmax = h%mmax
        nmax = h%nmax
        nt = 24*mmax + 16
        np = merge(1, 24*nmax + 16, h%axisym)
        do ir = 1, size(routes)
            call boozer_hybrid_init(h, wout, 8, 6, n_theta=nt, n_phi=np, route=routes(ir), &
                                    mnu=6*mmax, nnu=6*nmax)
            do j = 1, size(s_list)
                call hybrid_surface_diag(h, s_list(j), d)
                select case (routes(ir))
                case (SFL_JACOBIAN_B2)
                    prop = d%d_jacb2
                case (SFL_JACOBIAN_CONST)
                    prop = d%d_jac
                case (SFL_PEST)
                    call hybrid_surface(h, s_list(j), sv, with_spectrum=.false.)
                    prop = maxval(abs(sv%nu)) + maxval(abs(sv%nut)) + maxval(abs(sv%nup))
                case default
                    prop = 0.0_dp
                end select
                write (tag, '(a, f4.2, a, 4es9.2, a)') trim(sfl_route_name(routes(ir)))// &
                    ' enforced property at s = ', s_list(j), ' (jac jb2 cov cur', d%d_jac, &
                    d%d_jacb2, d%d_cov, d%d_cur, ')'
                call report(trim(tag), prop, 1.0e-6_dp)
            end do
        end do
    end subroutine check_properties

    !> Analytic radial derivatives of the Hamada II route against central
    !> differences of the covariant route.
    subroutine check_hamada2(wout)
        character(*), intent(in) :: wout
        type(boozer_hybrid_t) :: h4, h1
        type(hybrid_surface_t) :: sv, svp, svm
        type(hybrid_diag_t) :: d
        real(dp) :: ds, err_nus, err_ip, err_gp, err_eta, err_fb, nscale
        integer :: j

        print '(a)', ' Hamada II radial derivatives: li383 ('//wout//')'
        call boozer_hybrid_init(h4, wout, 8, 6, route=SFL_STRAIGHT_CURRENT)
        call boozer_hybrid_init(h1, wout, 8, 6, route=SFL_COVARIANT, n_theta=h4%nt, &
                                n_phi=h4%np)
        ds = 1.0e-4_dp
        err_nus = 0.0_dp
        err_ip = 0.0_dp
        err_gp = 0.0_dp
        err_eta = 0.0_dp
        err_fb = 0.0_dp
        do j = 1, size(s_list)
            call hybrid_surface(h4, s_list(j), sv)
            call hybrid_surface(h1, s_list(j) + ds, svp)
            call hybrid_surface(h1, s_list(j) - ds, svm)
            nscale = maxval(abs(sv%nus))
            err_nus = max(err_nus, maxval(abs((svp%nu - svm%nu)/(2*ds) - sv%nus))/nscale)
            err_ip = max(err_ip, abs((svp%ci - svm%ci)/(2*ds) - sv%cip)/abs(sv%cip))
            err_gp = max(err_gp, abs((svp%cg - svm%cg)/(2*ds) - sv%cgp)/abs(sv%cgp))
            err_eta = max(err_eta, abs(sum(sv%db*sv%eta))/sum(abs(sv%db*sv%eta)))
            call hybrid_surface_diag(h4, s_list(j), d)
            err_fb = max(err_fb, d%eps_fb)
            print '(a, f5.2, a, es10.3, a, es10.3)', '       s = ', s_list(j), &
                ': d_cur ', d%d_cur, ', eps_fb ', d%eps_fb
        end do
        call report('nu_s of the covariant route vs central differences', err_nus, 1.0e-6_dp)
        call report('I'' vs central differences of I', err_ip, 1.0e-6_dp)
        call report('G'' vs central differences of G', err_gp, 1.0e-6_dp)
        call report('<eta>_B = 0', err_eta, 1.0e-12_dp)
        ! li383 (ns = 49): eps_fb 6e-5 to 1.4e-3 (VMEC's radial resolution)
        call report('averaged force balance mu0 p'' <sqrt g> = -psi''(G'' + iota I'')', &
                    err_fb, 1.0e-2_dp)
    end subroutine check_hamada2

    subroutine check_ad(wout, mboz, nboz, label, routes)
        character(*), intent(in) :: wout, label
        integer, intent(in) :: mboz, nboz, routes(:)
        type(boozer_ad_context_t) :: ctx, cfine, cp, cm
        type(boozer_spectrum_t) :: y, yf, yd, yp, ym, yb, yv, yfd
        type(boozer_wout_tangent_t) :: t, tb
        real(dp) :: b00(size(s_list)), err, eps, lhs, rhs, best, fr(6), sj, err3, err4
        integer :: j, k, i, ir
        character(len=120) :: tag

        do ir = 1, size(routes)
            print '(a)', ' '//trim(sfl_route_name(routes(ir)))//': '//label//' ('//wout//')'
            call boozer_ad_init(ctx, wout, mboz, nboz, route=routes(ir), quality=coarse)
            call boozer_ad_primal(ctx, s_list, y)
            b00 = y%bmnc(1, :)
            call boozer_ad_init(cfine, wout, mboz, nboz, n_theta=2*ctx%grid%n_theta, &
                                n_phi=merge(1, 2*ctx%grid%n_phi, ctx%axisym), &
                                route=routes(ir), quality=coarse)
            call boozer_ad_primal(cfine, s_list, yf)
            do j = 1, size(s_list)
                err = max(maxval(abs(yf%bmnc(:, j) - y%bmnc(:, j))), &
                          maxval(abs(yf%bmns(:, j) - y%bmns(:, j))))/b00(j)
                write (tag, '(a, f4.2)') 'grid doubling, max |d b_mn|/B_00 at s = ', s_list(j)
                call report(trim(tag), err, merge(1.0e-10_dp, 1.0e-8_dp, s_list(j) <= 0.5_dp))
            end do

            ! Tangent proportional to the data with random factors, independent
            ! per radial node. Hamada II uses second radial derivatives and
            ! divides by G' + iota I', so for it these rough tangents make the
            ! map strongly nonlinear (fourth-order truncation and round-off
            ! cross at 4e-8 on li383); there the factors are smooth in s,
            ! a (1 + b s) per harmonic.
            call wout_tangent_alloc(ctx, t)
            if (routes(ir) == SFL_STRAIGHT_CURRENT) then
                do i = 1, ctx%nstrm
                    fr = [urand(), urand(), urand(), urand(), urand(), urand()]
                    do j = 0, ctx%kpar
                        sj = real(j, dp)/real(ctx%kpar, dp)
                        t%rmnc(i, j) = fr(1)*(1.0_dp + fr(2)*sj)*ctx%hyb%rmnc(i, j)
                        t%zmns(i, j) = fr(3)*(1.0_dp + fr(4)*sj)*ctx%hyb%zmns(i, j)
                        t%lmns(i, j) = fr(5)*(1.0_dp + fr(6)*sj)*ctx%hyb%lmns(i, j)
                    end do
                end do
                fr(1:2) = [urand(), urand()]
                do j = 0, ctx%kpar
                    t%iotas(j) = fr(1)*(1.0_dp + fr(2)*real(j, dp)/real(ctx%kpar, dp)) &
                                 *ctx%hyb%iotas(j)
                end do
            else
                do j = 0, ctx%kpar
                    do i = 1, ctx%nstrm
                        t%rmnc(i, j) = urand()*ctx%hyb%rmnc(i, j)
                        t%zmns(i, j) = urand()*ctx%hyb%zmns(i, j)
                        t%lmns(i, j) = urand()*ctx%hyb%lmns(i, j)
                    end do
                    t%iotas(j) = urand()*ctx%hyb%iotas(j)
                end do
            end if
            t%phiedge = 0.3_dp*ctx%phiedge
            t%rmajor = 0.1_dp
            call boozer_ad_jvp(ctx, s_list, t, yp, yd)
            call spectrum_alloc(ctx, size(s_list), yfd)
            err = max(maxval(abs(yp%bmnc - y%bmnc)), maxval(abs(yp%numns - y%numns)))
            call report('primal of jvp = primal', err, 0.0_dp)
            ! fourth-order central differences, (8 (y(e) - y(-e)) - (y(2e) - y(-2e)))/(12 e)
            best = huge(1.0_dp)
            do k = 2, 10
                eps = 10.0_dp**(-k)
                cp = ctx
                cm = ctx
                call boozer_ad_perturb(cp, t, eps)
                call boozer_ad_perturb(cm, t, -eps)
                call boozer_ad_primal(cp, s_list, yp)
                call boozer_ad_primal(cm, s_list, ym)
                call spec_axpy(8.0_dp/12.0_dp, yp, ym, yfd)
                cp = ctx
                cm = ctx
                call boozer_ad_perturb(cp, t, 2*eps)
                call boozer_ad_perturb(cm, t, -2*eps)
                call boozer_ad_primal(cp, s_list, yp)
                call boozer_ad_primal(cm, s_list, ym)
                call spec_axpy(-1.0_dp/12.0_dp, yp, ym, yfd)
                err = sqrt(diff2(yfd, yd, eps)/spec_dot(yd, yd))
                print '(a, es8.1, a, es10.3)', '       eps = ', eps, ': JVP vs FD ', err
                best = min(best, err)
                if (k == 3) err3 = err
                if (k == 4) err4 = err
            end do
            if (routes(ir) == SFL_STRAIGHT_CURRENT) then
                ! The second radial derivatives of the quintic rows (weights of order
                ! (ns - 1)^2 on O(1) data) and the division by G' + iota I' (K/25 on
                ! li383) leave a relative round-off of about 2e-11 in the Hamada II
                ! spectrum, so the V-curve bottoms out near 3e-8 at eps = 1e-4,
                ! where the fourth-order truncation is still visible (err(1e-3) /
                ! err(1e-4) = 1e4). A JVP defect of relative size d would stop that
                ! decrease at d.
                call report('JVP vs fourth-order central differences (V-curve minimum)', &
                            best, 1.0e-7_dp)
                call report('fourth-order decrease err(1e-4)/err(1e-3)', err4/err3, 1.0e-3_dp)
            else
                call report('JVP vs fourth-order central differences (V-curve minimum)', &
                            best, 1.0e-8_dp)
            end if

            call spectrum_alloc(ctx, size(s_list), yb)
            do j = 1, size(s_list)
                do i = 1, ctx%mnboz
                    yb%bmnc(i, j) = urand()/b00(j)
                    yb%bmns(i, j) = urand()/b00(j)
                    yb%numns(i, j) = urand()
                    yb%numnc(i, j) = urand()
                    yb%gradsmnc(i, j) = urand()
                    yb%gradsmns(i, j) = urand()
                end do
                yb%iota(j) = urand()
                yb%g(j) = urand()/abs(y%g(j))
                yb%i(j) = urand()/abs(y%g(j))
                yb%psi_edge(j) = urand()
                yb%rmajor(j) = urand()
            end do
            call boozer_ad_vjp(ctx, s_list, yb, tb, yv)
            lhs = spec_dot(yb, yd)
            rhs = sum(tb%rmnc*t%rmnc) + sum(tb%zmns*t%zmns) + sum(tb%lmns*t%lmns) &
                  + sum(tb%iotas*t%iotas) + tb%phiedge*t%phiedge + tb%rmajor*t%rmajor
            ! Dot-product test normalized by |yb| |J t| (Cauchy-Schwarz), so that an
            ! accidental cancellation in <yb, J t> does not count as an error.
            call report('adjoint duality |<yb, J t> - <J^T yb, t>|/(|yb| |J t|)', &
                        abs(lhs - rhs)/sqrt(spec_dot(yb, yb)*spec_dot(yd, yd)), 1.0e-13_dp)
            err = max(maxval(abs(yv%bmnc - y%bmnc)), abs(yv%g(1) - y%g(1)))
            call report('primal of vjp = primal', err, 0.0_dp)
        end do
    end subroutine check_ad

    real(dp) function spec_dot(a, b)
        type(boozer_spectrum_t), intent(in) :: a, b

        spec_dot = sum(a%bmnc*b%bmnc) + sum(a%bmns*b%bmns) + sum(a%numns*b%numns) &
                   + sum(a%numnc*b%numnc) + sum(a%gradsmnc*b%gradsmnc) &
                   + sum(a%gradsmns*b%gradsmns) + sum(a%iota*b%iota) + sum(a%g*b%g) &
                   + sum(a%i*b%i) + sum(a%psi_edge*b%psi_edge) + sum(a%rmajor*b%rmajor)
    end function spec_dot

    !> c accumulates w (a - b) over all outputs.
    subroutine spec_axpy(w, a, b, c)
        real(dp), intent(in) :: w
        type(boozer_spectrum_t), intent(in) :: a, b
        type(boozer_spectrum_t), intent(inout) :: c

        if (w > 0.0_dp) call zero_spec(c)
        c%bmnc = c%bmnc + w*(a%bmnc - b%bmnc)
        c%bmns = c%bmns + w*(a%bmns - b%bmns)
        c%numns = c%numns + w*(a%numns - b%numns)
        c%numnc = c%numnc + w*(a%numnc - b%numnc)
        c%gradsmnc = c%gradsmnc + w*(a%gradsmnc - b%gradsmnc)
        c%gradsmns = c%gradsmns + w*(a%gradsmns - b%gradsmns)
        c%iota = c%iota + w*(a%iota - b%iota)
        c%g = c%g + w*(a%g - b%g)
        c%i = c%i + w*(a%i - b%i)
        c%psi_edge = c%psi_edge + w*(a%psi_edge - b%psi_edge)
        c%rmajor = c%rmajor + w*(a%rmajor - b%rmajor)
    end subroutine spec_axpy

    subroutine zero_spec(c)
        type(boozer_spectrum_t), intent(inout) :: c

        c%bmnc = 0.0_dp
        c%bmns = 0.0_dp
        c%numns = 0.0_dp
        c%numnc = 0.0_dp
        c%gradsmnc = 0.0_dp
        c%gradsmns = 0.0_dp
        c%iota = 0.0_dp
        c%g = 0.0_dp
        c%i = 0.0_dp
        c%psi_edge = 0.0_dp
        c%rmajor = 0.0_dp
    end subroutine zero_spec

    !> |yfd/eps - yd|^2 over all outputs (yfd the difference stencil times eps).
    real(dp) function diff2(yfd, yd, eps)
        type(boozer_spectrum_t), intent(in) :: yfd, yd
        real(dp), intent(in) :: eps

        diff2 = sum((yfd%bmnc/eps - yd%bmnc)**2) + sum((yfd%bmns/eps - yd%bmns)**2) &
                + sum((yfd%numns/eps - yd%numns)**2) + sum((yfd%numnc/eps - yd%numnc)**2) &
                + sum((yfd%gradsmnc/eps - yd%gradsmnc)**2) &
                + sum((yfd%gradsmns/eps - yd%gradsmns)**2) &
                + sum((yfd%iota/eps - yd%iota)**2) + sum((yfd%g/eps - yd%g)**2) &
                + sum((yfd%i/eps - yd%i)**2) + sum((yfd%psi_edge/eps - yd%psi_edge)**2) &
                + sum((yfd%rmajor/eps - yd%rmajor)**2)
    end function diff2

end program test_sfl_routes
