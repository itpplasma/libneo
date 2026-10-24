!> Derivatives of the VMEC -> Boozer spectrum (boozer_vmec_ad).
!>
!> Checks, on the reference QA wout with the linear and the cubic lambda
!> half -> full mesh map and on li383 (net toroidal current, so the -I*lambda
!> term of the angle potential matters) with the cubic map:
!>  - the FortSym field kernel reproduces vmec_field and its Jacobian matches
!>    central differences of the kernel;
!>  - the surface core JVP matches central differences of the primal
!>    boozer_surface_transform, and core VJP/JVP satisfy adjoint duality;
!>  - the full JVP matches central differences of the complete libneo
!>    pipeline (wout harmonics perturbed in memory, VMEC data re-splined), and
!>    the full VJP is its adjoint.
program test_boozer_vmec_ad
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_surface, only: boozer_angle_grid_t, boozer_surface_lin_t, &
                              boozer_surface_lin_init, boozer_surface_transform, &
                              boozer_surface_jvp, boozer_surface_vjp
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, boozer_ad_primal, &
                              boozer_ad_jvp, boozer_ad_vjp, wout_tangent_alloc, &
                              spectrum_alloc, boozer_field_kernel
    implicit none

    type(boozer_ad_context_t) :: ctx
    real(dp), parameter :: s_list(3) = [0.15_dp, 0.5_dp, 0.85_dp]
    logical :: failed

    failed = .false.
    call boozer_ad_init(ctx, 'wout.nc', 8, 6, grid_refinement=3)
    call check_kernel()
    call check_core()
    call check_full('QA, linear lambda')
    call check_case('wout.nc', 'cubic', 'QA, cubic lambda')
    call check_case('wout_li383.nc', 'cubic', 'li383, cubic lambda')
    if (failed) error stop 'test_boozer_vmec_ad FAILED'
    print *, 'test_boozer_vmec_ad PASSED'

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

    subroutine check_kernel()
        use spline_vmec_sub, only: splint_vmec_data, vmec_field
        real(dp) :: x(10), xp(10), f(3), fp(3), fm(3), jac(3, 10), jd(3, 10)
        real(dp) :: A_phi, A_theta, aphis, athetas, aiota, R, Z, alam, ls, h
        real(dp) :: sqg, lt, lp, bct, bcp, bcs, bcvt, bcvp, err_f, err_j
        integer :: k, ip

        err_f = 0.0_dp
        err_j = 0.0_dp
        do ip = 1, 3
            call splint_vmec_data(s_list(ip), 0.3_dp*ip, 0.2_dp*ip, A_phi, A_theta, &
                                  aphis, athetas, aiota, R, Z, alam, x(2), x(3), &
                                  x(4), x(5), x(6), x(7), ls, x(8), x(9))
            x(1) = R
            x(10) = aphis
            call boozer_field_kernel(x, athetas, f, jac)
            call vmec_field(s_list(ip), 0.3_dp*ip, 0.2_dp*ip, A_theta, A_phi, &
                            athetas, aphis, aiota, sqg, alam, ls, lt, lp, bct, bcp, &
                            bcs, bcvt, bcvp)
            err_f = max(err_f, abs(f(1) - sqrt(bct*bcvt + bcp*bcvp))/abs(f(1)), &
                        abs(f(2) - bcvt*(1 + lt))/abs(f(2)), &
                        abs(f(3) - (bcvp + bcvt*lp))/abs(f(3)))
            do k = 1, 10
                h = 1.0e-6_dp*max(abs(x(k)), 1.0e-3_dp)
                xp = x
                xp(k) = x(k) + h
                call boozer_field_kernel(xp, athetas, fp, jd(:, 1))
                xp(k) = x(k) - h
                call boozer_field_kernel(xp, athetas, fm, jd(:, 1))
                jd(:, k) = (fp - fm)/(2*h)
            end do
            do k = 1, 3
                err_j = max(err_j, maxval(abs(jd(k, :) - jac(k, :))) &
                            /maxval(abs(jac(k, :))))
            end do
        end do
        call report('kernel value vs vmec_field', err_f, 1.0e-12_dp)
        call report('kernel Jacobian vs FD', err_j, 1.0e-6_dp)
    end subroutine check_kernel

    subroutine check_core()
        type(boozer_angle_grid_t) :: grid
        type(boozer_surface_lin_t) :: lin
        integer, parameter :: nt = 25, np = 17
        real(dp) :: lam(nt, np), bt(nt, np), bp(nt, np), car(1, nt, np), iota
        real(dp) :: lamd(nt, np), btd(nt, np), bpd(nt, np), card(1, nt, np), iotad
        real(dp) :: dthv(nt, np), dphv(nt, np), q2(3, nt, np), bct, bcp
        real(dp) :: dthv_d(nt, np), dphv_d(nt, np), q2_d(3, nt, np), bct_d, bcp_d
        real(dp) :: ybt, ybp, yth(nt, np), yph(nt, np), yq(3, nt, np)
        real(dp) :: xb_iota, xb_lam(nt, np), xb_bt(nt, np), xb_bp(nt, np)
        real(dp) :: xb_car(1, nt, np), lhs, rhs, eps, err
        real(dp) :: qp(3, nt, np), qm(3, nt, np), tp(nt, np), tm(nt, np)
        real(dp) :: pp(nt, np), pm(nt, np), ip_, im_, gp, gm
        integer :: a, b

        grid = boozer_angle_grid_t(n_theta=nt, n_phi=np, ns_tp=5, &
                                   h_theta=8*atan(1.0_dp)/(nt - 1), &
                                   h_phi=8*atan(1.0_dp)/(2*(np - 1)))
        call boozer_surface_lin_init(grid, lin)
        ! The phase 0.3 breaks stellarator symmetry on purpose: at symmetric points
        ! Delta theta = 0 puts a Boozer node exactly on a target, the stencil is a
        ! tie and the primal has a kink there (one-sided derivatives differ by the
        ! interpolation error), which a central difference cannot reproduce.
        iota = 0.42_dp
        do b = 1, np
            do a = 1, nt
                associate (t => (a - 1)*grid%h_theta, p => 2*(b - 1)*grid%h_phi)
                    lam(a, b) = 0.1_dp*sin(t + 0.3_dp) + 0.03_dp*sin(2*t - p)
                    bt(a, b) = 0.2_dp + 0.05_dp*cos(t) + 0.02_dp*cos(t - p)
                    bp(a, b) = 3.0_dp + 0.1_dp*cos(t) + 0.04_dp*cos(2*t - p)
                    car(1, a, b) = 1.0_dp + 0.1_dp*cos(t) + 0.05_dp*cos(t - p)
                    lamd(a, b) = cos(3*t + p) + 0.3_dp*sin(t)
                    btd(a, b) = sin(2*t) + cos(p)
                    bpd(a, b) = cos(t - p)
                    card(1, a, b) = sin(t + 2*p)
                    yth(a, b) = cos(t*1.3_dp)
                    yph(a, b) = sin(p + t)
                    yq(:, a, b) = [sin(t), cos(p), cos(t + p)]
                end associate
            end do
        end do
        iotad = 0.7_dp
        call boozer_surface_jvp(lin, iota, lam, bt, bp, car, iotad, lamd, btd, bpd, &
                                card, bct, bcp, dthv, dphv, q2, bct_d, bcp_d, &
                                dthv_d, dphv_d, q2_d)

        eps = 1.0e-6_dp
        call core_primal(grid, eps, iota, lam, bt, bp, car, iotad, lamd, btd, bpd, &
                         card, qp, tp, pp, ip_, gp)
        call core_primal(grid, -eps, iota, lam, bt, bp, car, iotad, lamd, btd, bpd, &
                         card, qm, tm, pm, im_, gm)
        err = max(maxval(abs((qp - qm)/(2*eps) - q2_d)), &
                  maxval(abs((tp - tm)/(2*eps) - dthv_d)), &
                  maxval(abs((pp - pm)/(2*eps) - dphv_d)), &
                  abs((ip_ - im_)/(2*eps) - bct_d), abs((gp - gm)/(2*eps) - bcp_d))
        call report('core JVP vs FD of primal', err, 1.0e-7_dp)

        ybt = 0.3_dp
        ybp = -1.1_dp
        xb_iota = 0.0_dp
        xb_lam = 0.0_dp
        xb_bt = 0.0_dp
        xb_bp = 0.0_dp
        xb_car = 0.0_dp
        call boozer_surface_vjp(lin, iota, lam, bt, bp, car, ybt, ybp, yth, yph, yq, &
                                xb_iota, xb_lam, xb_bt, xb_bp, xb_car)
        lhs = ybt*bct_d + ybp*bcp_d + sum(yth*dthv_d) + sum(yph*dphv_d) + sum(yq*q2_d)
        rhs = xb_iota*iotad + sum(xb_lam*lamd) + sum(xb_bt*btd) + sum(xb_bp*bpd) &
              + sum(xb_car*card)
        call report('core adjoint duality', abs(lhs - rhs)/abs(lhs), 1.0e-13_dp)

    end subroutine check_core

    subroutine core_primal(grid, e, iota, lam, bt, bp, car, iotad, lamd, btd, bpd, &
                           card, q, dth, dph, gi, gg)
        type(boozer_angle_grid_t), intent(in) :: grid
        real(dp), intent(in) :: e, iota, lam(:, :), bt(:, :), bp(:, :), car(:, :, :)
        real(dp), intent(in) :: iotad, lamd(:, :), btd(:, :), bpd(:, :), card(:, :, :)
        real(dp), intent(out) :: q(:, :, :), dth(:, :), dph(:, :), gi, gg

        q = 0.0_dp
        q(3, :, :) = car(1, :, :) + e*card(1, :, :)
        call boozer_surface_transform(grid, iota + e*iotad, lam + e*lamd, &
                                      bt + e*btd, bp + e*bpd, q, gi, gg, dth, dph)
    end subroutine core_primal

    subroutine check_case(wout, lambda_map, label)
        use new_vmec_stuff_mod, only: lambda_half_to_full
        character(*), intent(in) :: wout, lambda_map, label

        lambda_half_to_full = lambda_map
        call boozer_ad_init(ctx, wout, 8, 6, grid_refinement=3)
        call check_full(label)
        lambda_half_to_full = 'linear'
    end subroutine check_case

    subroutine check_full(label)
        character(*), intent(in) :: label
        type(boozer_wout_tangent_t) :: t, tb
        type(boozer_spectrum_t) :: y0, yd, yb, yp, ym
        real(dp) :: lhs, rhs, eps, scale, err
        integer :: i, j

        call wout_tangent_alloc(ctx, t)
        do j = 0, ctx%kpar
            associate (s => real(j, dp)/ctx%kpar)
                do i = 1, ctx%nstrm
                    t%rmnc(i, j) = 1.0e-3_dp*s*cos(0.7_dp*i)/(1 + ctx%xm(i))**2
                    t%zmns(i, j) = 1.0e-3_dp*s*sin(1.3_dp*i)/(1 + ctx%xm(i))**2
                    t%lmns(i, j) = 1.0e-2_dp*sin(4*s)*cos(2.1_dp*i)/(1 + ctx%xm(i))**2
                end do
                t%iotaf(j) = 0.05_dp*(1 - s**2)
            end associate
        end do
        call boozer_ad_jvp(ctx, s_list, t, y0, yd)

        call spectrum_alloc(ctx, size(s_list), yb)
        do j = 1, size(s_list)
            do i = 1, ctx%mnboz
                yb%bmnc(i, j) = cos(0.37_dp*i*j)
                yb%bmns(i, j) = 0.1_dp*sin(0.5_dp*i + j)
                yb%numns(i, j) = sin(0.21_dp*i*j)
                yb%numnc(i, j) = 0.2_dp*cos(0.9_dp*i)
            end do
            yb%iota(j) = 0.3_dp*j
            yb%g(j) = -0.2_dp*j
            yb%i(j) = 1.7_dp
        end do
        call boozer_ad_vjp(ctx, s_list, yb, tb)
        lhs = sum(yb%bmnc*yd%bmnc) + sum(yb%bmns*yd%bmns) + sum(yb%numns*yd%numns) &
              + sum(yb%numnc*yd%numnc) + sum(yb%iota*yd%iota) + sum(yb%g*yd%g) &
              + sum(yb%i*yd%i)
        rhs = sum(tb%rmnc*t%rmnc) + sum(tb%zmns*t%zmns) + sum(tb%lmns*t%lmns) &
              + sum(tb%iotaf*t%iotaf)
        print '(a)', '  '//label
        call report('full adjoint duality', abs(lhs - rhs)/abs(lhs), 1.0e-12_dp)

        eps = 1.0e-4_dp
        call respline(t, eps)
        call boozer_ad_primal(ctx, s_list, yp)
        call respline(t, -eps)
        call boozer_ad_primal(ctx, s_list, ym)
        call respline(t, 0.0_dp)
        scale = maxval(abs(yd%bmnc))
        err = maxval(abs((yp%bmnc - ym%bmnc)/(2*eps) - yd%bmnc))/scale
        call report('full JVP bmnc vs FD', err, 1.0e-5_dp)
        err = maxval(abs((yp%numns - ym%numns)/(2*eps) - yd%numns))/maxval(abs(yd%numns))
        call report('full JVP numns vs FD', err, 1.0e-5_dp)
        err = maxval(abs((yp%iota - ym%iota)/(2*eps) - yd%iota))/maxval(abs(yd%iota))
        call report('full JVP iota vs FD', err, 1.0e-7_dp)
        err = max(maxval(abs((yp%g - ym%g)/(2*eps) - yd%g))/maxval(abs(yd%g)), &
                  maxval(abs((yp%i - ym%i)/(2*eps) - yd%i))/maxval(abs(yd%i)))
        call report('full JVP G, I vs FD', err, 1.0e-6_dp)
        print '(a, 3es12.4)', '  primal iota, G, I (s=0.5): ', y0%iota(2), y0%g(2), y0%i(2)
    end subroutine check_full

    !> Re-spline the VMEC data with the wout harmonics shifted by eps*t, exactly
    !> as spline_vmec_data does after reading the file.
    subroutine respline(t, eps)
        use new_vmec_stuff_mod, only: rmnc, zmns, almns, aiota, kpar
        use vmec_alloc_sub, only: new_allocate_vmec_stuff, new_deallocate_vmec_stuff
        use vmecin_sub, only: lambda_half_to_full_map
        use spline_vmec_sub, only: initialize_vmec_data_and_arrays, &
                                   setup_poloidal_flux_splines, &
                                   setup_angular_grid_and_fourier_synthesis, &
                                   spline_over_phi, spline_over_theta, spline_over_s
        type(boozer_wout_tangent_t), intent(in) :: t
        real(dp), intent(in) :: eps
        real(dp), allocatable :: dlam(:, :)

        call new_allocate_vmec_stuff()
        call initialize_vmec_data_and_arrays()
        rmnc = rmnc + eps*ctx%fac_r*t%rmnc
        zmns = zmns + eps*ctx%fac_r*t%zmns
        allocate (dlam(size(t%lmns, 1), 0:kpar))
        call lambda_half_to_full_map(t%lmns, dlam, kpar)
        almns = almns + eps*dlam
        aiota = aiota + eps*t%iotaf
        call setup_poloidal_flux_splines()
        call setup_angular_grid_and_fourier_synthesis()
        call new_deallocate_vmec_stuff()
        call spline_over_phi()
        call spline_over_theta()
        call spline_over_s()
    end subroutine respline

end program test_boozer_vmec_ad
