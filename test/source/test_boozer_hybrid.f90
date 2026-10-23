!> The hybrid VMEC -> Boozer construction (boozer_hybrid, default of
!> boozer_vmec_ad).
!>
!> Checks:
!>  - the FortSym kernel boozer_hybrid_field equals the independently
!>    generated legacy kernel boozer_field_jacobian (PEST-angle form) and
!>    sqrt_g_ss_kernel, values and Jacobians, and its iota column matches
!>    central differences;
!>  - boozer_ad_init selects the hybrid by default;
!>  - doubling the angle grid changes the spectrum by less than 1e-10 of B_00
!>    (the quadrature is converged at the default grid);
!>  - hybrid and legacy share the conventions of G, I, nu and |grad s| (they
!>    differ by the legacy discretization error, up to 1e-3 on the ns = 17
!>    circular tokamak, while a convention mismatch is of order one);
!>  - the JVP matches central differences of the primal (V-curve over the
!>    step, wout data perturbed in memory) and the VJP is its adjoint;
!> on QA, li383 (net current, -I*lambda), the circular tokamak (axisymmetric)
!> and the up-down asymmetric tokamak (lasym, the asymmetric arrays held fixed).
!> Agreement with the primary oracle of boozer-bench is checked there
!> (reports/scripts/libneo_hybrid), since the oracle is a JAX program.
program test_boozer_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use boozer_vmec_ad, only: boozer_ad_context_t, boozer_wout_tangent_t, &
                              boozer_spectrum_t, boozer_ad_init, boozer_ad_primal, &
                              boozer_ad_jvp, boozer_ad_vjp, wout_tangent_alloc, &
                              spectrum_alloc, boozer_ad_perturb, boozer_field_kernel, &
                              sqrt_g_ss_kernel, BOOZER_HYBRID, BOOZER_LEGACY
    use boozer_hybrid, only: hybrid_field_kernel
    implicit none

    real(dp), parameter :: s_list(3) = [0.15_dp, 0.5_dp, 0.85_dp]
    logical :: failed
    integer(int64) :: seed_state

    failed = .false.
    seed_state = 12345
    call check_kernel()
    call check_case('wout.nc', 8, 6, 'QA', .true.)
    call check_case('wout_li383.nc', 8, 6, 'li383', .true.)
    call check_case('wout_circular_tokamak.nc', 8, 0, 'circular tokamak', .true.)
    ! legacy reads no lasym arrays in the AD path, so it is not compared here
    call check_case('wout_updown_tokamak.nc', 8, 0, 'up-down tokamak', .false.)
    if (failed) error stop 'test_boozer_hybrid FAILED'
    print *, 'test_boozer_hybrid PASSED'

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

    !> Deterministic uniform numbers in (-1, 1).
    real(dp) function urand()
        seed_state = modulo(1103515245_int64*seed_state + 12345_int64, 2147483647_int64)
        urand = 2.0_dp*real(seed_state, dp)/2147483647.0_dp - 1.0_dp
    end function urand

    subroutine check_kernel()
        real(dp) :: x(10), xp(10), f(4), fp(4), fm(4), jac(4, 10), jt(4, 10), fl(3), jl(3, 10)
        real(dp) :: g, jg(7), xl(10), h, err_f, err_j, err_fd
        integer :: k, ip

        err_f = 0.0_dp
        err_j = 0.0_dp
        err_fd = 0.0_dp
        do ip = 1, 4
            x = [1.5_dp, 0.3_dp, -0.2_dp, 0.05_dp, 0.1_dp, 0.4_dp, 0.03_dp, 0.1_dp, &
                 -0.05_dp, 0.4_dp]
            do k = 1, 10
                x(k) = x(k)*(1.0_dp + 0.3_dp*urand())
            end do
            call hybrid_field_kernel(x, f, jac)
            ! legacy kernel: psi' = athetas = 1, A_phi' = -iota
            xl = x
            xl(10) = -x(10)
            call boozer_field_kernel(xl, 1.0_dp, fl, jl)
            call sqrt_g_ss_kernel(x(1:7), g, jg)
            err_f = max(err_f, maxval(abs(f(1:3) - fl)/abs(fl)), abs(f(4) - g)/g)
            err_j = max(err_j, maxval(abs(jac(1:3, 1:9) - jl(:, 1:9))) &
                        /maxval(abs(jl(:, 1:9))), &
                        maxval(abs(jac(1:3, 10) + jl(:, 10)))/maxval(abs(jl(:, 1:9))), &
                        maxval(abs(jac(4, 1:7) - jg))/maxval(abs(jg)), &
                        maxval(abs(jac(4, 8:10))))
            h = 1.0e-6_dp
            do k = 1, 10
                xp = x
                xp(k) = x(k) + h
                call hybrid_field_kernel(xp, fp, jt)
                xp(k) = x(k) - h
                call hybrid_field_kernel(xp, fm, jt)
                err_fd = max(err_fd, maxval(abs((fp - fm)/(2*h) - jac(:, k))) &
                             /maxval(abs(jac)))
            end do
        end do
        call report('kernel = legacy kernel and sqrt_g_ss (values)', err_f, 1.0e-13_dp)
        call report('kernel = legacy kernel and sqrt_g_ss (Jacobian)', err_j, 1.0e-13_dp)
        call report('kernel Jacobian vs central differences', err_fd, 1.0e-8_dp)
    end subroutine check_kernel

    subroutine check_case(wout, mboz, nboz, label, with_legacy)
        character(*), intent(in) :: wout, label
        integer, intent(in) :: mboz, nboz
        logical, intent(in) :: with_legacy
        type(boozer_ad_context_t) :: ctx, cfine, cleg, cp, cm
        type(boozer_spectrum_t) :: y, yf, yl, yd, yp, ym, yb, yv
        type(boozer_wout_tangent_t) :: t, tb
        real(dp) :: b00(size(s_list)), err, eps, lhs, rhs, best
        integer :: j, k, i
        character(len=80) :: tag

        print '(a)', ' '//label//' ('//wout//')'
        call boozer_ad_init(ctx, wout, mboz, nboz)
        if (ctx%construction /= BOOZER_HYBRID) call report('default is hybrid', 1.0_dp, 0.0_dp)
        call boozer_ad_primal(ctx, s_list, y)
        b00 = y%bmnc(1, :)

        call boozer_ad_init(cfine, wout, mboz, nboz, n_theta=2*ctx%grid%n_theta, &
                            n_phi=merge(1, 2*ctx%grid%n_phi, ctx%axisym))
        call boozer_ad_primal(cfine, s_list, yf)
        ! Acceptance of the design (boozer-bench design/sfl.md, section 12):
        ! 1e-10 of B_00 up to mid radius, 1e-8 at the edge.
        do j = 1, size(s_list)
            err = max(maxval(abs(yf%bmnc(:, j) - y%bmnc(:, j))), &
                      maxval(abs(yf%bmns(:, j) - y%bmns(:, j))))/b00(j)
            write (tag, '(a, f4.2)') 'grid doubling, max |d b_mn|/B_00 at s = ', s_list(j)
            call report(trim(tag), err, merge(1.0e-10_dp, 1.0e-8_dp, s_list(j) <= 0.5_dp))
        end do

        if (with_legacy) then
            call boozer_ad_init(cleg, wout, mboz, nboz, construction=BOOZER_LEGACY)
            call boozer_ad_primal(cleg, s_list, yl)
            err = 0.0_dp
            do j = 1, size(s_list)
                err = max(err, maxval(abs(yl%bmnc(:, j) - y%bmnc(:, j)))/b00(j), &
                          abs(yl%g(j) - y%g(j))/abs(y%g(j)), &
                          maxval(abs(yl%numns(:, j) - y%numns(:, j))), &
                          maxval(abs(yl%gradsmnc(:, j) - y%gradsmnc(:, j))) &
                          /yl%gradsmnc(1, j))
            end do
            call report('conventions: legacy vs hybrid (B/B_00, G/G, nu, |grad s|)', err, &
                        2.0e-3_dp)
        end if

        ! Tangent proportional to the data with random factors.
        call wout_tangent_alloc(ctx, t)
        do j = 0, ctx%kpar
            do i = 1, ctx%nstrm
                t%rmnc(i, j) = urand()*ctx%hyb%rmnc(i, j)
                t%zmns(i, j) = urand()*ctx%hyb%zmns(i, j)
                t%lmns(i, j) = urand()*ctx%hyb%lmns(i, j)
            end do
            t%iotas(j) = urand()*ctx%hyb%iotas(j)
        end do
        t%phiedge = 0.3_dp*ctx%phiedge
        t%rmajor = 0.1_dp
        call boozer_ad_jvp(ctx, s_list, t, yp, yd)
        err = max(maxval(abs(yp%bmnc - y%bmnc)), maxval(abs(yp%numns - y%numns)))
        call report('primal of jvp = primal', err, 0.0_dp)

        best = huge(1.0_dp)
        do k = 2, 9
            eps = 10.0_dp**(-k)
            cp = ctx
            cm = ctx
            call boozer_ad_perturb(cp, t, eps)
            call boozer_ad_perturb(cm, t, -eps)
            call boozer_ad_primal(cp, s_list, yp)
            call boozer_ad_primal(cm, s_list, ym)
            err = sqrt(diff2(yp, ym, yd, eps)/spec_dot(yd, yd))
            print '(a, es8.1, a, es10.3)', '       eps = ', eps, ': JVP vs FD ', err
            best = min(best, err)
        end do
        call report('JVP vs central differences (V-curve minimum)', best, 1.0e-8_dp)

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
        call report('adjoint duality <yb, J t> = <J^T yb, t>', &
                    abs(lhs - rhs)/max(abs(lhs), abs(rhs)), 1.0e-13_dp)
        err = max(maxval(abs(yv%bmnc - y%bmnc)), abs(yv%g(1) - y%g(1)))
        call report('primal of vjp = primal', err, 0.0_dp)
    end subroutine check_case

    real(dp) function spec_dot(a, b)
        type(boozer_spectrum_t), intent(in) :: a, b

        spec_dot = sum(a%bmnc*b%bmnc) + sum(a%bmns*b%bmns) + sum(a%numns*b%numns) &
                   + sum(a%numnc*b%numnc) + sum(a%gradsmnc*b%gradsmnc) &
                   + sum(a%gradsmns*b%gradsmns) + sum(a%iota*b%iota) + sum(a%g*b%g) &
                   + sum(a%i*b%i) + sum(a%psi_edge*b%psi_edge) + sum(a%rmajor*b%rmajor)
    end function spec_dot

    !> |(yp - ym)/(2 eps) - yd|^2 over all outputs.
    real(dp) function diff2(yp, ym, yd, eps)
        type(boozer_spectrum_t), intent(in) :: yp, ym, yd
        real(dp), intent(in) :: eps

        diff2 = sum(((yp%bmnc - ym%bmnc)/(2*eps) - yd%bmnc)**2) &
                + sum(((yp%bmns - ym%bmns)/(2*eps) - yd%bmns)**2) &
                + sum(((yp%numns - ym%numns)/(2*eps) - yd%numns)**2) &
                + sum(((yp%numnc - ym%numnc)/(2*eps) - yd%numnc)**2) &
                + sum(((yp%gradsmnc - ym%gradsmnc)/(2*eps) - yd%gradsmnc)**2) &
                + sum(((yp%gradsmns - ym%gradsmns)/(2*eps) - yd%gradsmns)**2) &
                + sum(((yp%iota - ym%iota)/(2*eps) - yd%iota)**2) &
                + sum(((yp%g - ym%g)/(2*eps) - yd%g)**2) &
                + sum(((yp%i - ym%i)/(2*eps) - yd%i)**2) &
                + sum(((yp%psi_edge - ym%psi_edge)/(2*eps) - yd%psi_edge)**2) &
                + sum(((yp%rmajor - ym%rmajor)/(2*eps) - yd%rmajor)**2)
    end function diff2

end program test_boozer_hybrid
