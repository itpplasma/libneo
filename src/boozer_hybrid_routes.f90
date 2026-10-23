!> Straight-field-line routes of the hybrid construction beyond the covariant
!> Boozer route, with exact tangent and adjoint, and their diagnostics.
!>
!> Every route writes the angles with one potential nu on the VMEC grid,
!> theta* = theta + lambda + iota nu, phi* = phi + nu (Boozer 1981;
!> D'haeseleer et al. 1991, ch. 6), and differs only in the condition on nu.
!>
!> Jacobian routes (SFL_JACOBIAN_B2: J* B^2 flux function, Boozer;
!> SFL_JACOBIAN_CONST: J* flux function, Hamada I). With the PEST angle
!> theta_P = theta + lambda and J_P = sqrt(g)_V/(1 + lambda_t), J* = J_P/(1 +
!> (iota d_thetaP + d_phi) nu), so the condition is the magnetic differential
!> equation
!>     (iota d_thetaP + d_phi) nu = q/<q>_P - 1,   q = J_P B^2 or J_P,
!> solved mode by mode, nu_mn = f_mn/(i (iota m - n N)), on 0 <= m <= mnu,
!> |n| <= nnu. The PEST analysis f_mn = 2 <(1 + lambda_t) (q/<q>_P - 1)
!> exp(-i(m theta_P - n N phi))> is a quadrature on the uniform VMEC grid
!> (exact change of variables, no inversion), and nu, d_thetaP nu and d_phi nu
!> are synthesized at theta_P(theta, phi). Modes with |iota m - n N| <
!> res_exact max(1, m) are set to zero, as GVEC's sfl_hamada.F90 does; the
!> tangent of nu with respect to iota follows from L(iota) nu = f (implicit
!> function theorem) with the zeroed set held fixed.
!>
!> Hamada II (SFL_STRAIGHT_CURRENT, gvec#95): from the covariant Boozer
!> potential nu_B, with B = I grad theta_B + G grad phi_B + beta grad s,
!>     nu = nu_B + eta,  eta = -(beta - <beta>_B)/(G' + iota I'),
!> which makes the current potential I' theta_H + G' phi_H linear in the
!> angles. beta = B_s + B_theta d_s theta|_B + B_phi d_s phi|_B (kernel
!> sfl_beta), with nu_s = (w_s - I' lambda - I lambda_s - K' nu_B)/K from the
!> radial derivatives of the covariant components at fixed VMEC angles
!> (kernel sfl_radial_field, second radial derivatives of the quintic rows),
!> so no radial finite difference enters. eta is truncated to the nu modes
!> (|m| <= mnu, |n| <= nnu in VMEC angles) on the grid, and eta_theta, eta_phi
!> are the exact derivatives of the truncated series. Since
!> G' + iota I' = -mu0 p' <sqrt g>/psi' under averaged force balance, the
!> route needs p' /= 0 and stops with an error otherwise.
submodule(boozer_hybrid) boozer_hybrid_routes
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

contains

    include 'generated/sfl_route_kernels_leaves.inc'

    !> sqrt(g)_V and sqrt(g)_V |B|^2 with their Jacobian with respect to
    !> (R, R_s, R_t, Z_s, Z_t, |B|).
    module pure subroutine route_weight_kernel(x, f, jac)
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: f(2), jac(2, 6)

        call sfl_route_weight(x(1), x(2), x(3), x(4), x(5), x(6), &
                              f(1), jac(1, 1), jac(1, 2), jac(1, 3), jac(1, 4), jac(1, 5), &
                              jac(1, 6), &
                              f(2), jac(2, 1), jac(2, 2), jac(2, 3), jac(2, 4), jac(2, 5), &
                              jac(2, 6))
    end subroutine route_weight_kernel

    !> (B_s, d_s B_theta, d_s B_phi) at psi' = 1 and the Jacobian with respect to
    !> (x(1:10), R_ss, R_st, R_sp, Z_ss, Z_st, Z_sp, lambda_st, lambda_sp, iota').
    module pure subroutine route_radial_kernel(x, f, jac)
        real(dp), intent(in) :: x(NRIN)
        real(dp), intent(out) :: f(3), jac(3, NRIN)

        call sfl_radial_field(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), &
                              x(11), x(12), x(13), x(14), x(15), x(16), x(17), x(18), x(19), &
                              f(1), jac(1, 1), jac(1, 2), jac(1, 3), jac(1, 4), jac(1, 5), &
                              jac(1, 6), jac(1, 7), jac(1, 8), jac(1, 9), jac(1, 10), &
                              jac(1, 11), jac(1, 12), jac(1, 13), jac(1, 14), jac(1, 15), &
                              jac(1, 16), jac(1, 17), jac(1, 18), jac(1, 19), f(2), jac(2, 1), &
                              jac(2, 2), jac(2, 3), jac(2, 4), jac(2, 5), jac(2, 6), jac(2, 7), &
                              jac(2, 8), jac(2, 9), jac(2, 10), jac(2, 11), jac(2, 12), &
                              jac(2, 13), jac(2, 14), jac(2, 15), jac(2, 16), jac(2, 17), &
                              jac(2, 18), jac(2, 19), f(3), jac(3, 1), jac(3, 2), jac(3, 3), &
                              jac(3, 4), jac(3, 5), jac(3, 6), jac(3, 7), jac(3, 8), jac(3, 9), &
                              jac(3, 10), jac(3, 11), jac(3, 12), jac(3, 13), jac(3, 14), &
                              jac(3, 15), jac(3, 16), jac(3, 17), jac(3, 18), jac(3, 19))
    end subroutine route_radial_kernel

    !> beta and its gradient with respect to (B_s, B_theta, B_phi, lambda_t,
    !> lambda_p, lambda_s, nu, nu_t, nu_p, nu_s, iota, iota').
    module pure subroutine route_beta_kernel(x, f, jac)
        real(dp), intent(in) :: x(NBIN)
        real(dp), intent(out) :: f, jac(NBIN)

        call sfl_beta(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), &
                      x(12), f, jac(1), jac(2), jac(3), jac(4), jac(5), jac(6), jac(7), &
                      jac(8), jac(9), jac(10), jac(11), jac(12))
    end subroutine route_beta_kernel

    !> VMEC angles of the grid points (theta fastest).
    function grid_theta(h) result(t)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp) :: t(h%nt*h%np)

        t = reshape(spread(h%theta, 2, h%np), [h%nt*h%np])
    end function grid_theta

    function grid_phi(h) result(t)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp) :: t(h%nt*h%np)

        t = reshape(spread(h%phi, 1, h%nt), [h%nt*h%np])
    end function grid_phi

    !> Exact angular derivatives of the trigonometric interpolant of f on the
    !> uniform VMEC grid (numpy fftfreq ordering, as potential).
    subroutine grid_deriv(h, f, ft, fp)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: f(:)
        real(dp), intent(out) :: ft(:), fp(:)
        complex(dp) :: y(h%nt, h%np)

        y = dft2(h, cmplx(reshape(f, [h%nt, h%np]), 0.0_dp, kind=dp), .false.)
        ft = reshape(real(dft2(h, IU*spread(h%km, 2, h%np)*y, .true.)), [h%nt*h%np])
        fp = reshape(real(dft2(h, IU*spread(h%kn, 1, h%nt)*y, .true.)), [h%nt*h%np])
    end subroutine grid_deriv

    !> Truncation of f on the VMEC grid to |m| <= mnu, |n| <= nnu (field-period
    !> harmonics) and the exact angular derivatives of the truncated function:
    !> ff, ft, fp = Re IDFT(M (1, i k_m, i k_n) DFT f) with the mode mask M.
    subroutine grid_filter(h, f, ff, ft, fp)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: f(:)
        real(dp), intent(out) :: ff(:), ft(:), fp(:)
        complex(dp) :: y(h%nt, h%np)

        y = filter_mask(h)*dft2(h, cmplx(reshape(f, [h%nt, h%np]), 0.0_dp, kind=dp), .false.)
        ff = reshape(real(dft2(h, y, .true.)), [h%nt*h%np])
        ft = reshape(real(dft2(h, IU*spread(h%km, 2, h%np)*y, .true.)), [h%nt*h%np])
        fp = reshape(real(dft2(h, IU*spread(h%kn, 1, h%nt)*y, .true.)), [h%nt*h%np])
    end subroutine grid_filter

    !> Transpose of grid_filter, accumulated into fb.
    subroutine grid_filter_t(h, ffb, ftb, fpb, fb)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: ffb(:), ftb(:), fpb(:)
        real(dp), intent(inout) :: fb(:)
        complex(dp) :: v(h%nt, h%np)

        v = dft2(h, cmplx(reshape(ffb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.) &
            + IU*spread(h%km, 2, h%np) &
            *dft2(h, cmplx(reshape(ftb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.) &
            + IU*spread(h%kn, 1, h%nt) &
            *dft2(h, cmplx(reshape(fpb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.)
        fb = fb + reshape(real(dft2(h, filter_mask(h)*v, .false.)), [h%nt*h%np])
    end subroutine grid_filter_t

    function filter_mask(h) result(mk)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp) :: mk(h%nt, h%np)
        integer :: k, l

        do l = 1, h%np
            do k = 1, h%nt
                mk(k, l) = merge(1.0_dp, 0.0_dp, abs(h%km(k)) <= real(h%mnu, dp) .and. &
                                 abs(h%kn(l)) <= real(h%nnu*h%nfp, dp))
            end do
        end do
    end function filter_mask

    !> Transpose of grid_deriv, accumulated into fb.
    subroutine grid_deriv_t(h, ftb, fpb, fb)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: ftb(:), fpb(:)
        real(dp), intent(inout) :: fb(:)
        complex(dp) :: v(h%nt, h%np)

        v = IU*spread(h%km, 2, h%np) &
            *dft2(h, cmplx(reshape(ftb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.) &
            + IU*spread(h%kn, 1, h%nt) &
            *dft2(h, cmplx(reshape(fpb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.)
        fb = fb + reshape(real(dft2(h, v, .false.)), [h%nt*h%np])
    end subroutine grid_deriv_t

    pure real(dp) function relvar(f)
        real(dp), intent(in) :: f(:)

        relvar = (maxval(f) - minval(f))/abs(sum(f)/real(size(f), dp))
    end function relvar

    ! ------------------------------------------------------------------
    ! Jacobian routes: magnetic differential equation in PEST angles
    ! ------------------------------------------------------------------

    !> Phases of the PEST analysis and synthesis at theta_P = theta + lambda.
    subroutine pest_tables(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        complex(dp), allocatable :: e1(:)
        integer :: m, n, npt

        npt = h%nt*h%np
        allocate (sv%emp(0:h%mnu, npt), sv%enp(npt, -h%nnu:h%nnu), e1(npt))
        e1 = exp(IU*(grid_theta(h) + sv%lam))
        sv%emp(0, :) = (1.0_dp, 0.0_dp)
        do m = 1, h%mnu
            sv%emp(m, :) = sv%emp(m - 1, :)*e1
        end do
        e1 = exp(-IU*real(h%nfp, dp)*grid_phi(h))
        sv%enp(:, 0) = (1.0_dp, 0.0_dp)
        do n = 1, h%nnu
            sv%enp(:, n) = sv%enp(:, n - 1)*e1
            sv%enp(:, -n) = conjg(sv%enp(:, n))
        end do
    end subroutine pest_tables

    !> A(m, n) = sum_p v_p exp(i m theta_P) exp(-i n N phi).
    function pest_sum(h, sv, v) result(a)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: v(:)
        complex(dp) :: a(0:h%mnu, -h%nnu:h%nnu)

        a = matmul(sv%emp*spread(v, 1, h%mnu + 1), sv%enp)
    end function pest_sum

    !> Synthesis at the grid points of c (0:mnu, -nnu:nnu): v = Re sum c E,
    !> vth = Re sum i m c E, vph = Re sum (-i n N) c E and optionally the second
    !> derivatives vtt = Re sum (i m)^2 c E, vtp = Re sum (i m)(-i n N) c E,
    !> E = exp(i(m theta_P - n N phi)).
    subroutine pest_synth(h, sv, c, v, vth, vph, vtt, vtp)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        complex(dp), intent(in) :: c(0:, -h%nnu:)
        real(dp), intent(out) :: v(:), vth(:), vph(:)
        real(dp), intent(out), optional :: vtt(:), vtp(:)
        complex(dp), allocatable :: cn(:, :), t1(:, :), t2(:, :)
        real(dp) :: rm(0:h%mnu)
        integer :: m, n, p

        allocate (cn(0:h%mnu, -h%nnu:h%nnu), t1(0:h%mnu, h%nt*h%np))
        allocate (t2(0:h%mnu, h%nt*h%np))
        do n = -h%nnu, h%nnu
            cn(:, n) = -IU*real(n*h%nfp, dp)*c(:, n)
        end do
        do m = 0, h%mnu
            rm(m) = real(m, dp)
        end do
        t1 = matmul(c, transpose(sv%enp))
        t2 = matmul(cn, transpose(sv%enp))
        do p = 1, h%nt*h%np
            v(p) = real(sum(sv%emp(:, p)*t1(:, p)))
            vth(p) = real(sum(IU*rm*sv%emp(:, p)*t1(:, p)))
            vph(p) = real(sum(sv%emp(:, p)*t2(:, p)))
            if (present(vtt)) vtt(p) = -real(sum(rm**2*sv%emp(:, p)*t1(:, p)))
            if (present(vtp)) vtp(p) = real(sum(IU*rm*sv%emp(:, p)*t2(:, p)))
        end do
    end subroutine pest_synth

    module subroutine route_mde(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        real(dp) :: fw(2), jw(2, 6), dist
        integer :: p, m, n, npt, k

        npt = h%nt*h%np
        k = 1
        if (h%route == SFL_JACOBIAN_B2) k = 2
        allocate (sv%u(npt), sv%ju(6, npt), sv%g(npt))
        do p = 1, npt
            call route_weight_kernel([sv%x(1, p), sv%x(2, p), sv%x(3, p), sv%x(5, p), sv%x(6, p), &
                                      sv%f(1, p)], fw, jw)
            sv%u(p) = fw(k)
            sv%ju(:, p) = jw(k, :)
        end do
        sv%ubar = sum(sv%u)/real(npt, dp)
        sv%g = sv%u/sv%ubar - (1.0_dp + sv%x(8, :))
        call pest_tables(h, sv)
        allocate (sv%fmn(0:h%mnu, -h%nnu:h%nnu), sv%cmn(0:h%mnu, -h%nnu:h%nnu))
        allocate (sv%active(0:h%mnu, -h%nnu:h%nnu))
        sv%fmn = (2.0_dp/real(npt, dp))*conjg(pest_sum(h, sv, sv%g))
        sv%cmn = (0.0_dp, 0.0_dp)
        sv%active = .false.
        sv%n_zeroed = 0
        do n = -h%nnu, h%nnu
            do m = 0, h%mnu
                if (m == 0 .and. n <= 0) cycle
                dist = sv%iota*real(m, dp) - real(n*h%nfp, dp)
                if (abs(dist) < h%res_exact*max(1.0_dp, real(m, dp))) then
                    sv%n_zeroed = sv%n_zeroed + 1
                    cycle
                end if
                sv%active(m, n) = .true.
                sv%cmn(m, n) = sv%fmn(m, n)/(IU*dist)
            end do
        end do
        allocate (sv%nuth(npt), sv%nuph(npt), sv%nutt(npt), sv%nutp(npt))
        call pest_synth(h, sv, sv%cmn, sv%nu, sv%nuth, sv%nuph, sv%nutt, sv%nutp)
        sv%nut = sv%nuth*(1.0_dp + sv%x(8, :))
        sv%nup = sv%nuth*sv%x(9, :) + sv%nuph
    end subroutine route_mde

    module subroutine route_mde_jvp(h, sv, xd, lamd, fd, iota_d, nud, nutd, nupd)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: xd(:, :), lamd(:), fd(:, :), iota_d
        real(dp), intent(out) :: nud(:), nutd(:), nupd(:)
        real(dp), allocatable :: ud(:), gd(:), v(:), vth(:), vph(:)
        complex(dp), allocatable :: a(:, :), b(:, :), cd(:, :)
        real(dp) :: ubard, dist
        integer :: p, m, n, npt

        npt = h%nt*h%np
        allocate (ud(npt), gd(npt), v(npt), vth(npt), vph(npt))
        do p = 1, npt
            ud(p) = dot_product(sv%ju(:, p), [xd(1, p), xd(2, p), xd(3, p), xd(5, p), &
                                              xd(6, p), fd(1, p)])
        end do
        ubard = sum(ud)/real(npt, dp)
        gd = ud/sv%ubar - sv%u*ubard/sv%ubar**2 - xd(8, :)
        allocate (a(0:h%mnu, -h%nnu:h%nnu), b(0:h%mnu, -h%nnu:h%nnu))
        allocate (cd(0:h%mnu, -h%nnu:h%nnu))
        a = pest_sum(h, sv, gd)
        b = pest_sum(h, sv, sv%g*lamd)
        cd = (0.0_dp, 0.0_dp)
        do n = -h%nnu, h%nnu
            do m = 0, h%mnu
                if (.not. sv%active(m, n)) cycle
                dist = sv%iota*real(m, dp) - real(n*h%nfp, dp)
                cd(m, n) = ((2.0_dp/real(npt, dp))*conjg(a(m, n) + IU*real(m, dp)*b(m, n)) &
                            - IU*real(m, dp)*iota_d*sv%cmn(m, n))/(IU*dist)
            end do
        end do
        call pest_synth(h, sv, cd, v, vth, vph)
        nud = v + lamd*sv%nuth
        vth = vth + lamd*sv%nutt
        vph = vph + lamd*sv%nutp
        nutd = vth*(1.0_dp + sv%x(8, :)) + sv%nuth*xd(8, :)
        nupd = vth*sv%x(9, :) + sv%nuth*xd(9, :) + vph
    end subroutine route_mde_jvp

    module subroutine route_mde_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: nub(:), nutb(:), nupb(:)
        real(dp), intent(inout) :: xb(:, :), lamb(:), fb(:, :), iotab
        real(dp), allocatable :: nuthb(:), gb(:), ub(:), v(:), vth(:), vph(:)
        complex(dp), allocatable :: a(:, :), b(:, :), c(:, :), fbm(:, :)
        complex(dp) :: cb
        real(dp) :: ubarb, dist
        integer :: m, n, npt

        npt = h%nt*h%np
        allocate (nuthb(npt), gb(npt), ub(npt), v(npt), vth(npt), vph(npt))
        nuthb = nutb*(1.0_dp + sv%x(8, :)) + nupb*sv%x(9, :)
        xb(8, :) = xb(8, :) + nutb*sv%nuth
        xb(9, :) = xb(9, :) + nupb*sv%nuth
        lamb = lamb + nub*sv%nuth + nuthb*sv%nutt + nupb*sv%nutp
        allocate (a(0:h%mnu, -h%nnu:h%nnu), b(0:h%mnu, -h%nnu:h%nnu))
        allocate (c(0:h%mnu, -h%nnu:h%nnu), fbm(0:h%mnu, -h%nnu:h%nnu))
        a = pest_sum(h, sv, nub)
        b = pest_sum(h, sv, nuthb)
        c = pest_sum(h, sv, nupb)
        fbm = (0.0_dp, 0.0_dp)
        do n = -h%nnu, h%nnu
            do m = 0, h%mnu
                if (.not. sv%active(m, n)) cycle
                dist = sv%iota*real(m, dp) - real(n*h%nfp, dp)
                cb = conjg(a(m, n) + IU*real(m, dp)*b(m, n) - IU*real(n*h%nfp, dp)*c(m, n))
                fbm(m, n) = IU*cb/dist
                iotab = iotab + real(conjg(cb)*(-real(m, dp)*sv%cmn(m, n)/dist))
            end do
        end do
        call pest_synth(h, sv, fbm, v, vth, vph)
        gb = (2.0_dp/real(npt, dp))*v
        lamb = lamb + (2.0_dp/real(npt, dp))*sv%g*vth
        ub = gb/sv%ubar
        ubarb = -sum(gb*sv%u)/sv%ubar**2
        ub = ub + ubarb/real(npt, dp)
        xb(8, :) = xb(8, :) - gb
        xb(1, :) = xb(1, :) + ub*sv%ju(1, :)
        xb(2, :) = xb(2, :) + ub*sv%ju(2, :)
        xb(3, :) = xb(3, :) + ub*sv%ju(3, :)
        xb(5, :) = xb(5, :) + ub*sv%ju(4, :)
        xb(6, :) = xb(6, :) + ub*sv%ju(5, :)
        fb(1, :) = fb(1, :) + ub*sv%ju(6, :)
    end subroutine route_mde_vjp

    ! ------------------------------------------------------------------
    ! Hamada II: straight current lines
    ! ------------------------------------------------------------------

    subroutine radial_amp2(h, sv, c, a)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: c(:, 0:)
        real(dp), intent(out) :: a(:)
        integer :: i

        do i = 1, h%mnmax
            a(i) = dot_product(sv%rdd(:, ROW_FULL0 + mod(h%xm(i), 2)), c(i, :))
        end do
    end subroutine radial_amp2

    subroutine radial_amp2_t(h, sv, ab, cb)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: ab(:)
        real(dp), intent(inout) :: cb(:, 0:)
        integer :: i

        do i = 1, h%mnmax
            cb(i, :) = cb(i, :) + ab(i)*sv%rdd(:, ROW_FULL0 + mod(h%xm(i), 2))
        end do
    end subroutine radial_amp2_t

    !> Radial data y(1:9, point) = (R_ss, R_st, R_sp, Z_ss, Z_st, Z_sp,
    !> lambda_s, lambda_st, lambda_sp) and iota', linear in the harmonic arrays.
    subroutine field_data2(h, sv, rc, zs, ls, rs, zc, lc, io, y, iotap)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: rc(:, 0:), zs(:, 0:), ls(:, 0:), io(0:)
        real(dp), intent(in) :: rs(:, 0:), zc(:, 0:), lc(:, 0:)
        real(dp), intent(out) :: y(:, :), iotap
        real(dp), dimension(h%mnmax) :: a1, a1s, a2, a2s, a1ss, a2ss
        real(dp), dimension(h%nt, h%np) :: f, ft, fp
        integer :: npt

        npt = h%nt*h%np
        iotap = dot_product(sv%rd(:, ROW_HALF0), io)
        call radial_amp(h, sv, rc, .false., a1, a1s)
        call radial_amp(h, sv, rs, .false., a2, a2s)
        call radial_amp2(h, sv, rc, a1ss)
        call radial_amp2(h, sv, rs, a2ss)
        call synth(h, a1ss, a2ss, f)
        y(1, :) = reshape(f, [npt])
        call synth(h, a1s, a2s, f, ft, fp)
        y(2, :) = reshape(ft, [npt])
        y(3, :) = reshape(fp, [npt])
        call radial_amp(h, sv, zc, .false., a1, a1s)
        call radial_amp(h, sv, zs, .false., a2, a2s)
        call radial_amp2(h, sv, zc, a1ss)
        call radial_amp2(h, sv, zs, a2ss)
        call synth(h, a1ss, a2ss, f)
        y(4, :) = reshape(f, [npt])
        call synth(h, a1s, a2s, f, ft, fp)
        y(5, :) = reshape(ft, [npt])
        y(6, :) = reshape(fp, [npt])
        call radial_amp(h, sv, lc, .true., a1, a1s)
        call radial_amp(h, sv, ls, .true., a2, a2s)
        call synth(h, a1s, a2s, f, ft, fp)
        y(7, :) = reshape(f, [npt])
        y(8, :) = reshape(ft, [npt])
        y(9, :) = reshape(fp, [npt])
    end subroutine field_data2

    !> Transpose of field_data2 for the stellarator-symmetric arrays.
    subroutine field_data2_t(h, sv, yb, iotapb, rb, zb, lb, ib)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: yb(:, :), iotapb
        real(dp), intent(inout) :: rb(:, 0:), zb(:, 0:), lb(:, 0:), ib(0:)
        real(dp), dimension(h%mnmax) :: acb, asb, zero
        real(dp) :: z2(h%nt, h%np)
        integer :: nt, np

        nt = h%nt
        np = h%np
        z2 = 0.0_dp
        zero = 0.0_dp
        ib = ib + iotapb*sv%rd(:, ROW_HALF0)
        ! R: cos = rmnc
        acb = 0.0_dp
        asb = 0.0_dp
        call synth_t(h, reshape(yb(1, :), [nt, np]), acb=acb, asb=asb)
        call radial_amp2_t(h, sv, acb, rb)
        acb = 0.0_dp
        asb = 0.0_dp
        call synth_t(h, z2, reshape(yb(2, :), [nt, np]), reshape(yb(3, :), [nt, np]), acb, asb)
        call radial_amp_t(h, sv, zero, acb, .false., rb)
        ! Z: sin = zmns
        acb = 0.0_dp
        asb = 0.0_dp
        call synth_t(h, reshape(yb(4, :), [nt, np]), acb=acb, asb=asb)
        call radial_amp2_t(h, sv, asb, zb)
        acb = 0.0_dp
        asb = 0.0_dp
        call synth_t(h, z2, reshape(yb(5, :), [nt, np]), reshape(yb(6, :), [nt, np]), acb, asb)
        call radial_amp_t(h, sv, zero, asb, .false., zb)
        ! lambda: sin = lmns, half mesh
        acb = 0.0_dp
        asb = 0.0_dp
        call synth_t(h, reshape(yb(7, :), [nt, np]), reshape(yb(8, :), [nt, np]), &
                     reshape(yb(9, :), [nt, np]), acb, asb)
        call radial_amp_t(h, sv, zero, asb, .true., lb)
    end subroutine field_data2_t

    !> Radial data and radial field (B_s, d_s B_theta, d_s B_phi) of the primal
    !> surface, I' and G', scaled by psi'.
    subroutine radial_field(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        real(dp) :: psip
        integer :: p, npt

        npt = h%nt*h%np
        psip = -h%phiedge/TWOPI
        allocate (sv%y(9, npt), sv%kr(3, npt), sv%jr(3, NRIN, npt))
        call field_data2(h, sv, h%rmnc, h%zmns, h%lmns, h%rmns, h%zmnc, h%lmnc, h%iotas, &
                         sv%y, sv%iotap)
        do p = 1, npt
            call route_radial_kernel([sv%x(:, p), sv%y(1:6, p), sv%y(8, p), sv%y(9, p), sv%iotap], &
                                     sv%kr(:, p), sv%jr(:, :, p))
        end do
        sv%kr = psip*sv%kr
        sv%jr = psip*sv%jr
        sv%cip = sum(sv%kr(2, :))/real(npt, dp)
        sv%cgp = sum(sv%kr(3, :))/real(npt, dp)
    end subroutine radial_field

    real(dp) function pressure_derivative(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv

        pressure_derivative = dot_product(sv%rd(:, ROW_FULL0), h%presf)
    end function pressure_derivative

    module subroutine route_hamada2(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        real(dp), allocatable :: ws(:), wst(:), wsp(:), etaf(:), etat(:), etap(:)
        real(dp) :: pmax
        integer :: p, npt

        npt = h%nt*h%np
        sv%pprime = pressure_derivative(h, sv)
        pmax = maxval(abs(h%presf))
        if (pmax == 0.0_dp .or. abs(sv%pprime) <= 1.0e-8_dp*pmax) then
            write (error_unit, '(a, f8.5, a, es10.3, a, es10.3, a)') &
                'boozer_hybrid: route sfl_straight_current (Hamada II) at s = ', sv%s, &
                ' needs p'' /= 0, but dp/ds = ', sv%pprime, ' Pa (max |p| = ', pmax, &
                ' Pa). eta = -(beta - <beta>)/(G'' + iota I'') is undefined without '// &
                'pressure gradient; use sfl_jacobian_const (Hamada I) or sfl_covariant.'
            error stop 'boozer_hybrid: sfl_straight_current needs p'' /= 0'
        end if
        call radial_field(h, sv)
        sv%ckp = sv%iotap*sv%ci + sv%iota*sv%cip + sv%cgp
        allocate (ws(npt), wst(npt), wsp(npt), etaf(npt), etat(npt), etap(npt))
        call potential(h, sv%kr(2, :), sv%kr(3, :), ws, wst, wsp)
        allocate (sv%nus(npt), sv%beta(npt), sv%jb(NBIN, npt), sv%db(npt), sv%eta(npt))
        sv%nus = (ws - sv%cip*sv%lam - sv%ci*sv%y(7, :) - sv%ckp*sv%nuc)/sv%ck
        do p = 1, npt
            call route_beta_kernel([sv%kr(1, p), sv%f(2, p), sv%f(3, p), sv%x(8, p), sv%x(9, p), &
                                    sv%y(7, p), sv%nuc(p), sv%nutc(p), sv%nupc(p), sv%nus(p), &
                                    sv%iota, sv%iotap], sv%beta(p), sv%jb(:, p))
        end do
        sv%db = (1.0_dp + sv%x(8, :) + sv%iota*sv%nutc)*(1.0_dp + sv%nupc) &
                - (sv%x(9, :) + sv%iota*sv%nupc)*sv%nutc
        sv%bmean = sum(sv%db*sv%beta)/real(npt, dp)
        sv%den = sv%cgp + sv%iota*sv%cip
        if (abs(sv%den) <= 1.0e-12_dp*abs(sv%ck)) then
            write (error_unit, '(a, f8.5, a, es10.3, a, es10.3)') &
                'boozer_hybrid: route sfl_straight_current at s = ', sv%s, &
                ': G'' + iota I'' = ', sv%den, ' vanishes against iota I + G = ', sv%ck
            error stop 'boozer_hybrid: sfl_straight_current needs G'' + iota I'' /= 0'
        end if
        sv%eta = -(sv%beta - sv%bmean)/sv%den
        call grid_filter(h, sv%eta, etaf, etat, etap)
        sv%nu = sv%nuc + etaf
        sv%nut = sv%nutc + etat
        sv%nup = sv%nupc + etap
    end subroutine route_hamada2

    module subroutine route_hamada2_jvp(h, sv, rd, zd, ld, iod, xd, lamd, fd, iota_d, &
                                        ci_d, ckd, nud, nutd, nupd)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: rd(:, 0:), zd(:, 0:), ld(:, 0:), iod(0:)
        real(dp), intent(in) :: xd(:, :), lamd(:), fd(:, :), iota_d, ci_d, ckd
        real(dp), intent(inout) :: nud(:), nutd(:), nupd(:)
        real(dp), allocatable :: yd(:, :), krd(:, :), z0(:, :), wsd(:), wst(:), wsp(:)
        real(dp), allocatable :: nusd(:), betad(:), dbd(:), etad(:), etatd(:), etapd(:)
        real(dp), allocatable :: etafd(:)
        real(dp) :: iotapd, cipd, cgpd, ckpd, bmeand, dend
        integer :: p, npt

        npt = h%nt*h%np
        allocate (yd(9, npt), krd(3, npt), z0(h%mnmax, 0:h%ns - 1), wsd(npt), wst(npt))
        allocate (wsp(npt), nusd(npt), betad(npt), dbd(npt), etad(npt), etatd(npt))
        allocate (etapd(npt), etafd(npt))
        z0 = 0.0_dp
        call field_data2(h, sv, rd, zd, ld, z0, z0, z0, iod, yd, iotapd)
        do p = 1, npt
            krd(:, p) = matmul(sv%jr(:, :, p), [xd(:, p), yd(1:6, p), yd(8, p), yd(9, p), &
                                                iotapd])
        end do
        cipd = sum(krd(2, :))/real(npt, dp)
        cgpd = sum(krd(3, :))/real(npt, dp)
        ckpd = iotapd*sv%ci + sv%iotap*ci_d + iota_d*sv%cip + sv%iota*cipd + cgpd
        call potential(h, krd(2, :), krd(3, :), wsd, wst, wsp)
        nusd = (wsd - cipd*sv%lam - sv%cip*lamd - ci_d*sv%y(7, :) - sv%ci*yd(7, :) &
                - ckpd*sv%nuc - sv%ckp*nud - sv%nus*ckd)/sv%ck
        do p = 1, npt
            betad(p) = dot_product(sv%jb(:, p), [krd(1, p), fd(2, p), fd(3, p), xd(8, p), &
                                                 xd(9, p), yd(7, p), nud(p), nutd(p), &
                                                 nupd(p), nusd(p), iota_d, iotapd])
        end do
        dbd = (xd(8, :) + iota_d*sv%nutc + sv%iota*nutd)*(1.0_dp + sv%nupc) &
              + (1.0_dp + sv%x(8, :) + sv%iota*sv%nutc)*nupd &
              - (xd(9, :) + iota_d*sv%nupc + sv%iota*nupd)*sv%nutc &
              - (sv%x(9, :) + sv%iota*sv%nupc)*nutd
        bmeand = sum(dbd*sv%beta + sv%db*betad)/real(npt, dp)
        dend = cgpd + iota_d*sv%cip + sv%iota*cipd
        etad = -(betad - bmeand)/sv%den - sv%eta*dend/sv%den
        call grid_filter(h, etad, etafd, etatd, etapd)
        nud = nud + etafd
        nutd = nutd + etatd
        nupd = nupd + etapd
    end subroutine route_hamada2_jvp

    module subroutine route_hamada2_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab, &
                                        cib, ckb, rb, zb, lb, iob)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(inout) :: nub(:), nutb(:), nupb(:)
        real(dp), intent(inout) :: xb(:, :), lamb(:), fb(:, :), iotab, cib, ckb
        real(dp), intent(inout) :: rb(:, 0:), zb(:, 0:), lb(:, 0:), iob(0:)
        real(dp), allocatable :: etab(:), betab(:), dbb(:), ab(:), bqb(:), cb(:), bin(:, :)
        real(dp), allocatable :: nusb(:), krb(:, :), yb(:, :), zero(:), rin(:)
        real(dp) :: bmeanb, denb, cipb, cgpb, ckpb, iotapb
        integer :: p, npt

        npt = h%nt*h%np
        allocate (etab(npt), betab(npt), dbb(npt), ab(npt), bqb(npt), cb(npt))
        allocate (bin(NBIN, npt), nusb(npt), krb(3, npt), yb(9, npt), zero(npt), rin(NRIN))
        zero = 0.0_dp
        ! nu = nu_B + eta and the derivatives
        etab = 0.0_dp
        call grid_filter_t(h, nub, nutb, nupb, etab)
        ! eta = -(beta - <beta>)/den
        betab = -etab/sv%den
        bmeanb = sum(etab)/sv%den
        denb = -sum(etab*sv%eta)/sv%den
        ! <beta> = mean(D_B beta)
        betab = betab + bmeanb*sv%db/real(npt, dp)
        dbb = bmeanb*sv%beta/real(npt, dp)
        ! den = G' + iota I'
        cgpb = denb
        cipb = denb*sv%iota
        iotab = iotab + denb*sv%cip
        ! D_B = A B - C nu_t
        ab = dbb*(1.0_dp + sv%nupc)
        bqb = dbb*(1.0_dp + sv%x(8, :) + sv%iota*sv%nutc)
        cb = -dbb*sv%nutc
        nutb = nutb - dbb*(sv%x(9, :) + sv%iota*sv%nupc) + sv%iota*ab
        xb(8, :) = xb(8, :) + ab
        iotab = iotab + sum(ab*sv%nutc) + sum(cb*sv%nupc)
        nupb = nupb + bqb + sv%iota*cb
        xb(9, :) = xb(9, :) + cb
        ! beta kernel
        do p = 1, npt
            bin(:, p) = betab(p)*sv%jb(:, p)
        end do
        krb = 0.0_dp
        yb = 0.0_dp
        krb(1, :) = bin(1, :)
        fb(2, :) = fb(2, :) + bin(2, :)
        fb(3, :) = fb(3, :) + bin(3, :)
        xb(8, :) = xb(8, :) + bin(4, :)
        xb(9, :) = xb(9, :) + bin(5, :)
        yb(7, :) = bin(6, :)
        nub = nub + bin(7, :)
        nutb = nutb + bin(8, :)
        nupb = nupb + bin(9, :)
        nusb = bin(10, :)
        iotab = iotab + sum(bin(11, :))
        iotapb = sum(bin(12, :))
        ! nu_s = (w_s - I' lambda - I lambda_s - K' nu_B)/K
        cipb = cipb - sum(nusb*sv%lam)/sv%ck
        lamb = lamb - sv%cip*nusb/sv%ck
        cib = cib - sum(nusb*sv%y(7, :))/sv%ck
        yb(7, :) = yb(7, :) - sv%ci*nusb/sv%ck
        ckpb = -sum(nusb*sv%nuc)/sv%ck
        nub = nub - sv%ckp*nusb/sv%ck
        ckb = ckb - sum(nusb*sv%nus)/sv%ck
        ! K' = iota' I + iota I' + G'
        iotapb = iotapb + ckpb*sv%ci
        cib = cib + ckpb*sv%iotap
        iotab = iotab + ckpb*sv%cip
        cipb = cipb + ckpb*sv%iota
        cgpb = cgpb + ckpb
        ! w_s = potential of (d_s B_theta, d_s B_phi); I', G' their means
        call potential_t(h, nusb/sv%ck, zero, zero, krb(2, :), krb(3, :))
        krb(2, :) = krb(2, :) + cipb/real(npt, dp)
        krb(3, :) = krb(3, :) + cgpb/real(npt, dp)
        ! radial field kernel
        do p = 1, npt
            rin = matmul(krb(:, p), sv%jr(:, :, p))
            xb(:, p) = xb(:, p) + rin(1:NHIN)
            yb(1:6, p) = yb(1:6, p) + rin(11:16)
            yb(8, p) = yb(8, p) + rin(17)
            yb(9, p) = yb(9, p) + rin(18)
            iotapb = iotapb + rin(19)
        end do
        call field_data2_t(h, sv, yb, iotapb, rb, zb, lb, iob)
    end subroutine route_hamada2_vjp

    ! ------------------------------------------------------------------
    ! Diagnostics
    ! ------------------------------------------------------------------

    module subroutine hybrid_surface_diag(h, s, d, ds_fd)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: s
        type(hybrid_diag_t), intent(out) :: d
        real(dp), intent(in), optional :: ds_fd
        type(hybrid_surface_t) :: sv, svp, svm
        real(dp), allocatable :: nus(:), bst(:), bsp(:), btt(:), btp(:), bpt(:), bpp(:)
        real(dp), allocatable :: a(:), b(:), c(:), e(:), sqgv(:), jst(:), bts(:), bps(:)
        real(dp), allocatable :: jjp(:)
        real(dp) :: hfd, sp, sm, psip, t1, t2, dist
        integer :: npt, m, n, k, nres

        npt = h%nt*h%np
        hfd = 1.0e-4_dp
        if (present(ds_fd)) hfd = ds_fd
        sp = min(s + hfd, 1.0_dp)
        sm = sp - 2.0_dp*hfd
        call hybrid_surface(h, s, sv, with_spectrum=.false.)
        call hybrid_surface(h, sp, svp, with_spectrum=.false.)
        call hybrid_surface(h, sm, svm, with_spectrum=.false.)
        allocate (nus(npt), a(npt), b(npt), c(npt), e(npt), sqgv(npt), jst(npt))
        allocate (bts(npt), bps(npt), jjp(npt))
        nus = (svp%nu - svm%nu)/(sp - sm)
        if (.not. allocated(sv%kr)) call radial_field(h, sv)
        allocate (bst(npt), bsp(npt), btt(npt), btp(npt), bpt(npt), bpp(npt))
        call grid_deriv(h, sv%kr(1, :), bst, bsp)
        call grid_deriv(h, sv%f(2, :), btt, btp)
        call grid_deriv(h, sv%f(3, :), bpt, bpp)
        associate (x => sv%x, f => sv%f)
            a = 1.0_dp + x(8, :) + sv%iota*sv%nut
            b = x(9, :) + sv%iota*sv%nup
            c = sv%nut
            e = 1.0_dp + sv%nup
            sqgv = x(1, :)*(x(3, :)*x(5, :) - x(2, :)*x(6, :))
            jst = sqgv/sv%dj
            bts = (f(2, :)*e - f(3, :)*c)/sv%dj
            bps = (f(3, :)*a - f(2, :)*b)/sv%dj
            ! sqrt(g)_V J^s, J^theta, J^phi (mu0 = 1) and J* J^phi*
            jjp = ((bpt - btp)*nus + (bsp - sv%kr(3, :))*sv%nut &
                   + (sv%kr(2, :) - bst)*(1.0_dp + sv%nup))/sv%dj
            d%d_jac = relvar(jst)
            d%d_jacb2 = relvar(jst*f(1, :)**2)
            d%d_cov = max(maxval(bts) - minval(bts), maxval(bps) - minval(bps)) &
                      /(abs(sv%cg) + abs(sv%iota*sv%ci))
            d%d_cur = relvar(jjp)
            d%d_cur_n = (maxval(jjp) - minval(jjp))/(abs(sv%cip) + abs(sv%cgp))
        end associate
        d%s = s
        d%iota = sv%iota
        d%ci = sv%ci
        d%cg = sv%cg
        d%cip = sv%cip
        d%cgp = sv%cgp
        d%route = h%route
        d%pprime = pressure_derivative(h, sv)
        if (d%pprime /= 0.0_dp) then
            psip = -h%phiedge/TWOPI
            t1 = MU0*d%pprime*sum(sqgv)/real(npt, dp)
            t2 = psip*(sv%cgp + sv%iota*sv%cip)
            d%eps_fb = abs(t1 + t2)/(abs(t1) + abs(t2))
        end if
        select case (h%route)
        case (SFL_COVARIANT)
            d%defect = d%d_jacb2
        case (SFL_JACOBIAN_B2)
            d%defect = d%d_cov
        case (SFL_JACOBIAN_CONST)
            d%defect = d%d_cur
        case (SFL_STRAIGHT_CURRENT)
            d%defect = d%d_jac
        case default
            d%defect = 0.0_dp
        end select
        if (h%route /= SFL_JACOBIAN_B2 .and. h%route /= SFL_JACOBIAN_CONST) then
            allocate (d%res_m(0), d%res_n(0), d%res_dist(0), d%res_f(0), d%res_amp(0))
            return
        end if
        d%n_zeroed = sv%n_zeroed
        do k = 1, 2
            nres = 0
            do n = -h%nnu, h%nnu
                do m = 0, h%mnu
                    if (m == 0 .and. n <= 0) cycle
                    dist = sv%iota*real(m, dp) - real(n*h%nfp, dp)
                    if (abs(dist) >= h%res_report) cycle
                    nres = nres + 1
                    if (k == 1) cycle
                    d%res_m(nres) = m
                    d%res_n(nres) = n
                    d%res_dist(nres) = dist
                    d%res_f(nres) = abs(sv%fmn(m, n))
                    d%res_amp(nres) = abs(sv%cmn(m, n))
                end do
            end do
            if (k == 1) allocate (d%res_m(nres), d%res_n(nres), d%res_dist(nres), &
                                  d%res_f(nres), d%res_amp(nres))
        end do
    end subroutine hybrid_surface_diag

end submodule boozer_hybrid_routes
