!> Hybrid VMEC -> Boozer construction from VMEC's primary quantities, with
!> exact tangent (JVP) and adjoint (VJP).
!>
!> Only R, Z (full mesh), lambda (half-mesh lmns/lmnc, as stored), iota
!> (half-mesh iotas) and the edge toroidal flux are read. The derived wout
!> arrays bmnc, gmnc, bsub*, bsup*, buco, bvco and the full-mesh iotaf are
!> never used. The construction is
!>  radial   each harmonic c(s) is the quintic spline (spl_reg, libneo's end
!>           conditions) of c_j/s_j^(p/2), p = m mod 2, times s^(p/2): R, Z on
!>           the full mesh (the axis node dropped for odd m), lambda on the
!>           half mesh, iota the same spline of iotas;
!>  angular  exact Fourier synthesis and exact angular derivatives on the
!>           uniform grid theta_a = 2 pi a/n_theta, phi_b = 2 pi b/(N n_phi);
!>  field    |B|, the covariant B_theta, B_phi and |grad s| by the FortSym
!>           kernel boozer_hybrid_field (src/generated), scaled by psi';
!>  Boozer   covariant route: w the Fourier potential of (B_theta, B_phi)
!>           (w_mn = B_theta,mn/(i m), w_0n = B_phi,0n/(-i n N), as
!>           booz_xform), nu = (w - I lambda)/(iota I + G) including the
!>           -I lambda term, and the booz_xform quadrature
!>           f_mn = 2 < D f exp(i(m theta_B - n N phi_B)) > over the VMEC grid
!>           with D the Jacobian of the angle map, so nothing is
!>           re-interpolated.
!>  routes   the covariant route above is SFL_COVARIANT, the default; the
!>           submodule boozer_hybrid_routes adds SFL_JACOBIAN_B2 and
!>           SFL_JACOBIAN_CONST (magnetic differential equation in PEST
!>           angles), SFL_STRAIGHT_CURRENT (Hamada II from the covariant
!>           route and beta) and SFL_PEST, each with tangent and adjoint, and
!>           the route diagnostics (hybrid_surface_diag). They change only
!>           nu; the quadrature is the same. presf (full mesh, an input
!>           profile) is read for the p' /= 0 check of Hamada II.
!> Every stage is linear except the pointwise kernel, the quotient by
!> iota I + G and the phases and D of the quadrature, whose tangents are
!> written out; the linear stages have hand-written transposes.
!>
!> Sign: psi' = -phi_edge/(2 pi) is libneo's torflux (vmecinm_m), so G and I
!> carry the sign of the legacy construction; |B|, nu and |grad s| do not
!> depend on it. The reference is boozer_bench.primary (oracle and hybrid).
!> The grid of the default is n_theta = 4 max(mpol-1, mboz) + 8 and
!> n_phi = 4 max(ntor, nboz) + 8 (n_phi = 1 for axisymmetric files); the
!> other routes (submodule boozer_hybrid_routes) use 4 max(mpol-1, mboz, mnu)
!> + 8 and 4 max(ntor, nboz, nnu) + 8.
module boozer_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_reg

    implicit none
    private

    public :: boozer_hybrid_t, hybrid_surface_t, boozer_hybrid_init
    public :: hybrid_surface, hybrid_surface_jvp, hybrid_surface_vjp
    public :: hybrid_grid, hybrid_covariant_radial
    public :: hybrid_field_kernel, hybrid_radial_row
    public :: hybrid_axis_row, boozer_axis_continuation_weights
    public :: hybrid_radial_second, hybrid_radial_second_t
    public :: hybrid_diag_t, hybrid_surface_diag
    public :: route_weight_kernel, route_radial_kernel, route_beta_kernel
    !> Shared with the routes submodule. GCC 12 otherwise omits their symbols
    !> from the parent object at -O3, leaving the submodule unable to link.
    public :: synth, synth_t, radial_amp, radial_amp_t, dft2, potential, potential_t
    public :: NSPEC, NHIN
    public :: SFL_COVARIANT, SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST, SFL_STRAIGHT_CURRENT, &
              SFL_PEST, sfl_route_name

    integer, parameter :: NHIN = 10, NHOUT = 4, NSPEC = 6
    !> Straight-field-line routes, named by the condition they enforce:
    !> SFL_COVARIANT        Boozer, B_theta and B_phi flux functions ('boozer'),
    !> SFL_JACOBIAN_B2      Boozer, sqrt(g) B^2 flux function,
    !> SFL_JACOBIAN_CONST   Hamada I, sqrt(g) flux function ('hamada'),
    !> SFL_STRAIGHT_CURRENT Hamada II, current lines straight (needs p' /= 0),
    !> SFL_PEST             nu = 0.
    integer, parameter :: SFL_COVARIANT = 1, SFL_JACOBIAN_B2 = 2, SFL_JACOBIAN_CONST = 3, &
                          SFL_STRAIGHT_CURRENT = 4, SFL_PEST = 5
    !> Inputs of the radial field kernel sfl_radial_field beyond the ten of
    !> boozer_hybrid_field, and of sfl_beta.
    integer, parameter :: NRIN = 19, NBIN = 12
    real(dp), parameter :: MU0 = 16.0e-7_dp*atan(1.0_dp)
    integer, parameter :: ROW_FULL0 = 1, ROW_FULL1 = 2, ROW_HALF0 = 3, ROW_HALF1 = 4
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    complex(dp), parameter :: IU = (0.0_dp, 1.0_dp)

    !> Primaries of one wout file and the static tables of the transform.
    type :: boozer_hybrid_t
        integer :: ns = 0, mnmax = 0, nfp = 1, mmax = 0, nmax = 0
        integer :: nt = 0, np = 0, mboz = 0, nboz = 0, mnboz = 0
        logical :: lasym = .false., axisym = .false.
        real(dp) :: phiedge = 0.0_dp, rmajor = 0.0_dp
        !> Route and the nu truncation (0 <= m <= mnu, |n| <= nnu field-period
        !> harmonics) of the magnetic differential equation of the Jacobian
        !> routes. A mode with |iota m - n N| < res_exact max(1, m) is set to
        !> zero (as GVEC), one with |iota m - n N| < res_report is reported.
        integer :: route = SFL_COVARIANT, mnu = 0, nnu = 0
        real(dp) :: res_exact = 1.0e-10_dp, res_report = 1.0e-2_dp
        !> xn and xn_b include nfp; nmax counts field periods.
        integer, allocatable :: xm(:), xn(:), xm_b(:), xn_b(:)
        !> wout shapes (mnmax, 0:ns-1)
        real(dp), allocatable :: rmnc(:, :), zmns(:, :), lmns(:, :)
        real(dp), allocatable :: rmns(:, :), zmnc(:, :), lmnc(:, :), iotas(:)
        !> full-mesh pressure (Pa), an input profile; only Hamada II uses p'
        real(dp), allocatable :: presf(:)
        real(dp), allocatable :: theta(:), phi(:), km(:), kn(:)
        !> exp(i m theta_a) (nt, 0:mmax), i m exp(i m theta_a), exp(-i n N phi_b)
        !> (np, -nmax:nmax), and the DFT matrices exp(-2 pi i k a/n).
        complex(dp), allocatable :: et(:, :), det(:, :), ep(:, :)
        complex(dp), allocatable :: ft(:, :), fp(:, :)
    end type boozer_hybrid_t

    !> Primal state of one surface, kept for the tangent and the adjoint.
    !> nu, nu_t, nu_p are the potential of the route on the VMEC grid,
    !> theta* = theta + lambda + iota nu, phi* = phi + nu; nuc, nutc, nupc the
    !> covariant (Boozer) potential that the routes SFL_COVARIANT and
    !> SFL_STRAIGHT_CURRENT start from.
    type :: hybrid_surface_t
        real(dp) :: s = 0.0_dp, iota = 0.0_dp, ci = 0.0_dp, cg = 0.0_dp, ck = 0.0_dp
        !> rows rv, rd (0:ns-1, 4) and the second derivative rows rdd of the
        !> full mesh (0:ns-1, 2)
        real(dp), allocatable :: rv(:, :), rd(:, :), rdd(:, :)
        !> Scalar-primary continuation below the first positive source knot, indexed
        !> by m >= 2 and full/half mesh. The original four rows stay intact.
        real(dp), allocatable :: axis_rv(:, :, :), axis_rd(:, :, :), axis_rdd(:, :)
        real(dp), allocatable :: x(:, :), lam(:), f(:, :), jac(:, :, :)
        real(dp), allocatable :: nu(:), nut(:), nup(:), thb(:), phb(:), dj(:)
        real(dp), allocatable :: nuc(:), nutc(:), nupc(:)
        complex(dp), allocatable :: em(:, :), en(:, :)
        !> spec(1:6, k): bmnc, bmns, numnc, numns, gradsmnc, gradsmns
        real(dp), allocatable :: spec(:, :)
        !> Jacobian routes: density u = sqrt(g)_V (B^2), its Jacobian ju with
        !> respect to (R, R_s, R_t, Z_s, Z_t, |B|), the mean ubar, the right-hand
        !> side g = u/ubar - (1 + lambda_t) in the VMEC measure, the PEST phases
        !> emp (0:mnu, npt) = exp(i m theta_P), enp (npt, -nnu:nnu) =
        !> exp(-i n N phi), the analysis fmn and solution cmn (0:mnu, -nnu:nnu),
        !> the active-mode mask, d_theta_P and d_phi of nu and their derivatives.
        real(dp) :: ubar = 0.0_dp
        real(dp), allocatable :: u(:), ju(:, :), g(:), nuth(:), nuph(:), nutt(:), nutp(:)
        complex(dp), allocatable :: emp(:, :), enp(:, :), fmn(:, :), cmn(:, :)
        logical, allocatable :: active(:, :)
        integer :: n_zeroed = 0
        !> Hamada II: radial data y (9, npt) = (R_ss, R_st, R_sp, Z_ss, Z_st,
        !> Z_sp, lambda_s, lambda_st, lambda_sp), iota', the radial field kr
        !> (3, npt) = (B_s, d_s B_theta, d_s B_phi) with Jacobian jr, I', G',
        !> K' = (iota I + G)', nu_s of the covariant route, beta with Jacobian
        !> jb, the Boozer angle-map Jacobian db, <beta>_B, G' + iota I' and eta.
        real(dp) :: iotap = 0.0_dp, cip = 0.0_dp, cgp = 0.0_dp, ckp = 0.0_dp
        real(dp) :: bmean = 0.0_dp, den = 0.0_dp, pprime = 0.0_dp
        real(dp), allocatable :: y(:, :), kr(:, :), jr(:, :, :), nus(:), beta(:), jb(:, :)
        real(dp), allocatable :: db(:), eta(:)
    end type hybrid_surface_t

    !> Diagnostics of one surface (primal only, report). Each d_* is the
    !> relative variation (max - min)/|mean| of the property on the VMEC grid:
    !> d_jac of J* = sqrt(g)_V/D, d_jacb2 of J* B^2, d_cur of J* J^phi* (current
    !> lines straight: J^theta*/J^phi* = iota_J, and J* J^phi* = I'/mu0 in
    !> Hamada coordinates); d_cov = max(ptp B_theta*, ptp B_phi*)/(|G| + |iota I|).
    !> d_cur_n = ptp(J* J^phi*)/(|I'| + |G'|) (mu0 = 1) normalizes the current
    !> defect by the flux-function currents instead of their mean, which stays
    !> meaningful when the net toroidal current I' is small (quasi-isodynamic
    !> configurations); Hamada coordinates have J* J^phi* = I', J* J^theta* = -G'.
    !> defect is the property the route does not enforce (SFL_COVARIANT d_jacb2,
    !> SFL_JACOBIAN_B2 d_cov, SFL_JACOBIAN_CONST d_cur, SFL_STRAIGHT_CURRENT
    !> d_jac, SFL_PEST 0). nu_s in J^phi* = J^s nu_s + J^theta nu_t + J^phi
    !> (1 + nu_p) is a central difference of the route potential (step ds_fd).
    !> eps_fb is the averaged radial force balance of the covariant fluxes,
    !> |mu0 p' <sqrt g> + psi'(G' + iota I')| over the sum of the moduli (only
    !> with p' /= 0). The Jacobian routes list the modes with |iota m - n N| <
    !> res_report: res_m, res_n, res_dist = iota m - n N, res_f = |f_mn| (the
    !> resonant residual of the normalized right-hand side) and res_amp =
    !> |f_mn/(iota m - n N)| (zero for zeroed modes).
    type :: hybrid_diag_t
        real(dp) :: s = 0.0_dp, iota = 0.0_dp, ci = 0.0_dp, cg = 0.0_dp
        real(dp) :: cip = 0.0_dp, cgp = 0.0_dp, pprime = 0.0_dp, eps_fb = -1.0_dp
        real(dp) :: d_jac = 0.0_dp, d_jacb2 = 0.0_dp, d_cov = 0.0_dp, d_cur = 0.0_dp
        real(dp) :: d_cur_n = 0.0_dp
        real(dp) :: defect = 0.0_dp
        integer :: route = SFL_COVARIANT, n_zeroed = 0
        integer, allocatable :: res_m(:), res_n(:)
        real(dp), allocatable :: res_dist(:), res_f(:), res_amp(:)
    end type hybrid_diag_t

    interface
        !> Exact radial covector of the covariant route, shared by Hamada II
        !> and the compatibility grid adapter.
        module subroutine hybrid_covariant_radial(h, sv)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(inout) :: sv
        end subroutine hybrid_covariant_radial
        !> FortSym kernels of the routes (gen/gen_sfl_route_kernels.f90):
        !> sqrt(g)_V and sqrt(g)_V |B|^2 with the Jacobian with respect to
        !> (R, R_s, R_t, Z_s, Z_t, |B|);
        module pure subroutine route_weight_kernel(x, f, jac)
            real(dp), intent(in) :: x(6)
            real(dp), intent(out) :: f(2), jac(2, 6)
        end subroutine route_weight_kernel
        !> (B_s, d_s B_theta, d_s B_phi) at psi' = 1 with the Jacobian with
        !> respect to (R, R_s, R_t, R_p, Z_s, Z_t, Z_p, lambda_t, lambda_p, iota,
        !> R_ss, R_st, R_sp, Z_ss, Z_st, Z_sp, lambda_st, lambda_sp, iota');
        module pure subroutine route_radial_kernel(x, f, jac)
            real(dp), intent(in) :: x(NRIN)
            real(dp), intent(out) :: f(3), jac(3, NRIN)
        end subroutine route_radial_kernel
        !> beta with the gradient with respect to (B_s, B_theta, B_phi,
        !> lambda_t, lambda_p, lambda_s, nu, nu_t, nu_p, nu_s, iota, iota').
        module pure subroutine route_beta_kernel(x, f, jac)
            real(dp), intent(in) :: x(NBIN)
            real(dp), intent(out) :: f, jac(NBIN)
        end subroutine route_beta_kernel
        !> Jacobian routes (SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST): nu from the
        !> magnetic differential equation (iota d_thetaP + d_phi) nu = q/<q>_P - 1
        !> in PEST angles (submodule boozer_hybrid_routes).
        module subroutine route_mde(h, sv)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(inout) :: sv
        end subroutine route_mde
        module subroutine route_mde_jvp(h, sv, xd, lamd, fd, iota_d, nud, nutd, nupd)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(in) :: sv
            real(dp), intent(in) :: xd(:, :), lamd(:), fd(:, :), iota_d
            real(dp), intent(out) :: nud(:), nutd(:), nupd(:)
        end subroutine route_mde_jvp
        module subroutine route_mde_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(in) :: sv
            real(dp), intent(in) :: nub(:), nutb(:), nupb(:)
            real(dp), intent(inout) :: xb(:, :), lamb(:), fb(:, :), iotab
        end subroutine route_mde_vjp
        !> Hamada II (SFL_STRAIGHT_CURRENT): nu = nu_B + eta with
        !> eta = -(beta - <beta>_B)/(G' + iota I').
        module subroutine route_hamada2(h, sv)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(inout) :: sv
        end subroutine route_hamada2
        module subroutine route_hamada2_jvp(h, sv, rd, zd, ld, iod, xd, lamd, fd, iota_d, &
                                            ci_d, ckd, nud, nutd, nupd)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(in) :: sv
            real(dp), intent(in) :: rd(:, 0:), zd(:, 0:), ld(:, 0:), iod(0:)
            real(dp), intent(in) :: xd(:, :), lamd(:), fd(:, :), iota_d, ci_d, ckd
            real(dp), intent(inout) :: nud(:), nutd(:), nupd(:)
        end subroutine route_hamada2_jvp
        module subroutine route_hamada2_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab, &
                                            cib, ckb, rb, zb, lb, iob)
            type(boozer_hybrid_t), intent(in) :: h
            type(hybrid_surface_t), intent(in) :: sv
            real(dp), intent(inout) :: nub(:), nutb(:), nupb(:)
            real(dp), intent(inout) :: xb(:, :), lamb(:), fb(:, :), iotab, cib, ckb
            real(dp), intent(inout) :: rb(:, 0:), zb(:, 0:), lb(:, 0:), iob(0:)
        end subroutine route_hamada2_vjp
        !> Diagnostics of surface s for the route of h.
        module subroutine hybrid_surface_diag(h, s, d, ds_fd)
            type(boozer_hybrid_t), intent(in) :: h
            real(dp), intent(in) :: s
            type(hybrid_diag_t), intent(out) :: d
            real(dp), intent(in), optional :: ds_fd
        end subroutine hybrid_surface_diag
    end interface

contains

    include 'generated/boozer_hybrid_field_leaves.inc'
    include 'generated/boozer_axis_continuation_leaves.inc'

    !> |B|, B_theta, B_phi (VMEC covariant) and |grad s| at psi' = 1 and their
    !> Jacobian with respect to x = (R, R_s, R_t, R_p, Z_s, Z_t, Z_p,
    !> lambda_t, lambda_p, iota).
    pure subroutine hybrid_field_kernel(x, f, jac)
        real(dp), intent(in) :: x(NHIN)
        real(dp), intent(out) :: f(NHOUT), jac(NHOUT, NHIN)

        call boozer_hybrid_field(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
                                 x(9), x(10), &
                                 f(1), jac(1, 1), jac(1, 2), jac(1, 3), jac(1, 4), &
                                 jac(1, 5), jac(1, 6), jac(1, 7), jac(1, 8), &
                                 jac(1, 9), jac(1, 10), &
                                 f(2), jac(2, 1), jac(2, 2), jac(2, 3), jac(2, 4), &
                                 jac(2, 5), jac(2, 6), jac(2, 7), jac(2, 8), &
                                 jac(2, 9), jac(2, 10), &
                                 f(3), jac(3, 1), jac(3, 2), jac(3, 3), jac(3, 4), &
                                 jac(3, 5), jac(3, 6), jac(3, 7), jac(3, 8), &
                                 jac(3, 9), jac(3, 10), &
                                 f(4), jac(4, 1), jac(4, 2), jac(4, 3), jac(4, 4), &
                                 jac(4, 5), jac(4, 6), jac(4, 7), jac(4, 8), &
                                 jac(4, 9), jac(4, 10))
    end subroutine hybrid_field_kernel

    !> Read the primaries and build the tables. mboz, nboz as booz_xform
    !> (m = 0..mboz-1, |n| <= nboz field-period harmonics). n_theta and n_phi
    !> override the default grid.
    !> route selects the straight-field-line route (default SFL_COVARIANT);
    !> mnu, nnu (> 0) the nu truncation of the Jacobian routes (default
    !> 4 (mpol - 1) and 4 ntor), res_exact and res_report their resonance
    !> thresholds.
    subroutine boozer_hybrid_init(h, wout_file, mboz, nboz, n_theta, n_phi, route, &
                                  mnu, nnu, res_exact, res_report)
        type(boozer_hybrid_t), intent(out) :: h
        character(len=*), intent(in) :: wout_file
        integer, intent(in) :: mboz, nboz
        integer, intent(in), optional :: n_theta, n_phi, route, mnu, nnu
        real(dp), intent(in), optional :: res_exact, res_report
        integer :: a, m, n, k

        call read_primaries(h, wout_file)
        if (present(route)) h%route = route
        if (h%route < SFL_COVARIANT .or. h%route > SFL_PEST) &
            error stop 'boozer_hybrid_init: unknown route'
        ! 4 (mpol - 1), 4 ntor: the truncation error of nu is at round-off
        ! (1e-14 of B_00 on li383, QA and the sheared ellipstell; 2 (mpol - 1)
        ! leaves 3e-10 to 6e-8)
        h%mnu = 4*h%mmax
        h%nnu = 4*h%nmax
        if (present(mnu)) then
            if (mnu > 0) h%mnu = mnu
        end if
        if (present(nnu)) then
            if (nnu > 0) h%nnu = nnu
        end if
        if (h%axisym) h%nnu = 0
        if (present(res_exact)) h%res_exact = res_exact
        if (present(res_report)) h%res_report = res_report
        h%nt = 4*max(h%mmax, mboz) + 8
        h%np = 4*max(h%nmax, nboz) + 8
        if (h%route /= SFL_COVARIANT) then
            ! The PEST analysis and the truncated eta resolve nu up to mnu, nnu,
            ! and the spectrum in PEST angles is broader than in Boozer angles.
            ! With the default mnu, nnu, doubling this grid changes the spectra
            ! of all routes by 1e-15 of B_00 (test_sfl_routes).
            h%nt = 4*max(h%mmax, mboz, h%mnu) + 8
            h%np = 4*max(h%nmax, nboz, h%nnu) + 8
        end if
        if (present(n_theta)) then
            if (n_theta > 0) h%nt = n_theta
        end if
        if (present(n_phi)) then
            if (n_phi > 0) h%np = n_phi
        end if
        if (h%axisym) h%np = 1
        if (2*(mboz - 1) >= h%nt) error stop 'boozer_hybrid_init: mboz not resolved'
        if (.not. h%axisym .and. 2*nboz >= h%np) &
            error stop 'boozer_hybrid_init: nboz not resolved'
        if (h%route /= SFL_COVARIANT .and. (2*h%mnu >= h%nt .or. &
                                            (.not. h%axisym .and. 2*h%nnu >= h%np))) &
            error stop 'boozer_hybrid_init: mnu, nnu not resolved by the grid'
        h%mboz = mboz
        h%nboz = nboz
        h%mnboz = nboz + 1 + (mboz - 1)*(2*nboz + 1)
        allocate (h%xm_b(h%mnboz), h%xn_b(h%mnboz))
        k = 0
        do m = 0, mboz - 1
            do n = -nboz, nboz
                if (m == 0 .and. n < 0) cycle
                k = k + 1
                h%xm_b(k) = m
                h%xn_b(k) = n*h%nfp
            end do
        end do

        allocate (h%theta(h%nt), h%phi(h%np), h%km(h%nt), h%kn(h%np))
        do a = 1, h%nt
            h%theta(a) = TWOPI*real(a - 1, dp)/real(h%nt, dp)
            h%km(a) = real(fftfreq(a - 1, h%nt), dp)
        end do
        do a = 1, h%np
            h%phi(a) = TWOPI*real(a - 1, dp)/real(h%np*h%nfp, dp)
            h%kn(a) = real(fftfreq(a - 1, h%np)*h%nfp, dp)
        end do
        allocate (h%et(h%nt, 0:h%mmax), h%det(h%nt, 0:h%mmax))
        allocate (h%ep(h%np, -h%nmax:h%nmax))
        do m = 0, h%mmax
            h%et(:, m) = exp(IU*real(m, dp)*h%theta)
            h%det(:, m) = IU*real(m, dp)*h%et(:, m)
        end do
        do n = -h%nmax, h%nmax
            h%ep(:, n) = exp(-IU*real(n*h%nfp, dp)*h%phi)
        end do
        call dft_matrix(h%nt, h%ft)
        call dft_matrix(h%np, h%fp)
    end subroutine boozer_hybrid_init

    pure integer function fftfreq(k, n)
        integer, intent(in) :: k, n

        if (k < (n + 1)/2) then
            fftfreq = k
        else
            fftfreq = k - n
        end if
    end function fftfreq

    subroutine dft_matrix(n, f)
        integer, intent(in) :: n
        complex(dp), allocatable, intent(out) :: f(:, :)
        integer :: j, k

        allocate (f(n, n))
        do k = 1, n
            do j = 1, n
                f(k, j) = exp(-IU*TWOPI*real(mod((k - 1)*(j - 1), n), dp)/real(n, dp))
            end do
        end do
    end subroutine dft_matrix

    subroutine read_primaries(h, wout_file)
        use nctools_module, only: nc_open, nc_close, nc_get, nc_inq_dim
        type(boozer_hybrid_t), intent(inout) :: h
        character(len=*), intent(in) :: wout_file
        integer :: ncid, lens(2), lasym_int
        real(dp), allocatable :: axm(:), axn(:), phi(:), buf(:, :)

        call nc_open(trim(wout_file), ncid)
        call nc_inq_dim(ncid, 'rmnc', lens)
        h%mnmax = lens(1)
        h%ns = lens(2)
        call nc_get(ncid, 'nfp', h%nfp)
        call nc_get(ncid, 'lasym__logical__', lasym_int)
        h%lasym = lasym_int == 1
        call nc_get(ncid, 'Rmajor_p', h%rmajor)
        allocate (axm(h%mnmax), axn(h%mnmax), phi(h%ns), buf(h%mnmax, h%ns))
        call nc_get(ncid, 'xm', axm)
        call nc_get(ncid, 'xn', axn)
        call nc_get(ncid, 'phi', phi)
        h%phiedge = phi(h%ns)
        allocate (h%xm(h%mnmax), h%xn(h%mnmax))
        h%xm = nint(axm)
        h%xn = nint(axn)
        h%mmax = maxval(abs(h%xm))
        h%nmax = maxval(abs(h%xn))/h%nfp
        h%axisym = h%nmax == 0
        allocate (h%iotas(0:h%ns - 1), h%presf(0:h%ns - 1))
        call nc_get(ncid, 'iotas', h%iotas)
        call nc_get(ncid, 'presf', h%presf)
        call get2('rmnc', h%rmnc)
        call get2('zmns', h%zmns)
        call get2('lmns', h%lmns)
        if (h%lasym) then
            call get2('rmns', h%rmns)
            call get2('zmnc', h%zmnc)
            call get2('lmnc', h%lmnc)
        else
            allocate (h%rmns(h%mnmax, 0:h%ns - 1), h%zmnc(h%mnmax, 0:h%ns - 1))
            allocate (h%lmnc(h%mnmax, 0:h%ns - 1))
            h%rmns = 0.0_dp
            h%zmnc = 0.0_dp
            h%lmnc = 0.0_dp
        end if
        call nc_close(ncid)

    contains

        subroutine get2(name, arr)
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(out) :: arr(:, :)

            call nc_get(ncid, name, buf)
            allocate (arr(h%mnmax, 0:h%ns - 1))
            arr = buf
        end subroutine get2
    end subroutine read_primaries

    !> Rows (value, d/ds) of the quintic s-spline of c_j/s_j^(p/2) times
    !> s^(p/2), acting on the wout radial index 0:ns-1. half selects the half
    !> mesh (nodes (j - 1/2) h, j = 1..ns-1), otherwise the full mesh
    !> (nodes j h, with j = 0 dropped for p = 1).
    subroutine hybrid_radial_row(ns, s, half, p, rv, rd, rdd)
        integer, intent(in) :: ns, p
        real(dp), intent(in) :: s
        logical, intent(in) :: half
        real(dp), intent(out) :: rv(0:ns - 1), rd(0:ns - 1)
        real(dp), intent(out), optional :: rdd(0:ns - 1)
        real(dp) :: hs, x0, t, v, dv, ddv, r
        real(dp), allocatable :: coef(:, :)
        integer :: j0, nn, i, j, k

        hs = 1.0_dp/real(ns - 1, dp)
        if (half) then
            j0 = 1
            x0 = 0.5_dp*hs
        else if (p == 1) then
            j0 = 1
            x0 = hs
        else
            j0 = 0
            x0 = 0.0_dp
        end if
        nn = ns - j0
        i = max(0, min(nn - 2, floor((s - x0)/hs)))
        t = s - (x0 + real(i, dp)*hs)
        r = sqrt(s)
        allocate (coef(0:5, nn))
        rv = 0.0_dp
        rd = 0.0_dp
        if (present(rdd)) rdd = 0.0_dp
        do j = 1, nn
            coef = 0.0_dp
            coef(0, j) = 1.0_dp/sqrt(x0 + real(j - 1, dp)*hs)**p
            call spl_reg(5, nn, hs, coef)
            v = 0.0_dp
            dv = 0.0_dp
            do k = 5, 0, -1
                v = coef(k, i + 1) + t*v
            end do
            do k = 5, 1, -1
                dv = real(k, dp)*coef(k, i + 1) + t*dv
            end do
            rv(j0 + j - 1) = v*r**p
            rd(j0 + j - 1) = dv*r**p
            if (p == 1) rd(j0 + j - 1) = rd(j0 + j - 1) + 0.5_dp*v/r
            if (present(rdd)) then
                ddv = 0.0_dp
                do k = 5, 2, -1
                    ddv = real(k*(k - 1), dp)*coef(k, i + 1) + t*ddv
                end do
                rdd(j0 + j - 1) = ddv*r**p
                if (p == 1) rdd(j0 + j - 1) = rdd(j0 + j - 1) + dv/r - 0.25_dp*v/r**3
            end if
        end do
    end subroutine hybrid_radial_row

    !> Scalar Fourier primaries R, Z and lambda must scale as rho**|m|.
    !> Below the first positive source knot, continue in that regular basis with a
    !> quadratic in s matching the parity spline's value and first two
    !> derivatives. Above the join the original operator is unchanged.
    subroutine hybrid_axis_row(ns, s, half, m, rv, rd, rdd)
        integer, intent(in) :: ns, m
        real(dp), intent(in) :: s
        logical, intent(in) :: half
        real(dp), intent(out) :: rv(0:ns - 1), rd(0:ns - 1)
        real(dp), intent(out), optional :: rdd(0:ns - 1)
        real(dp) :: anchor, w(3, 3), va(0:ns - 1), da(0:ns - 1), dda(0:ns - 1)

        if (s <= 0.0_dp) error stop 'hybrid_axis_row: s must be positive'
        anchor = 1.0_dp/real(ns - 1, dp)
        if (half) anchor = 0.5_dp*anchor
        if (abs(m) < 2 .or. s >= anchor) then
            call hybrid_radial_row(ns, s, half, mod(abs(m), 2), rv, rd, rdd)
            return
        end if
        call hybrid_radial_row(ns, anchor, half, mod(abs(m), 2), va, da, dda)
        call boozer_axis_continuation_weights(abs(m), s/anchor, w)
        rv = w(1, 1)*va + w(1, 2)*anchor*da + w(1, 3)*anchor**2*dda
        rd = w(2, 1)*va/anchor + w(2, 2)*da + w(2, 3)*anchor*dda
        if (present(rdd)) &
            rdd = w(3, 1)*va/anchor**2 + w(3, 2)*da/anchor + w(3, 3)*dda
    end subroutine hybrid_axis_row

    !> An even-parity spline reproduces constants, whose second derivative
    !> vanishes. Center its data before contraction to avoid cancellation of
    !> O(ns**2) weights on the major-radius component. Odd parity retains its
    !> s**(1/2) factor and the original row application.
    pure real(dp) function hybrid_radial_second(row, c, parity) result(value)
        real(dp), intent(in) :: row(0:), c(0:)
        integer, intent(in) :: parity
        integer :: anchor

        if (parity == 0) then
            anchor = maxloc(abs(row), dim=1) - 1
            value = dot_product(row, c - c(anchor))
        else
            value = dot_product(row, c)
        end if
    end function hybrid_radial_second

    !> Transpose of the same centered linear operator, including its anchor.
    pure subroutine hybrid_radial_second_t(row, parity, ab, cb)
        real(dp), intent(in) :: row(0:), ab
        integer, intent(in) :: parity
        real(dp), intent(inout) :: cb(0:)
        integer :: anchor

        cb = cb + ab*row
        if (parity == 0) then
            anchor = maxloc(abs(row), dim=1) - 1
            cb(anchor) = cb(anchor) - ab*sum(row)
        end if
    end subroutine hybrid_radial_second_t

    !> Fourier synthesis f = Re sum_i (ac_i - i as_i) exp(i(m theta - n phi)) and
    !> its exact angular derivatives on the grid.
    subroutine synth(h, ac, as, f, ft, fp)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: ac(:), as(:)
        real(dp), intent(out) :: f(:, :)
        real(dp), intent(out), optional :: ft(:, :), fp(:, :)
        complex(dp) :: g0(h%np, 0:h%mmax), gp(h%np, 0:h%mmax), a
        integer :: i, m, n

        g0 = (0.0_dp, 0.0_dp)
        gp = (0.0_dp, 0.0_dp)
        do i = 1, h%mnmax
            m = h%xm(i)
            n = h%xn(i)/h%nfp
            a = cmplx(ac(i), -as(i), kind=dp)
            g0(:, m) = g0(:, m) + a*h%ep(:, n)
            gp(:, m) = gp(:, m) - IU*real(h%xn(i), dp)*a*h%ep(:, n)
        end do
        f = real(matmul(h%et, transpose(g0)))
        if (present(ft)) ft = real(matmul(h%det, transpose(g0)))
        if (present(fp)) fp = real(matmul(h%et, transpose(gp)))
    end subroutine synth

    !> Transpose of synth: amplitude cotangents ab = (Re H, Im H) accumulated,
    !> H_i = sum_ab E_i(a, b) (fb + i m ftb - i xn fpb), E_i = exp(i(m theta - n phi)).
    subroutine synth_t(h, fb, ftb, fpb, acb, asb)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: fb(:, :)
        real(dp), intent(in), optional :: ftb(:, :), fpb(:, :)
        real(dp), intent(inout) :: acb(:), asb(:)
        complex(dp) :: h0(h%np, 0:h%mmax), hp(h%np, 0:h%mmax), z
        integer :: i, m, n

        h0 = matmul(transpose(fb), h%et)
        if (present(ftb)) h0 = h0 + matmul(transpose(ftb), h%det)
        hp = (0.0_dp, 0.0_dp)
        if (present(fpb)) hp = matmul(transpose(fpb), h%et)
        do i = 1, h%mnmax
            m = h%xm(i)
            n = h%xn(i)/h%nfp
            z = sum(h0(:, m)*h%ep(:, n)) - IU*real(h%xn(i), dp)*sum(hp(:, m)*h%ep(:, n))
            acb(i) = acb(i) + real(z)
            asb(i) = asb(i) + aimag(z)
        end do
    end subroutine synth_t

    !> Radial amplitudes of one harmonic array with the rows of the surface.
    subroutine radial_amp(h, sv, c, half, a, as)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: c(:, 0:)
        logical, intent(in) :: half
        real(dp), intent(out) :: a(:), as(:)
        integer :: i, r, mesh, m

        mesh = merge(2, 1, half)
        do i = 1, h%mnmax
            m = abs(h%xm(i))
            if (allocated(sv%axis_rv) .and. m >= 2) then
                a(i) = dot_product(sv%axis_rv(:, m, mesh), c(i, :))
                as(i) = dot_product(sv%axis_rd(:, m, mesh), c(i, :))
                cycle
            end if
            r = row_index(half, mod(h%xm(i), 2))
            a(i) = dot_product(sv%rv(:, r), c(i, :))
            as(i) = dot_product(sv%rd(:, r), c(i, :))
        end do
    end subroutine radial_amp

    subroutine radial_amp_t(h, sv, ab, asb, half, cb)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: ab(:), asb(:)
        logical, intent(in) :: half
        real(dp), intent(inout) :: cb(:, 0:)
        integer :: i, r, mesh, m

        mesh = merge(2, 1, half)
        do i = 1, h%mnmax
            m = abs(h%xm(i))
            if (allocated(sv%axis_rv) .and. m >= 2) then
                cb(i, :) = cb(i, :) + ab(i)*sv%axis_rv(:, m, mesh) &
                            + asb(i)*sv%axis_rd(:, m, mesh)
                cycle
            end if
            r = row_index(half, mod(h%xm(i), 2))
            cb(i, :) = cb(i, :) + ab(i)*sv%rv(:, r) + asb(i)*sv%rd(:, r)
        end do
    end subroutine radial_amp_t

    pure integer function row_index(half, p)
        logical, intent(in) :: half
        integer, intent(in) :: p

        if (half) then
            row_index = ROW_HALF0 + p
        else
            row_index = ROW_FULL0 + p
        end if
    end function row_index

    !> Field data x(1:10, point), lambda on the grid from harmonic arrays
    !> (linear in rc, zs, ls, rs, zc, lc for fixed rows) and iota.
    subroutine field_data(h, sv, rc, zs, ls, rs, zc, lc, io, x, lam, iota)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: rc(:, 0:), zs(:, 0:), ls(:, 0:), io(0:)
        real(dp), intent(in) :: rs(:, 0:), zc(:, 0:), lc(:, 0:)
        real(dp), intent(out) :: x(:, :), lam(:), iota
        real(dp), dimension(h%mnmax) :: a1, a1s, a2, a2s
        real(dp), dimension(h%nt, h%np) :: f, ft, fp

        iota = dot_product(sv%rv(:, ROW_HALF0), io)
        call radial_amp(h, sv, rc, .false., a1, a1s)
        call radial_amp(h, sv, rs, .false., a2, a2s)
        call synth(h, a1, a2, f, ft, fp)
        x(1, :) = reshape(f, [h%nt*h%np])
        x(3, :) = reshape(ft, [h%nt*h%np])
        x(4, :) = reshape(fp, [h%nt*h%np])
        call synth(h, a1s, a2s, f)
        x(2, :) = reshape(f, [h%nt*h%np])
        call radial_amp(h, sv, zc, .false., a1, a1s)
        call radial_amp(h, sv, zs, .false., a2, a2s)
        call synth(h, a1, a2, f, ft, fp)
        x(6, :) = reshape(ft, [h%nt*h%np])
        x(7, :) = reshape(fp, [h%nt*h%np])
        call synth(h, a1s, a2s, f)
        x(5, :) = reshape(f, [h%nt*h%np])
        call radial_amp(h, sv, lc, .true., a1, a1s)
        call radial_amp(h, sv, ls, .true., a2, a2s)
        call synth(h, a1, a2, f, ft, fp)
        lam = reshape(f, [h%nt*h%np])
        x(8, :) = reshape(ft, [h%nt*h%np])
        x(9, :) = reshape(fp, [h%nt*h%np])
        x(10, :) = iota
    end subroutine field_data

    !> Transpose of field_data for the stellarator-symmetric arrays (the
    !> lasym arrays are held fixed): cotangents accumulated into rb, zb, lb, ib.
    subroutine field_data_t(h, sv, xb, lamb, rb, zb, lb, ib)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: xb(:, :), lamb(:)
        real(dp), intent(inout) :: rb(:, 0:), zb(:, 0:), lb(:, 0:), ib(0:)
        real(dp), dimension(h%mnmax) :: acb, asb, acsb, assb
        integer :: nt, np

        nt = h%nt
        np = h%np
        ib = ib + sum(xb(10, :))*sv%rv(:, ROW_HALF0)
        call zero4()
        call synth_t(h, reshape(xb(1, :), [nt, np]), reshape(xb(3, :), [nt, np]), &
                     reshape(xb(4, :), [nt, np]), acb, asb)
        call synth_t(h, reshape(xb(2, :), [nt, np]), acb=acsb, asb=assb)
        ! R: cos = rmnc (acb), sin = rmns (held fixed)
        call radial_amp_t(h, sv, acb, acsb, .false., rb)
        call zero4()
        call synth_t(h, reshape(0.0_dp*xb(1, :), [nt, np]), reshape(xb(6, :), [nt, np]), &
                     reshape(xb(7, :), [nt, np]), acb, asb)
        call synth_t(h, reshape(xb(5, :), [nt, np]), acb=acsb, asb=assb)
        ! Z: sin = zmns (asb)
        call radial_amp_t(h, sv, asb, assb, .false., zb)
        call zero4()
        call synth_t(h, reshape(lamb, [nt, np]), reshape(xb(8, :), [nt, np]), &
                     reshape(xb(9, :), [nt, np]), acb, asb)
        acsb = 0.0_dp
        call radial_amp_t(h, sv, asb, acsb, .true., lb)

    contains

        subroutine zero4()
            acb = 0.0_dp
            asb = 0.0_dp
            acsb = 0.0_dp
            assb = 0.0_dp
        end subroutine zero4
    end subroutine field_data_t

    !> Covariant potential (R1) of (bt, bp): w, w_theta, w_phi on the grid,
    !> w = Re IDFT(D1 DFT(bt) + D2 DFT(bp)), D1 = 1/(i k_m) for k_m /= 0,
    !> D2 = 1/(i k_n) for k_m = 0, k_n /= 0 (numpy fftfreq ordering).
    subroutine potential(h, bt, bp, w, wt, wp)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: bt(:), bp(:)
        real(dp), intent(out) :: w(:), wt(:), wp(:)
        complex(dp) :: y1(h%nt, h%np), y2(h%nt, h%np), wf(h%nt, h%np)
        integer :: k, l

        y1 = dft2(h, cmplx(reshape(bt, [h%nt, h%np]), 0.0_dp, kind=dp), .false.)
        y2 = dft2(h, cmplx(reshape(bp, [h%nt, h%np]), 0.0_dp, kind=dp), .false.)
        do l = 1, h%np
            do k = 1, h%nt
                wf(k, l) = d1(h, k)*y1(k, l) + d2(h, k, l)*y2(k, l)
            end do
        end do
        w = reshape(real(dft2(h, wf, .true.)), [h%nt*h%np])
        wt = reshape(real(dft2(h, IU*spread(h%km, 2, h%np)*wf, .true.)), [h%nt*h%np])
        wp = reshape(real(dft2(h, IU*spread(h%kn, 1, h%nt)*wf, .true.)), [h%nt*h%np])
    end subroutine potential

    !> Transpose of potential: Re DFT(D (IDFT wb + i k_m IDFT wtb + i k_n IDFT wpb)).
    subroutine potential_t(h, wb, wtb, wpb, btb, bpb)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: wb(:), wtb(:), wpb(:)
        real(dp), intent(inout) :: btb(:), bpb(:)
        complex(dp) :: v(h%nt, h%np), v1(h%nt, h%np), v2(h%nt, h%np)
        integer :: k, l

        v = dft2(h, cmplx(reshape(wb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.) &
            + IU*spread(h%km, 2, h%np) &
            *dft2(h, cmplx(reshape(wtb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.) &
            + IU*spread(h%kn, 1, h%nt) &
            *dft2(h, cmplx(reshape(wpb, [h%nt, h%np]), 0.0_dp, kind=dp), .true.)
        do l = 1, h%np
            do k = 1, h%nt
                v1(k, l) = d1(h, k)*v(k, l)
                v2(k, l) = d2(h, k, l)*v(k, l)
            end do
        end do
        btb = btb + reshape(real(dft2(h, v1, .false.)), [h%nt*h%np])
        bpb = bpb + reshape(real(dft2(h, v2, .false.)), [h%nt*h%np])
    end subroutine potential_t

    pure complex(dp) function d1(h, k)
        type(boozer_hybrid_t), intent(in) :: h
        integer, intent(in) :: k

        d1 = (0.0_dp, 0.0_dp)
        if (h%km(k) /= 0.0_dp) d1 = 1.0_dp/(IU*h%km(k))
    end function d1

    pure complex(dp) function d2(h, k, l)
        type(boozer_hybrid_t), intent(in) :: h
        integer, intent(in) :: k, l

        d2 = (0.0_dp, 0.0_dp)
        if (h%km(k) == 0.0_dp .and. h%kn(l) /= 0.0_dp) d2 = 1.0_dp/(IU*h%kn(l))
    end function d2

    !> Two-dimensional DFT (numpy fft2) or its inverse (ifft2).
    function dft2(h, x, inverse) result(y)
        type(boozer_hybrid_t), intent(in) :: h
        complex(dp), intent(in) :: x(:, :)
        logical, intent(in) :: inverse
        complex(dp) :: y(h%nt, h%np)

        if (inverse) then
            y = matmul(matmul(conjg(h%ft), x), conjg(h%fp))/real(h%nt*h%np, dp)
        else
            y = matmul(matmul(h%ft, x), h%fp)
        end if
    end function dft2

    !> Primal transform of surface s; sv keeps the state for jvp/vjp. With
    !> with_spectrum = .false. only the potential of the route and the angle map are
    !> computed (no spectrum).
    subroutine hybrid_surface(h, s, sv, with_spectrum)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: s
        type(hybrid_surface_t), intent(out) :: sv
        logical, intent(in), optional :: with_spectrum
        real(dp) :: psip
        real(dp), allocatable :: w(:), wt(:), wp(:)
        integer :: p, npt, m

        if (s <= 0.0_dp) error stop 'hybrid_surface: s must be positive'
        npt = h%nt*h%np
        sv%s = s
        allocate (sv%rv(0:h%ns - 1, 4), sv%rd(0:h%ns - 1, 4), sv%rdd(0:h%ns - 1, 2))
        call hybrid_radial_row(h%ns, s, .false., 0, sv%rv(:, ROW_FULL0), sv%rd(:, ROW_FULL0), &
                               sv%rdd(:, ROW_FULL0))
        call hybrid_radial_row(h%ns, s, .false., 1, sv%rv(:, ROW_FULL1), sv%rd(:, ROW_FULL1), &
                               sv%rdd(:, ROW_FULL1))
        call hybrid_radial_row(h%ns, s, .true., 0, sv%rv(:, ROW_HALF0), sv%rd(:, ROW_HALF0))
        call hybrid_radial_row(h%ns, s, .true., 1, sv%rv(:, ROW_HALF1), sv%rd(:, ROW_HALF1))
        if (s < 1.0_dp/real(h%ns - 1, dp) .and. h%mmax >= 2) then
            allocate (sv%axis_rv(0:h%ns - 1, 2:h%mmax, 2), &
                      sv%axis_rd(0:h%ns - 1, 2:h%mmax, 2), &
                      sv%axis_rdd(0:h%ns - 1, 2:h%mmax))
            do m = 2, h%mmax
                call hybrid_axis_row(h%ns, s, .false., m, &
                    sv%axis_rv(:, m, 1), sv%axis_rd(:, m, 1), sv%axis_rdd(:, m))
                call hybrid_axis_row(h%ns, s, .true., m, &
                    sv%axis_rv(:, m, 2), sv%axis_rd(:, m, 2))
            end do
        end if
        allocate (sv%x(NHIN, npt), sv%lam(npt), sv%f(NHOUT, npt), sv%jac(NHOUT, NHIN, npt))
        call field_data(h, sv, h%rmnc, h%zmns, h%lmns, h%rmns, h%zmnc, h%lmnc, h%iotas, &
                        sv%x, sv%lam, sv%iota)
        psip = -h%phiedge/TWOPI
        do p = 1, npt
            call hybrid_field_kernel(sv%x(:, p), sv%f(:, p), sv%jac(:, :, p))
        end do
        sv%f(1, :) = abs(psip)*sv%f(1, :)
        sv%jac(1, :, :) = abs(psip)*sv%jac(1, :, :)
        sv%f(2:3, :) = psip*sv%f(2:3, :)
        sv%jac(2:3, :, :) = psip*sv%jac(2:3, :, :)

        sv%ci = sum(sv%f(2, :))/real(npt, dp)
        sv%cg = sum(sv%f(3, :))/real(npt, dp)
        sv%ck = sv%iota*sv%ci + sv%cg
        allocate (sv%nu(npt), sv%nut(npt), sv%nup(npt), sv%thb(npt), sv%phb(npt))
        allocate (sv%dj(npt))
        select case (h%route)
        case (SFL_COVARIANT, SFL_STRAIGHT_CURRENT)
            allocate (w(npt), wt(npt), wp(npt))
            call potential(h, sv%f(2, :), sv%f(3, :), w, wt, wp)
            allocate (sv%nuc(npt), sv%nutc(npt), sv%nupc(npt))
            sv%nuc = (w - sv%ci*sv%lam)/sv%ck
            sv%nutc = (wt - sv%ci*sv%x(8, :))/sv%ck
            sv%nupc = (wp - sv%ci*sv%x(9, :))/sv%ck
            sv%nu = sv%nuc
            sv%nut = sv%nutc
            sv%nup = sv%nupc
            if (h%route == SFL_STRAIGHT_CURRENT) call route_hamada2(h, sv)
        case (SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST)
            call route_mde(h, sv)
        case default
            sv%nu = 0.0_dp
            sv%nut = 0.0_dp
            sv%nup = 0.0_dp
        end select
        sv%thb = reshape(spread(h%theta, 2, h%np), [npt]) + sv%lam + sv%iota*sv%nu
        sv%phb = reshape(spread(h%phi, 1, h%nt), [npt]) + sv%nu
        sv%dj = (1.0_dp + sv%x(8, :) + sv%iota*sv%nut)*(1.0_dp + sv%nup) &
                - (sv%x(9, :) + sv%iota*sv%nup)*sv%nut
        if (present(with_spectrum)) then
            if (.not. with_spectrum) return
        end if
        call phase_tables(h, sv)
        allocate (sv%spec(NSPEC, h%mnboz))
        call project(h, sv, reshape([sv%dj*sv%f(1, :), sv%dj*sv%nu, sv%dj*sv%f(4, :)], &
                                    [npt, 3]), sv%spec)
    end subroutine hybrid_surface

    !> Covariant-route values on an arbitrary uniform consumer angle grid.
    !> The rows are |B|, |grad s|, B_s, theta_B-theta and phi_B-phi, in SI.
    !> Both the moving-node projection and synthesis are spectral. The source
    !> shifts are separately synthesized in VMEC angles, for the forward map.
    subroutine hybrid_grid(h, sv, theta, phi, values, shifts_v, with_radial)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        real(dp), intent(in) :: theta(:), phi(:)
        real(dp), intent(out) :: values(:, :, :), shifts_v(:, :, :)
        logical, intent(in) :: with_radial
        real(dp), allocatable :: u(:, :), spec(:, :)
        complex(dp), allocatable :: et(:, :), ep(:, :), c(:, :)
        integer :: nf, k, j, m, n, npt

        if (h%route /= SFL_COVARIANT) &
            error stop 'hybrid_grid: the compatibility grid supports only sfl_covariant'
        npt = h%nt*h%np
        allocate (u(npt, 5), spec(10, h%mnboz))
        u(:, 1) = sv%f(1, :)
        u(:, 2) = sv%f(4, :)
        u(:, 3) = 0.0_dp
        if (with_radial) then
            if (.not. allocated(sv%beta)) call hybrid_covariant_radial(h, sv)
            u(:, 3) = sv%beta
        end if
        u(:, 4) = sv%lam + sv%iota*sv%nu
        u(:, 5) = sv%nu
        do j = 1, 5
            u(:, j) = sv%dj*u(:, j)
        end do
        if (.not. allocated(sv%em)) call phase_tables(h, sv)
        call project(h, sv, u, spec)
        allocate (et(size(theta), 0:h%mboz - 1))
        allocate (ep(-h%nboz:h%nboz, size(phi)))
        allocate (c(0:h%mboz - 1, -h%nboz:h%nboz))
        do m = 0, h%mboz - 1
            et(:, m) = exp(IU*real(m, dp)*theta)
        end do
        do n = -h%nboz, h%nboz
            ep(n, :) = exp(-IU*real(n*h%nfp, dp)*phi)
        end do
        do nf = 1, 5
            c = 0.0_dp
            do k = 1, h%mnboz
                c(h%xm_b(k), h%xn_b(k)/h%nfp) = &
                    cmplx(spec(2*nf - 1, k), -spec(2*nf, k), dp)
            end do
            values(nf, :, :) = real(matmul(matmul(et, c), ep), dp)
        end do
        call resample_vmec(h, sv%lam + sv%iota*sv%nu, theta, phi, shifts_v(1, :, :))
        call resample_vmec(h, sv%nu, theta, phi, shifts_v(2, :, :))
    end subroutine hybrid_grid

    !> Fourier resampling on the source chart; no moving-node interpolation.
    subroutine resample_vmec(h, v, theta, phi, f)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: v(:), theta(:), phi(:)
        real(dp), intent(out) :: f(:, :)
        complex(dp), allocatable :: et(:, :), ep(:, :), c(:, :)
        integer :: k

        allocate (et(size(theta), h%nt), ep(h%np, size(phi)))
        do k = 1, h%nt
            et(:, k) = exp(IU*h%km(k)*theta)
        end do
        do k = 1, h%np
            ep(k, :) = exp(IU*h%kn(k)*phi)
        end do
        c = matmul(matmul(h%ft, reshape(v, [h%nt, h%np])), transpose(h%fp)) &
            /real(h%nt*h%np, dp)
        f = real(matmul(matmul(et, c), ep), dp)
    end subroutine resample_vmec

    subroutine phase_tables(h, sv)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(inout) :: sv
        complex(dp), allocatable :: e1(:)
        integer :: m, n, npt

        npt = h%nt*h%np
        allocate (sv%em(0:h%mboz - 1, npt), sv%en(npt, -h%nboz:h%nboz), e1(npt))
        e1 = exp(IU*sv%thb)
        sv%em(0, :) = (1.0_dp, 0.0_dp)
        do m = 1, h%mboz - 1
            sv%em(m, :) = sv%em(m - 1, :)*e1
        end do
        e1 = exp(-IU*real(h%nfp, dp)*sv%phb)
        sv%en(:, 0) = (1.0_dp, 0.0_dp)
        do n = 1, h%nboz
            sv%en(:, n) = sv%en(:, n - 1)*e1
            sv%en(:, -n) = conjg(sv%en(:, n))
        end do
    end subroutine phase_tables

    !> Mode weight of the quadrature: 2/npt, 1/npt for (0, 0), 0 for n /= 0 of
    !> an axisymmetric file.
    pure real(dp) function mode_weight(h, k)
        type(boozer_hybrid_t), intent(in) :: h
        integer, intent(in) :: k

        mode_weight = 2.0_dp/real(h%nt*h%np, dp)
        if (h%xm_b(k) == 0 .and. h%xn_b(k) == 0) mode_weight = 0.5_dp*mode_weight
        if (h%axisym .and. h%xn_b(k) /= 0) mode_weight = 0.0_dp
    end function mode_weight

    !> S_k(m, n) = sum_p u_pk exp(i m theta_B) exp(-i n N phi_B) for the three
    !> projected fields k (|B| D, nu D, |grad s| D) in one product, and their
    !> cos/sin parts spec(2k-1:2k, :).
    subroutine project(h, sv, u, spec)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: u(:, :)
        real(dp), intent(out) :: spec(:, :)
        complex(dp), allocatable :: w(:, :), sm(:, :)
        integer :: k, mb, nf

        mb = h%mboz
        nf = size(u, 2)
        allocate (w(nf*mb, h%nt*h%np), sm(nf*mb, -h%nboz:h%nboz))
        do k = 1, nf
            w((k - 1)*mb + 1:k*mb, :) = sv%em*spread(u(:, k), 1, mb)
        end do
        sm = matmul(w, sv%en)
        do k = 1, nf
            call select_modes(h, sm((k - 1)*mb + 1:k*mb, :), spec(2*k - 1, :), spec(2*k, :))
        end do
    end subroutine project

    subroutine select_modes(h, sm, fc, fs)
        type(boozer_hybrid_t), intent(in) :: h
        complex(dp), intent(in) :: sm(0:, -h%nboz:)
        real(dp), intent(out) :: fc(:), fs(:)
        integer :: k
        complex(dp) :: z

        do k = 1, h%mnboz
            z = mode_weight(h, k)*sm(h%xm_b(k), h%xn_b(k)/h%nfp)
            fc(k) = real(z)
            fs(k) = aimag(z)
            if (h%xm_b(k) == 0 .and. h%xn_b(k) == 0) fs(k) = 0.0_dp
        end do
    end subroutine select_modes

    !> Tangent of project: u, u_d the integrands and their tangents, thb_d,
    !> phb_d the angle tangents,
    !> dS = sum_p [u_d + i m u thb_d] E + (-i n N) sum_p u phb_d E.
    subroutine project_jvp(h, sv, u, u_d, thb_d, phb_d, spec_d)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: u(:, :), u_d(:, :), thb_d(:), phb_d(:)
        real(dp), intent(out) :: spec_d(:, :)
        complex(dp), allocatable :: w1(:, :), w2(:, :), a1(:, :), a2(:, :)
        integer :: k, m, n, mb

        mb = h%mboz
        allocate (w1(3*mb, h%nt*h%np), w2(3*mb, h%nt*h%np))
        allocate (a1(3*mb, -h%nboz:h%nboz), a2(3*mb, -h%nboz:h%nboz))
        do k = 1, 3
            do m = 0, mb - 1
                w1((k - 1)*mb + m + 1, :) = sv%em(m, :)*cmplx(u_d(:, k), &
                                            real(m, dp)*u(:, k)*thb_d, kind=dp)
            end do
            w2((k - 1)*mb + 1:k*mb, :) = sv%em*spread(u(:, k)*phb_d, 1, mb)
        end do
        a1 = matmul(w1, sv%en)
        a2 = matmul(w2, sv%en)
        do n = -h%nboz, h%nboz
            a1(:, n) = a1(:, n) - IU*real(n*h%nfp, dp)*a2(:, n)
        end do
        do k = 1, 3
            call select_modes(h, a1((k - 1)*mb + 1:k*mb, :), spec_d(2*k - 1, :), &
                              spec_d(2*k, :))
        end do
    end subroutine project_jvp

    !> Adjoint of project: cotangents of the integrands (ub) and of the angles
    !> (thbb, phbb, accumulated) from the mode cotangents spec_b.
    subroutine project_vjp(h, sv, u, spec_b, ub, thbb, phbb)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: u(:, :), spec_b(:, :)
        real(dp), intent(out) :: ub(:, :)
        real(dp), intent(inout) :: thbb(:), phbb(:)
        complex(dp), allocatable :: sb(:, :), sbn(:, :), q(:, :), qn(:, :)
        complex(dp) :: z
        real(dp) :: rm(0:h%mboz - 1)
        integer :: i, k, m, n, p, mb, r

        mb = h%mboz
        allocate (sb(3*mb, -h%nboz:h%nboz), sbn(3*mb, -h%nboz:h%nboz))
        allocate (q(3*mb, h%nt*h%np), qn(3*mb, h%nt*h%np))
        sb = (0.0_dp, 0.0_dp)
        do k = 1, 3
            do i = 1, h%mnboz
                z = mode_weight(h, i)*cmplx(spec_b(2*k - 1, i), -spec_b(2*k, i), kind=dp)
                if (h%xm_b(i) == 0 .and. h%xn_b(i) == 0) &
                    z = mode_weight(h, i)*spec_b(2*k - 1, i)
                r = (k - 1)*mb + h%xm_b(i) + 1
                sb(r, h%xn_b(i)/h%nfp) = sb(r, h%xn_b(i)/h%nfp) + z
            end do
        end do
        do n = -h%nboz, h%nboz
            sbn(:, n) = -IU*real(n*h%nfp, dp)*sb(:, n)
        end do
        q = matmul(sb, transpose(sv%en))
        qn = matmul(sbn, transpose(sv%en))
        do m = 0, mb - 1
            rm(m) = real(m, dp)
        end do
        do k = 1, 3
            do p = 1, h%nt*h%np
                associate (qk => q((k - 1)*mb + 1:k*mb, p), qnk => qn((k - 1)*mb + 1:k*mb, p))
                    ub(p, k) = real(sum(sv%em(:, p)*qk))
                    thbb(p) = thbb(p) - u(p, k)*aimag(sum(rm*sv%em(:, p)*qk))
                    phbb(p) = phbb(p) + u(p, k)*real(sum(sv%em(:, p)*qnk))
                end associate
            end do
        end do
    end subroutine project_vjp

    !> Tangent of the surface spectrum for tangents of rmnc, zmns, lmns (half
    !> mesh) and iotas at fixed phi_edge (the flux scaling is applied by the
    !> caller). spec_d as sv%spec; iota_d, ci_d, cg_d the profile tangents.
    subroutine hybrid_surface_jvp(h, sv, rd, zd, ld, iod, spec_d, iota_d, ci_d, cg_d)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: rd(:, 0:), zd(:, 0:), ld(:, 0:), iod(0:)
        real(dp), intent(out) :: spec_d(:, :), iota_d, ci_d, cg_d
        real(dp), allocatable :: xd(:, :), lamd(:), fd(:, :), z0(:, :)
        real(dp), allocatable :: w(:), wt(:), wp(:), nud(:), nutd(:), nupd(:)
        real(dp), allocatable :: thbd(:), phbd(:), djd(:)
        real(dp) :: ckd
        integer :: p, npt

        npt = h%nt*h%np
        allocate (xd(NHIN, npt), lamd(npt), fd(NHOUT, npt), z0(h%mnmax, 0:h%ns - 1))
        allocate (w(npt), wt(npt), wp(npt), nud(npt), nutd(npt), nupd(npt))
        allocate (thbd(npt), phbd(npt), djd(npt))
        z0 = 0.0_dp
        call field_data(h, sv, rd, zd, ld, z0, z0, z0, iod, xd, lamd, iota_d)
        do p = 1, npt
            fd(:, p) = matmul(sv%jac(:, :, p), xd(:, p))
        end do
        ci_d = sum(fd(2, :))/real(npt, dp)
        cg_d = sum(fd(3, :))/real(npt, dp)
        ckd = iota_d*sv%ci + sv%iota*ci_d + cg_d
        select case (h%route)
        case (SFL_COVARIANT, SFL_STRAIGHT_CURRENT)
            call potential(h, fd(2, :), fd(3, :), w, wt, wp)
            nud = (w - ci_d*sv%lam - sv%ci*lamd - sv%nuc*ckd)/sv%ck
            nutd = (wt - ci_d*sv%x(8, :) - sv%ci*xd(8, :) - sv%nutc*ckd)/sv%ck
            nupd = (wp - ci_d*sv%x(9, :) - sv%ci*xd(9, :) - sv%nupc*ckd)/sv%ck
            if (h%route == SFL_STRAIGHT_CURRENT) &
                call route_hamada2_jvp(h, sv, rd, zd, ld, iod, xd, lamd, fd, iota_d, ci_d, &
                                       ckd, nud, nutd, nupd)
        case (SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST)
            call route_mde_jvp(h, sv, xd, lamd, fd, iota_d, nud, nutd, nupd)
        case default
            nud = 0.0_dp
            nutd = 0.0_dp
            nupd = 0.0_dp
        end select
        thbd = lamd + iota_d*sv%nu + sv%iota*nud
        phbd = nud
        djd = (xd(8, :) + iota_d*sv%nut + sv%iota*nutd)*(1.0_dp + sv%nup) &
              + (1.0_dp + sv%x(8, :) + sv%iota*sv%nut)*nupd &
              - (xd(9, :) + iota_d*sv%nup + sv%iota*nupd)*sv%nut &
              - (sv%x(9, :) + sv%iota*sv%nup)*nutd
        call project_jvp(h, sv, reshape([sv%dj*sv%f(1, :), sv%dj*sv%nu, sv%dj*sv%f(4, :)], &
                                        [npt, 3]), &
                         reshape([djd*sv%f(1, :) + sv%dj*fd(1, :), djd*sv%nu + sv%dj*nud, &
                                  djd*sv%f(4, :) + sv%dj*fd(4, :)], [npt, 3]), &
                         thbd, phbd, spec_d)
    end subroutine hybrid_surface_jvp

    !> Adjoint of hybrid_surface_jvp: cotangents accumulated into rb, zb, lb, iob.
    subroutine hybrid_surface_vjp(h, sv, spec_b, iota_b, ci_b, cg_b, rb, zb, lb, iob)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: spec_b(:, :), iota_b, ci_b, cg_b
        real(dp), intent(inout) :: rb(:, 0:), zb(:, 0:), lb(:, 0:), iob(0:)
        real(dp), allocatable :: xb(:, :), lamb(:), fb(:, :), ub(:, :), thbb(:), phbb(:)
        real(dp), allocatable :: djb(:), nub(:), nutb(:), nupb(:), ab(:), bqb(:), cb(:)
        real(dp) :: iotab, cib, cgb, ckb
        integer :: p, npt

        npt = h%nt*h%np
        allocate (xb(NHIN, npt), lamb(npt), fb(NHOUT, npt), ub(npt, 3), thbb(npt))
        allocate (phbb(npt), djb(npt), nub(npt), nutb(npt), nupb(npt))
        allocate (ab(npt), bqb(npt), cb(npt))
        xb = 0.0_dp
        fb = 0.0_dp
        thbb = 0.0_dp
        phbb = 0.0_dp
        call project_vjp(h, sv, reshape([sv%dj*sv%f(1, :), sv%dj*sv%nu, sv%dj*sv%f(4, :)], &
                                        [npt, 3]), spec_b, ub, thbb, phbb)
        djb = ub(:, 1)*sv%f(1, :) + ub(:, 2)*sv%nu + ub(:, 3)*sv%f(4, :)
        fb(1, :) = ub(:, 1)*sv%dj
        nub = ub(:, 2)*sv%dj
        fb(4, :) = ub(:, 3)*sv%dj
        ! theta_B = theta + lambda + iota nu, phi_B = phi + nu
        lamb = thbb
        iotab = iota_b + sum(thbb*sv%nu)
        nub = nub + sv%iota*thbb + phbb
        ! D = A B - C nu_t
        ab = djb*(1.0_dp + sv%nup)
        bqb = djb*(1.0_dp + sv%x(8, :) + sv%iota*sv%nut)
        cb = -djb*sv%nut
        nutb = -djb*(sv%x(9, :) + sv%iota*sv%nup) + sv%iota*ab
        xb(8, :) = xb(8, :) + ab
        iotab = iotab + sum(ab*sv%nut) + sum(cb*sv%nup)
        nupb = bqb + sv%iota*cb
        xb(9, :) = xb(9, :) + cb
        cib = ci_b
        cgb = cg_b
        ckb = 0.0_dp
        select case (h%route)
        case (SFL_COVARIANT, SFL_STRAIGHT_CURRENT)
            if (h%route == SFL_STRAIGHT_CURRENT) &
                call route_hamada2_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab, cib, ckb, &
                                       rb, zb, lb, iob)
            ! nu = (w - I lambda)/K and its derivatives
            cib = cib - (sum(nub*sv%lam) + sum(nutb*sv%x(8, :)) &
                         + sum(nupb*sv%x(9, :)))/sv%ck
            lamb = lamb - sv%ci*nub/sv%ck
            xb(8, :) = xb(8, :) - sv%ci*nutb/sv%ck
            xb(9, :) = xb(9, :) - sv%ci*nupb/sv%ck
            ckb = ckb - (sum(nub*sv%nuc) + sum(nutb*sv%nutc) + sum(nupb*sv%nupc))/sv%ck
            call potential_t(h, nub/sv%ck, nutb/sv%ck, nupb/sv%ck, fb(2, :), fb(3, :))
        case (SFL_JACOBIAN_B2, SFL_JACOBIAN_CONST)
            call route_mde_vjp(h, sv, nub, nutb, nupb, xb, lamb, fb, iotab)
        end select
        iotab = iotab + ckb*sv%ci
        cib = cib + ckb*sv%iota
        cgb = cgb + ckb
        fb(2, :) = fb(2, :) + cib/real(npt, dp)
        fb(3, :) = fb(3, :) + cgb/real(npt, dp)
        do p = 1, npt
            xb(:, p) = xb(:, p) + matmul(fb(:, p), sv%jac(:, :, p))
        end do
        xb(10, :) = xb(10, :) + iotab/real(npt, dp)
        call field_data_t(h, sv, xb, lamb, rb, zb, lb, iob)
    end subroutine hybrid_surface_vjp

    !> Name of a route as in the design (sfl_covariant, ...), blank padded.
    pure function sfl_route_name(route) result(name)
        integer, intent(in) :: route
        character(len=20) :: name

        select case (route)
        case (SFL_COVARIANT)
            name = 'sfl_covariant'
        case (SFL_JACOBIAN_B2)
            name = 'sfl_jacobian_b2'
        case (SFL_JACOBIAN_CONST)
            name = 'sfl_jacobian_const'
        case (SFL_STRAIGHT_CURRENT)
            name = 'sfl_straight_current'
        case (SFL_PEST)
            name = 'sfl_pest'
        case default
            name = 'unknown'
        end select
    end function sfl_route_name

end module boozer_hybrid
