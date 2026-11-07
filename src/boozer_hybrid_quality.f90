!> Equilibrium-quality checks of a VMEC wout file for the hybrid Boozer
!> construction, from the primaries only (R, Z, lambda, iota, phi, presf).
!>
!> Metrics (rho window [rho_min, rho_max], n_check surfaces uniform in rho):
!>  eps_rad  radial resolution: every second full-mesh node of R, Z held out,
!>           the quintic spline (as the construction uses it) through the
!>           others evaluated there, max error over modes and held-out nodes in
!>           the window over the minor radius |rmnc(1, 0)| at the edge;
!>  eps_fb   averaged radial force balance (design/sfl.md H3),
!>           |a + b + c|/(|a| + |b| + |c|), a = mu0 p' <sqrt(g)_V>,
!>           b = psi' <d_s B_phi (1 + lambda_t)>, c = psi' <d_s B_theta (iota -
!>           lambda_p)>, d_s at fixed VMEC angles; median over the window;
!>  eps_js   surface current, max |d_theta B_phi - d_phi B_theta|/|G + iota I|
!>           on the VMEC grid (sqrt(g) J^s); median over the window;
!>  eps_lam  lambda-B consistency: lambda recomputed on the surface from the
!>           weak form of J^s = 0 on VMEC's lambda modes (as GVEC
!>           recompute_lambda, lambda_solve.F90), max |lambda_rec - lambda|
!>           (rad) on the grid; max over the window;
!>  d_fam    covariant (sfl_covariant) against Jacobian (sfl_jacobian_b2)
!>           Boozer route, max |b_mn(1) - b_mn(2)|/|b_00|, the within-family
!>           difference that vanishes for J^s = 0; median over the window;
!>  res_f    resonant residual |f_mn| of the Jacobian route at the rational
!>           surfaces iota = n N/m (m <= mpol - 1, |n| <= ntor) in the window,
!>           max over the harmonics (k m, k n); max over the window.
!> A metric above its tolerance is an error (sfl_raise: message and error
!> stop) unless allow_poor_equilibrium is set, which keeps the status with
!> overridden = .true., prints the message once with the prefix
!> 'sfl: OVERRIDE (<reason>)' and lets the caller record it. A negative
!> tolerance makes the metric report-only. Calibration: boozer-bench
!> reports/libneo-quality.md.
module boozer_hybrid_quality
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, boozer_hybrid_init, &
                             hybrid_surface, hybrid_radial_row, SFL_COVARIANT, &
                             SFL_JACOBIAN_B2
    use spl_three_to_five_sub, only: spl_reg

    implicit none
    private

    public :: sfl_quality_config_t, sfl_status_t, sfl_quality_t
    public :: sfl_quality_check, sfl_raise, sfl_quality_message
    public :: SFL_OK, SFL_ERR_RESOLUTION, SFL_ERR_FORCE_BALANCE, SFL_ERR_SURFACE_CURRENT, &
              SFL_ERR_LAMBDA, SFL_ERR_FAMILY, SFL_ERR_RESONANCE, SFL_ERR_INPUT, NMETRIC, &
              sfl_metric_name

    integer, parameter :: SFL_OK = 0, SFL_ERR_RESOLUTION = 1, SFL_ERR_FORCE_BALANCE = 2, &
                          SFL_ERR_SURFACE_CURRENT = 3, SFL_ERR_LAMBDA = 4, &
                          SFL_ERR_FAMILY = 5, SFL_ERR_RESONANCE = 6, SFL_ERR_INPUT = 7
    integer, parameter :: NMETRIC = 6
    real(dp), parameter :: MU0 = 16.0e-7_dp*atan(1.0_dp)
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    complex(dp), parameter :: IU = (0.0_dp, 1.0_dp)

    !> tol(k) of metric k = SFL_ERR_RESOLUTION .. SFL_ERR_RESONANCE (eps_rad,
    !> eps_fb, eps_js, eps_lam, d_fam, res_f); negative: report only.
    type :: sfl_quality_config_t
        real(dp) :: rho_window(2) = [0.2_dp, 0.95_dp]
        integer :: n_check = 16
        !> Calibrated on VMEC++ solves of six public inputs (boozer-bench
        !> reports/libneo-quality.md): about twice the largest value at ns 136,
        !> 201, 406 (ftol 1e-15, also mpol + 4, ntor + 2); every file at ns <=
        !> 46 whose Boozer |B| spectrum differs from ns 406 by more than 3.5e-4
        !> of B_00 fails. eps_js is a gross check (mpol-limited on li383).
        real(dp) :: tol(NMETRIC) = [1.0e-5_dp, 1.0e-3_dp, 5.0e-2_dp, 6.0e-4_dp, 1.0e-5_dp, &
                                    2.0e-4_dp]
        real(dp) :: ds_fd = 1.0e-4_dp
        logical :: allow_poor_equilibrium = .false.
        character(len=256) :: override_reason = ''
    end type sfl_quality_config_t

    type :: sfl_status_t
        integer :: code = SFL_OK
        logical :: overridden = .false.
        integer :: m_res = 0, n_res = 0
        real(dp) :: s = -1.0_dp, value = 0.0_dp, tol = 0.0_dp, res_dist = 0.0_dp
        character(len=:), allocatable :: message
    end type sfl_status_t

    !> Per-surface metrics on the check surfaces, the file-level eps_rad, the
    !> rational surfaces with their resonant residual, the aggregated value of
    !> each metric (value(k)) with its worst surface (s_worst(k)), the wout
    !> convergence record and the status.
    type :: sfl_quality_t
        integer :: ns = 0, mpol = 0, ntor = 0, nfp = 0
        character(len=256) :: wout = ''
        real(dp) :: fsqr = -1.0_dp, fsqz = -1.0_dp, fsql = -1.0_dp, ftolv = -1.0_dp
        real(dp) :: rho_window(2) = 0.0_dp, tol(NMETRIC) = 0.0_dp
        real(dp) :: eps_rad = 0.0_dp, s_rad = -1.0_dp
        real(dp), allocatable :: s(:), rho(:), iota(:), eps_fb(:), eps_js(:), eps_lam(:)
        real(dp), allocatable :: d_fam(:)
        integer, allocatable :: res_m(:), res_n(:)
        real(dp), allocatable :: res_s(:), res_f(:)
        real(dp) :: value(NMETRIC) = 0.0_dp, s_worst(NMETRIC) = -1.0_dp
        logical :: overridden = .false.
        character(len=256) :: override_reason = ''
        type(sfl_status_t) :: status
    end type sfl_quality_t

    !> The last checked file and configuration with its record: a repeated
    !> check of the same file name and configuration in one process returns the
    !> record without recomputing or printing (the override is printed once per
    !> file). The cache assumes the file does not change during the run.
    character(len=1024), save :: cache_file = ''
    type(sfl_quality_config_t), save :: cache_cfg
    type(sfl_quality_t), save :: cache_q

contains

    pure function sfl_metric_name(k) result(name)
        integer, intent(in) :: k
        character(len=8) :: name

        select case (k)
        case (SFL_ERR_RESOLUTION)
            name = 'eps_rad'
        case (SFL_ERR_FORCE_BALANCE)
            name = 'eps_fb'
        case (SFL_ERR_SURFACE_CURRENT)
            name = 'eps_js'
        case (SFL_ERR_LAMBDA)
            name = 'eps_lam'
        case (SFL_ERR_FAMILY)
            name = 'd_fam'
        case (SFL_ERR_RESONANCE)
            name = 'res_f'
        case default
            name = 'input'
        end select
    end function sfl_metric_name

    !> Run all checks on wout_file. Without status, a failing check stops the
    !> program (sfl_raise) unless cfg%allow_poor_equilibrium.
    subroutine sfl_quality_check(wout_file, q, cfg, status)
        character(len=*), intent(in) :: wout_file
        type(sfl_quality_t), intent(out) :: q
        type(sfl_quality_config_t), intent(in), optional :: cfg
        type(sfl_status_t), intent(out), optional :: status
        type(sfl_quality_config_t) :: c
        type(boozer_hybrid_t) :: h1, h2
        integer :: mb, nb

        if (present(cfg)) c = cfg
        if (len_trim(cache_file) > 0 .and. trim(wout_file) == trim(cache_file) .and. &
            same_config(c, cache_cfg)) then
            q = cache_q
            call finish()
            return
        end if
        q%wout = wout_file
        q%rho_window = c%rho_window
        q%tol = c%tol
        q%overridden = .false.
        q%override_reason = c%override_reason
        call read_convergence(wout_file, q)
        call boozer_hybrid_init(h1, wout_file, 1, 0, route=SFL_COVARIANT)
        mb = h1%mmax + 1
        nb = h1%nmax
        call boozer_hybrid_init(h1, wout_file, mb, nb, route=SFL_COVARIANT)
        call boozer_hybrid_init(h2, wout_file, mb, nb, route=SFL_JACOBIAN_B2)
        q%ns = h1%ns
        q%nfp = h1%nfp
        call radial_holdout(h1, c%rho_window, q%eps_rad, q%s_rad)
        call surface_metrics(h1, h2, c, q)
        call rational_surfaces(h2, c, q)
        call aggregate(q)
        call evaluate(q, c)
        cache_file = wout_file
        cache_cfg = c
        cache_q = q
        call finish()

    contains

        subroutine finish()
            if (present(status)) then
                status = q%status
            else if (q%status%code /= SFL_OK) then
                call sfl_raise(q%status)
            end if
        end subroutine finish
    end subroutine sfl_quality_check

    pure logical function same_config(a, b)
        type(sfl_quality_config_t), intent(in) :: a, b

        same_config = all(a%rho_window == b%rho_window) .and. a%n_check == b%n_check &
                      .and. all(a%tol == b%tol) .and. a%ds_fd == b%ds_fd .and. &
                      (a%allow_poor_equilibrium .eqv. b%allow_poor_equilibrium) .and. &
                      a%override_reason == b%override_reason
    end function same_config

    !> Print the message to error_unit and stop, unless overridden (then the
    !> message is printed once with the OVERRIDE prefix by the check).
    subroutine sfl_raise(status)
        type(sfl_status_t), intent(in) :: status

        if (status%code == SFL_OK .or. status%overridden) return
        write (error_unit, '(a)') status%message
        error stop 'sfl: inconsistent equilibrium (see message above)'
    end subroutine sfl_raise

    function sfl_quality_message(q) result(msg)
        type(sfl_quality_t), intent(in) :: q
        character(len=:), allocatable :: msg

        if (allocated(q%status%message)) then
            call append(msg, q%status%message)
        else
            call append(msg, '')
        end if
    end function sfl_quality_message

    !> s = s // v for a deferred-length string (libneo builds with
    !> -fno-realloc-lhs); an unallocated s starts empty.
    subroutine append(s, v)
        character(len=:), allocatable, intent(inout) :: s
        character(len=*), intent(in) :: v
        character(len=:), allocatable :: tmp
        integer :: n

        n = 0
        if (allocated(s)) n = len(s)
        allocate (character(len=n + len(v)) :: tmp)
        if (n > 0) tmp(1:n) = s
        tmp(n + 1:) = v
        call move_alloc(tmp, s)
    end subroutine append

    subroutine read_convergence(wout_file, q)
        use netcdf
        character(len=*), intent(in) :: wout_file
        type(sfl_quality_t), intent(inout) :: q
        integer :: ncid

        if (nf90_open(trim(wout_file), NF90_NOWRITE, ncid) /= NF90_NOERR) return
        call get_r('fsqr', q%fsqr)
        call get_r('fsqz', q%fsqz)
        call get_r('fsql', q%fsql)
        call get_r('ftolv', q%ftolv)
        call get_i('mpol', q%mpol)
        call get_i('ntor', q%ntor)
        if (nf90_close(ncid) /= NF90_NOERR) return

    contains

        subroutine get_r(name, v)
            character(len=*), intent(in) :: name
            real(dp), intent(inout) :: v
            integer :: id

            if (nf90_inq_varid(ncid, name, id) /= NF90_NOERR) return
            if (nf90_get_var(ncid, id, v) /= NF90_NOERR) v = -1.0_dp
        end subroutine get_r

        subroutine get_i(name, v)
            character(len=*), intent(in) :: name
            integer, intent(inout) :: v
            integer :: id

            if (nf90_inq_varid(ncid, name, id) /= NF90_NOERR) return
            if (nf90_get_var(ncid, id, v) /= NF90_NOERR) v = 0
        end subroutine get_i
    end subroutine read_convergence

    !> eps_rad: the spline of every second full-mesh node (j = 0, 2, 4, ...,
    !> uniform with spacing 2 h, the axis dropped for odd m) at the held-out
    !> odd nodes inside the window and the kept mesh.
    subroutine radial_holdout(h, win, eps, s_worst)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: win(2)
        real(dp), intent(out) :: eps, s_worst
        real(dp), allocatable :: rv(:, :), rd(:), sub(:, :)
        real(dp) :: smax, sj, a, e
        integer :: nsub, j, k, p, kk, ia

        nsub = (h%ns + 1)/2
        smax = real(2*(nsub - 1), dp)/real(h%ns - 1, dp)
        allocate (rv(0:nsub - 1, 0:1), rd(0:nsub - 1), sub(h%mnmax, 0:nsub - 1))
        ia = 0
        do k = 1, h%mnmax
            if (h%xm(k) == 1 .and. h%xn(k) == 0) ia = k
        end do
        a = 1.0_dp
        if (ia > 0) a = abs(h%rmnc(ia, h%ns - 1))
        eps = 0.0_dp
        s_worst = -1.0_dp
        do j = 1, h%ns - 1, 2
            sj = real(j, dp)/real(h%ns - 1, dp)
            if (sj > smax .or. sqrt(sj) < win(1) .or. sqrt(sj) > win(2)) cycle
            do p = 0, 1
                call hybrid_radial_row(nsub, sj/smax, .false., p, rv(:, p), rd)
            end do
            do kk = 1, 4
                select case (kk)
                case (1)
                    sub = h%rmnc(:, 0:h%ns - 1:2)
                case (2)
                    sub = h%zmns(:, 0:h%ns - 1:2)
                case (3)
                    if (.not. h%lasym) cycle
                    sub = h%rmns(:, 0:h%ns - 1:2)
                case (4)
                    if (.not. h%lasym) cycle
                    sub = h%zmnc(:, 0:h%ns - 1:2)
                end select
                do k = 1, h%mnmax
                    p = mod(h%xm(k), 2)
                    e = abs(dot_product(rv(:, p), sub(k, :)) - field(kk, k, j))
                    if (e > eps) then
                        eps = e
                        s_worst = sj
                    end if
                end do
            end do
        end do
        eps = eps/a

    contains

        real(dp) function field(kk, k, j)
            integer, intent(in) :: kk, k, j

            select case (kk)
            case (1)
                field = h%rmnc(k, j)
            case (2)
                field = h%zmns(k, j)
            case (3)
                field = h%rmns(k, j)
            case default
                field = h%zmnc(k, j)
            end select
        end function field
    end subroutine radial_holdout

    subroutine surface_metrics(h1, h2, c, q)
        type(boozer_hybrid_t), intent(in) :: h1, h2
        type(sfl_quality_config_t), intent(in) :: c
        type(sfl_quality_t), intent(inout) :: q
        type(hybrid_surface_t) :: sv, svp, svm, sv2
        real(dp) :: rho, b00
        integer :: k, nc

        nc = max(c%n_check, 1)
        allocate (q%s(nc), q%rho(nc), q%iota(nc), q%eps_fb(nc), q%eps_js(nc))
        allocate (q%eps_lam(nc), q%d_fam(nc))
        do k = 1, nc
            rho = c%rho_window(1)
            if (nc > 1) rho = rho + (c%rho_window(2) - c%rho_window(1)) &
                              *real(k - 1, dp)/real(nc - 1, dp)
            q%rho(k) = rho
            q%s(k) = rho**2
            call hybrid_surface(h1, q%s(k), sv, with_spectrum=.true.)
            call hybrid_surface(h1, min(q%s(k) + c%ds_fd, 1.0_dp), svp, with_spectrum=.false.)
            call hybrid_surface(h1, min(q%s(k) + c%ds_fd, 1.0_dp) - 2.0_dp*c%ds_fd, svm, &
                                with_spectrum=.false.)
            call hybrid_surface(h2, q%s(k), sv2, with_spectrum=.true.)
            q%iota(k) = sv%iota
            q%eps_fb(k) = force_balance(h1, sv, svp, svm)
            q%eps_js(k) = surface_current(h1, sv)
            q%eps_lam(k) = lambda_defect(h1, sv)
            b00 = abs(sv%spec(1, 1))
            q%d_fam(k) = maxval(abs(sv%spec(1:2, :) - sv2%spec(1:2, :)))/b00
        end do
    end subroutine surface_metrics

    real(dp) function force_balance(h, sv, svp, svm) result(eps)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv, svp, svm
        real(dp), allocatable :: dbt(:), dbp(:), sqgv(:)
        real(dp) :: rv(0:h%ns - 1), rd(0:h%ns - 1), pp, psip, a, b, c
        integer :: npt

        npt = h%nt*h%np
        allocate (dbt(npt), dbp(npt), sqgv(npt))
        dbt = (svp%f(2, :) - svm%f(2, :))/(svp%s - svm%s)
        dbp = (svp%f(3, :) - svm%f(3, :))/(svp%s - svm%s)
        associate (x => sv%x)
            sqgv = x(1, :)*(x(3, :)*x(5, :) - x(2, :)*x(6, :))
            call hybrid_radial_row(h%ns, sv%s, .false., 0, rv, rd)
            pp = dot_product(rd, h%presf)
            psip = -h%phiedge/TWOPI
            a = MU0*pp*sum(sqgv)/real(npt, dp)
            b = psip*sum(dbp*(1.0_dp + x(8, :)))/real(npt, dp)
            c = psip*sum(dbt*(sv%iota - x(9, :)))/real(npt, dp)
        end associate
        eps = abs(a + b + c)/(abs(a) + abs(b) + abs(c))
    end function force_balance

    real(dp) function surface_current(h, sv) result(eps)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), allocatable :: d1(:), d2(:)

        allocate (d1(h%nt*h%np), d2(h%nt*h%np))
        d1 = grid_derivative(h, sv%f(3, :), 1)
        d2 = grid_derivative(h, sv%f(2, :), 2)
        eps = maxval(abs(d1 - d2))/abs(sv%ck)
    end function surface_current

    !> Exact angular derivative (dir 1: theta, 2: phi) of the trigonometric
    !> interpolant of f on the grid, the Nyquist harmonic dropped.
    function grid_derivative(h, f, dir) result(df)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: f(:)
        integer, intent(in) :: dir
        real(dp), allocatable :: df(:)
        complex(dp), allocatable :: y(:, :)
        real(dp), allocatable :: k(:)
        integer :: l

        allocate (y(h%nt, h%np), df(h%nt*h%np))
        y = matmul(h%ft, matmul(cmplx(reshape(f, [h%nt, h%np]), 0.0_dp, kind=dp), &
                                transpose(h%fp)))
        if (dir == 1) then
            allocate (k(h%nt))
            k = h%km
            if (mod(h%nt, 2) == 0) k(h%nt/2 + 1) = 0.0_dp
            do l = 1, h%np
                y(:, l) = IU*k*y(:, l)
            end do
        else
            allocate (k(h%np))
            k = h%kn
            if (mod(h%np, 2) == 0) k(h%np/2 + 1) = 0.0_dp
            do l = 1, h%nt
                y(l, :) = IU*k*y(l, :)
            end do
        end if
        y = matmul(conjg(h%ft), matmul(y, conjg(transpose(h%fp))))
        df = reshape(real(y), [h%nt*h%np])/real(h%nt*h%np, dp)
    end function grid_derivative

    !> lambda from the weak form of J^s = 0 on the surface (VMEC's modes, the
    !> (0, 0) mode excluded): for every test mode e_k,
    !> int (B_theta d_phi e_k - B_phi d_theta e_k) = 0 with
    !> sqrt(g) B_theta = g_tt (iota - lambda_p) + g_tp (1 + lambda_t),
    !> sqrt(g) B_phi = g_tp (iota - lambda_p) + g_pp (1 + lambda_t).
    !> Returns max |lambda_rec - lambda| on the grid (rad), constants removed.
    real(dp) function lambda_defect(h, sv) result(eps)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), allocatable :: e(:, :), a(:, :), rhs(:, :), lrec(:), g(:, :)
        real(dp), allocatable :: rm(:), rn(:), sg(:)
        integer, allocatable :: mi(:), ni(:), dc(:)
        complex(dp), allocatable :: gh(:, :, :), et(:, :), ep(:, :)
        real(dp) :: ct(h%nt), st(h%nt), cp(h%np), sp(h%np), stt, stp, spp
        integer :: npt, nb, k, i, j, info, par, l, mm, nn

        npt = h%nt*h%np
        nb = 0
        do par = 1, merge(2, 1, h%lasym)
            do k = 1, h%mnmax
                if (h%xm(k) == 0 .and. h%xn(k) == 0) cycle
                nb = nb + 1
            end do
        end do
        ! Mode k: e = sin (par 1) or cos (par 2) of m theta - n N phi; its
        ! gradient is sg (m, -n N) times dc = cos (1) or sin (2).
        allocate (e(npt, nb), mi(nb), ni(nb), rm(nb), rn(nb), sg(nb), dc(nb))
        j = 0
        do par = 1, merge(2, 1, h%lasym)
            do k = 1, h%mnmax
                if (h%xm(k) == 0 .and. h%xn(k) == 0) cycle
                j = j + 1
                mi(j) = h%xm(k)
                ni(j) = h%xn(k)/h%nfp
                rm(j) = real(h%xm(k), dp)
                rn(j) = real(h%xn(k), dp)
                sg(j) = merge(1.0_dp, -1.0_dp, par == 1)
                dc(j) = par
                ct = cos(rm(j)*h%theta)
                st = sin(rm(j)*h%theta)
                cp = cos(rn(j)*h%phi)
                sp = sin(rn(j)*h%phi)
                do l = 1, h%np
                    if (par == 1) then
                        e((l - 1)*h%nt + 1:l*h%nt, j) = st*cp(l) - ct*sp(l)
                    else
                        e((l - 1)*h%nt + 1:l*h%nt, j) = ct*cp(l) + st*sp(l)
                    end if
                end do
            end do
        end do
        allocate (g(npt, 5))
        associate (x => sv%x)
            g(:, 1) = 1.0_dp/(x(1, :)*(x(3, :)*x(5, :) - x(2, :)*x(6, :)))
            g(:, 2) = g(:, 1)*(x(3, :)*x(4, :) + x(6, :)*x(7, :))
            g(:, 3) = g(:, 1)*(x(1, :)**2 + x(4, :)**2 + x(7, :)**2)
            g(:, 1) = g(:, 1)*(x(3, :)**2 + x(6, :)**2)
        end associate
        g(:, 4) = g(:, 1)*sv%iota + g(:, 2)
        g(:, 5) = g(:, 2)*sv%iota + g(:, 3)
        ! Grid sums of products of two modes are sums of the metric DFT at the
        ! sum and difference modes: O(npt (M + N)) instead of O(npt nb**2).
        mm = 2*h%mmax
        nn = 2*h%nmax
        allocate (et(-mm:mm, h%nt), ep(h%np, -nn:nn), gh(-mm:mm, -nn:nn, 5))
        do i = -mm, mm
            et(i, :) = exp(IU*real(i, dp)*h%theta)
        end do
        do i = -nn, nn
            ep(:, i) = exp(-IU*real(i*h%nfp, dp)*h%phi)
        end do
        do k = 1, 5
            gh(:, :, k) = matmul(et, matmul(cmplx(reshape(g(:, k), [h%nt, h%np]), &
                                                  0.0_dp, kind=dp), ep))
        end do
        allocate (a(nb, nb), rhs(nb, 1), lrec(npt))
        do j = 1, nb
            do i = 1, nb
                stt = pair_sum(1, i, j)
                stp = pair_sum(2, i, j)
                spp = pair_sum(3, i, j)
                a(i, j) = -sg(i)*sg(j)*(rn(i)*rn(j)*stt &
                                        + (rn(i)*rm(j) + rm(i)*rn(j))*stp &
                                        + rm(i)*rm(j)*spp)
            end do
            rhs(j, 1) = sg(j)*(rn(j)*single_sum(4, j) + rm(j)*single_sum(5, j))
        end do
        call dense_solve(a, rhs(:, 1), info)
        if (info /= 0) then
            eps = huge(1.0_dp)
            return
        end if
        lrec = matmul(e, rhs(:, 1))
        lrec = lrec - sv%lam
        eps = maxval(abs(lrec - sum(lrec)/real(npt, dp)))

    contains

        !> sum_p g_k D_i D_j, D = cos or sin of the mode phase (dc).
        real(dp) function pair_sum(kf, i, j) result(r)
            integer, intent(in) :: kf, i, j
            complex(dp) :: zs, zd

            zs = gh(mi(i) + mi(j), ni(i) + ni(j), kf)
            zd = gh(mi(i) - mi(j), ni(i) - ni(j), kf)
            if (dc(i) == 1 .and. dc(j) == 1) then
                r = 0.5_dp*real(zd + zs)
            else if (dc(i) == 2 .and. dc(j) == 2) then
                r = 0.5_dp*real(zd - zs)
            else if (dc(i) == 2) then
                r = 0.5_dp*aimag(zs + zd)
            else
                r = 0.5_dp*aimag(zs - zd)
            end if
        end function pair_sum

        !> sum_p g_k D_j.
        real(dp) function single_sum(kf, j) result(r)
            integer, intent(in) :: kf, j

            if (dc(j) == 1) then
                r = real(gh(mi(j), ni(j), kf))
            else
                r = aimag(gh(mi(j), ni(j), kf))
            end if
        end function single_sum
    end function lambda_defect

    !> Gaussian elimination with partial pivoting in a fixed operation order:
    !> Accelerate's dgesv rounds differently depending on the memory layout, so
    !> repeated checks of one file would not reproduce bitwise. info = k if the
    !> k-th pivot is zero.
    pure subroutine dense_solve(a, b, info)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: info
        real(dp) :: row(size(a, 2)), t
        integer :: n, k, i, p

        n = size(a, 1)
        info = 0
        do k = 1, n
            p = k - 1 + maxloc(abs(a(k:n, k)), 1)
            if (a(p, k) == 0.0_dp) then
                info = k
                return
            end if
            if (p /= k) then
                row = a(k, :)
                a(k, :) = a(p, :)
                a(p, :) = row
                t = b(k)
                b(k) = b(p)
                b(p) = t
            end if
            a(k + 1:n, k) = a(k + 1:n, k)/a(k, k)
            do i = k + 1, n
                a(k + 1:n, i) = a(k + 1:n, i) - a(k + 1:n, k)*a(k, i)
            end do
            b(k + 1:n) = b(k + 1:n) - a(k + 1:n, k)*b(k)
        end do
        do k = n, 1, -1
            b(k) = (b(k) - dot_product(a(k, k + 1:n), b(k + 1:n)))/a(k, k)
        end do
    end subroutine dense_solve

    !> Resonant residual at the rational surfaces iota = n N/m in the window
    !> (m <= mpol - 1, |n| <= ntor, m and n coprime): max over the harmonics
    !> (k m, k n), k m <= mnu, of |f_mn| of the sfl_jacobian_b2 route.
    subroutine rational_surfaces(h, c, q)
        type(boozer_hybrid_t), intent(in) :: h
        type(sfl_quality_config_t), intent(in) :: c
        type(sfl_quality_t), intent(inout) :: q
        integer, parameter :: NSCAN = 400, MAXRES = 512
        real(dp), allocatable :: coef(:, :), sg(:), ig(:), rs(:), rf(:)
        integer, allocatable :: rm(:), rn(:)
        type(hybrid_surface_t) :: sv
        real(dp) :: s0, s1, target, fbest, f
        integer :: m, n, i, nr, kk, mbest, nbest, nn

        nn = h%ns - 1
        allocate (coef(0:5, nn))
        coef = 0.0_dp
        coef(0, :) = h%iotas(1:nn)
        call spl_reg(5, nn, 1.0_dp/real(h%ns - 1, dp), coef)
        allocate (sg(0:NSCAN), ig(0:NSCAN))
        do i = 0, NSCAN
            sg(i) = (c%rho_window(1) + (c%rho_window(2) - c%rho_window(1)) &
                     *real(i, dp)/real(NSCAN, dp))**2
            ig(i) = iota_at(sg(i))
        end do
        nr = 0
        allocate (rs(MAXRES), rf(MAXRES), rm(MAXRES), rn(MAXRES))
        do m = 1, h%mmax
            do n = -h%nmax, h%nmax
                if (gcd(m, abs(n)) /= 1) cycle
                target = real(n*h%nfp, dp)/real(m, dp)
                do i = 1, NSCAN
                    if ((ig(i - 1) - target)*(ig(i) - target) > 0.0_dp) cycle
                    if (ig(i - 1) == target .and. i > 1) cycle
                    s0 = sg(i - 1)
                    s1 = sg(i)
                    call bisect(s0, s1, target)
                    call hybrid_surface(h, 0.5_dp*(s0 + s1), sv, with_spectrum=.false.)
                    fbest = 0.0_dp
                    mbest = m
                    nbest = n
                    do kk = 1, h%mnu/m
                        if (abs(kk*n) > h%nnu) exit
                        f = abs(sv%fmn(kk*m, kk*n))
                        if (f > fbest) then
                            fbest = f
                            mbest = kk*m
                            nbest = kk*n
                        end if
                    end do
                    if (nr == MAXRES) cycle
                    nr = nr + 1
                    rs(nr) = 0.5_dp*(s0 + s1)
                    rf(nr) = fbest
                    rm(nr) = mbest
                    rn(nr) = nbest
                end do
            end do
        end do
        allocate (q%res_s(nr), q%res_f(nr), q%res_m(nr), q%res_n(nr))
        q%res_s = rs(1:nr)
        q%res_f = rf(1:nr)
        q%res_m = rm(1:nr)
        q%res_n = rn(1:nr)

    contains

        real(dp) function iota_at(s) result(v)
            real(dp), intent(in) :: s
            real(dp) :: hs, x0, t
            integer :: j, k

            hs = 1.0_dp/real(h%ns - 1, dp)
            x0 = 0.5_dp*hs
            j = max(0, min(nn - 2, floor((s - x0)/hs)))
            t = s - (x0 + real(j, dp)*hs)
            v = 0.0_dp
            do k = 5, 0, -1
                v = coef(k, j + 1) + t*v
            end do
        end function iota_at

        subroutine bisect(a, b, y)
            real(dp), intent(inout) :: a, b
            real(dp), intent(in) :: y
            real(dp) :: mid, fa
            integer :: it

            fa = iota_at(a) - y
            do it = 1, 60
                mid = 0.5_dp*(a + b)
                if ((iota_at(mid) - y)*fa > 0.0_dp) then
                    a = mid
                    fa = iota_at(a) - y
                else
                    b = mid
                end if
            end do
        end subroutine bisect

        pure integer function gcd(i, j)
            integer, intent(in) :: i, j
            integer :: a, b, t

            a = i
            b = j
            do while (b /= 0)
                t = mod(a, b)
                a = b
                b = t
            end do
            gcd = a
        end function gcd
    end subroutine rational_surfaces

    subroutine aggregate(q)
        type(sfl_quality_t), intent(inout) :: q

        q%value(SFL_ERR_RESOLUTION) = q%eps_rad
        q%s_worst(SFL_ERR_RESOLUTION) = q%s_rad
        call median_of(q%eps_fb, SFL_ERR_FORCE_BALANCE)
        call median_of(q%eps_js, SFL_ERR_SURFACE_CURRENT)
        call max_of(q%eps_lam, q%s, SFL_ERR_LAMBDA)
        call median_of(q%d_fam, SFL_ERR_FAMILY)
        if (size(q%res_f) > 0) then
            call max_of(q%res_f, q%res_s, SFL_ERR_RESONANCE)
        else
            q%value(SFL_ERR_RESONANCE) = 0.0_dp
        end if

    contains

        !> Median over the check surfaces; the worst surface is the maximum.
        subroutine median_of(v, k)
            real(dp), intent(in) :: v(:)
            integer, intent(in) :: k
            real(dp) :: w(size(v)), t
            integer :: i, j, n

            n = size(v)
            w = v
            do i = 2, n
                t = w(i)
                j = i - 1
                do while (j >= 1)
                    if (w(j) <= t) exit
                    w(j + 1) = w(j)
                    j = j - 1
                end do
                w(j + 1) = t
            end do
            if (mod(n, 2) == 1) then
                q%value(k) = w((n + 1)/2)
            else
                q%value(k) = 0.5_dp*(w(n/2) + w(n/2 + 1))
            end if
            q%s_worst(k) = q%s(maxloc(v, 1))
        end subroutine median_of

        subroutine max_of(v, s, k)
            real(dp), intent(in) :: v(:), s(:)
            integer, intent(in) :: k

            q%value(k) = maxval(v)
            q%s_worst(k) = s(maxloc(v, 1))
        end subroutine max_of
    end subroutine aggregate

    !> Status and message: one line per failing metric (surface, cause, nearest
    !> resonance, action), the first failing metric sets the code.
    subroutine evaluate(q, c)
        type(sfl_quality_t), intent(inout) :: q
        type(sfl_quality_config_t), intent(in) :: c
        character(len=:), allocatable :: msg
        character(len=1024) :: line
        integer :: k, mr, nr
        real(dp) :: dist, io

        call append(msg, '')
        q%status = sfl_status_t()
        do k = 1, NMETRIC
            if (q%tol(k) < 0.0_dp) cycle
            if (q%value(k) <= q%tol(k)) cycle
            io = iota_near(q, q%s_worst(k))
            call nearest_resonance(q, io, mr, nr, dist)
            if (k == SFL_ERR_RESONANCE) then
                mr = q%res_m(maxloc(q%res_f, 1))
                nr = q%res_n(maxloc(q%res_f, 1))
                dist = 0.0_dp
            end if
            write (line, '(a, f6.4, a, f5.3, a, a, a, es8.2, a, es8.2, a, a, a, i0, a, i0, &
                   &a, es9.2, a)') &
                'sfl: inconsistent equilibrium '//trim(q%wout)//' at s = ', q%s_worst(k), &
                ' (rho ', sqrt(max(q%s_worst(k), 0.0_dp)), ', ns '//trim(itoa(q%ns))//'): ', &
                trim(sfl_metric_name(k))//trim(aggregation(k)), ' = ', q%value(k), &
                ' > tol ', q%tol(k), &
                ' [cause: ', trim(cause(k)), ']; nearest resonance (m,n) = (', mr, ',', nr, &
                '), iota m - n N = ', dist, '.'
            call append(msg, trim(line)//' '//trim(action(k))//new_line('a'))
            if (q%status%code == SFL_OK) then
                q%status%code = k
                q%status%s = q%s_worst(k)
                q%status%value = q%value(k)
                q%status%tol = q%tol(k)
                q%status%m_res = mr
                q%status%n_res = nr
                q%status%res_dist = dist
            end if
        end do
        if (q%status%code == SFL_OK) then
            call append(q%status%message, '')
            return
        end if
        call append(msg, 'Override: allow_poor_equilibrium = .true. with a reason '// &
                    '(logged and written to the output).')
        call append(q%status%message, msg)
        if (c%allow_poor_equilibrium) then
            q%status%overridden = .true.
            q%overridden = .true.
            write (error_unit, '(a)') 'sfl: OVERRIDE ('//trim(c%override_reason)//')'
            write (error_unit, '(a)') msg
        end if

    contains

        function aggregation(k) result(t)
            integer, intent(in) :: k
            character(len=24) :: t

            select case (k)
            case (SFL_ERR_FORCE_BALANCE, SFL_ERR_SURFACE_CURRENT, SFL_ERR_FAMILY)
                t = ' (median over rho)'
            case (SFL_ERR_RESOLUTION)
                t = ''
            case default
                t = ' (max over rho)'
            end select
        end function aggregation

        function cause(k) result(t)
            integer, intent(in) :: k
            character(len=64) :: t

            select case (k)
            case (SFL_ERR_RESOLUTION)
                t = 'radial resolution of R, Z'
            case (SFL_ERR_FORCE_BALANCE)
                t = 'radial force balance (resolution or convergence)'
            case (SFL_ERR_SURFACE_CURRENT)
                t = 'surface current J^s (angular resolution)'
            case (SFL_ERR_LAMBDA)
                t = 'lambda inconsistent with R, Z (J^s = 0 not satisfied)'
            case (SFL_ERR_FAMILY)
                t = 'Boozer covariant and Jacobian routes differ (J^s /= 0)'
            case default
                t = 'resonant residual at a rational surface'
            end select
        end function cause

        function action(k) result(t)
            integer, intent(in) :: k
            character(len=200) :: t
            character(len=200) :: conv

            write (conv, '(a, es9.2, a, es9.2, a)') '(fsqr+fsqz+fsql = ', &
                q%fsqr + q%fsqz + q%fsql, ', ftolv = ', q%ftolv, ')'
            select case (k)
            case (SFL_ERR_RESOLUTION)
                t = 'Raise ns.'
            case (SFL_ERR_FORCE_BALANCE, SFL_ERR_LAMBDA)
                t = 'Raise ns and converge to a smaller ftol '//trim(conv)//'.'
            case (SFL_ERR_SURFACE_CURRENT)
                t = 'Raise mpol and ntor (now '//trim(itoa(q%mpol))//', '// &
                    trim(itoa(q%ntor))//').'
            case default
                t = 'Raise ns, mpol and ntor and converge to a smaller ftol '// &
                    trim(conv)//'.'
            end select
        end function action
    end subroutine evaluate

    real(dp) function iota_near(q, s) result(io)
        type(sfl_quality_t), intent(in) :: q
        real(dp), intent(in) :: s
        integer :: k

        k = minloc(abs(q%s - s), 1)
        io = q%iota(k)
    end function iota_near

    subroutine nearest_resonance(q, io, mr, nr, dist)
        type(sfl_quality_t), intent(in) :: q
        real(dp), intent(in) :: io
        integer, intent(out) :: mr, nr
        real(dp), intent(out) :: dist
        integer :: m, n, mmax, nmax
        real(dp) :: d

        mmax = max(q%mpol - 1, 1)
        nmax = q%ntor
        dist = huge(1.0_dp)
        mr = 0
        nr = 0
        do m = 1, mmax
            do n = -nmax, nmax
                d = io*real(m, dp) - real(n*q%nfp, dp)
                if (abs(d) < abs(dist)) then
                    dist = d
                    mr = m
                    nr = n
                end if
            end do
        end do
    end subroutine nearest_resonance

    pure function itoa(i) result(t)
        integer, intent(in) :: i
        character(len=12) :: t

        write (t, '(i0)') i
    end function itoa

end module boozer_hybrid_quality
