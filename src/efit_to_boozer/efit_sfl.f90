!> Straight-field-line coordinates of an axisymmetric G-EQDSK equilibrium from
!> its primaries psi(R, Z) and F(psi), by flux-surface contours.
!>
!> Primaries: the psirz grid (tensor quintic spline, libneo spl_reg) and the
!> fpol column (quintic spline in x_file = (psi - simag)/(sibry - simag)).
!> qpsi, pprime, ffprim and pres are derived columns and enter only the
!> consistency checks. Gaussian units (cm, G, G cm^2 per radian) as in
!> efit_to_boozer.
!>
!> Construction. The magnetic axis is the Newton root of grad psi. A surface
!> psi_t is found on rays R = R_0 + rho c, Z = Z_0 + rho s at uniform geometric
!> angles theta_g (Newton on the psi spline). With psi oriented to grow outward
!> and f = |F|, the kernel K7 (gen/gen_efit_sfl.f90) gives the integrands
!> h = f rho/(R psi_rho) (mean q), bg = |grad psi|^2 rho/(R psi_rho) (mean |I|),
!> |B|, jc = R rho/psi_rho and jb = h + bg/f. With A the spectral antiderivative
!> on the theta_g grid (zero at theta_g = 0), the routes give nu = A g with
!>   covariant (Boozer, B_theta and B_phi flux functions):
!>       g = (bg - h <bg>/q)/(f + <bg>/q)
!>   jacobian_b2 (Boozer, J B^2 flux function): g = q jb/<jb> - h
!>   jacobian_const (Hamada, J flux function):  g = q jc/<jc> - h
!>   straight_current (Hamada II, current lines straight):
!>       g = (q a - h <a>)/(<a> - q sig_f F'),  a = -jc Delta* psi/R^2
!>   pest: g = 0,
!> and theta* = theta_g + A(h + g)/q, phi* = phi + nu, d theta*/d theta_g =
!> (h + g)/q. straight_current is R4 of the design, phi_H = phi_B + eta with
!> eta = -(beta - <beta>_B)/(G' + iota I'), in tangential form: in axisymmetry
!> Ampere's law gives d beta/d theta_B = I' - J_B (curl B).grad phi_B, so beta
!> needs no radial derivative, and the straightness of curl B in (theta*,
!> phi*), with (curl B).grad phi = -Delta* psi/R^2 and (curl B).grad theta_g =
!> sig_f F'/jc, is solved for g directly (same angles up to the gauge
!> nu(theta_g = 0) = 0). The divisor is |<a> - q sig_f F'| = q |G' + iota I'|
!> ds/|dpsi| and vanishes with p' (EFIT_SFL_ERR_NO_PRESSURE). In equilibrium the
!> route equals jacobian_const (d_family = d_H). Both Boozer routes agree
!> analytically in axisymmetry (J.grad psi
!> = 0 identically), so their difference measures round-off only. Spectra are
!> quadratures over theta_g with the angle-map Jacobian, without
!> re-interpolation, in the convention of fromefit.bc (nu(theta_g = 0) = 0, theta
!> along the poloidal motion of the field line for increasing phi).
!>
!> Surfaces are labelled by s, the toroidal flux normalized to psi_edge. q is
!> computed on n_q contours at x_j = (j - 1/2)/n_q, x = (psi - psi_axis)/(psi_edge
!> - psi_axis), splined (quintic) and integrated.
!>
!> Derivatives (efit_sfl_jvp, efit_sfl_vjp) are exact for the discrete map
!> with respect to psirz and fpol: implicit-function theorem for the axis, the
!> contours and the s labels, the generated kernel Jacobian, and hand-written
!> transposes of the spline rows and spectral operators.
!>
!> Checks (section 6.7 of the boozer-bench SFL design): eps_gs, the relative
!> L2 Grad-Shafranov residual of psi with the file's pprime and ffprim columns,
!> over the rho = sqrt(s) window; eps_q = max |q_contour/qpsi - 1| in the window;
!> eps_fb, the flux-surface averaged force balance with the file's p'; the
!> consistency of the pprime and ffprim columns with pres and fpol. eps_gs and
!> eps_q above tolerance are hard errors unless allow_poor_equilibrium is set.
module efit_sfl
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use geqdsk_tools, only: geqdsk_t, geqdsk_read, geqdsk_deinit
    use spl_three_to_five_sub, only: spl_reg
    implicit none
    private

    public :: efit_sfl_config_t, efit_sfl_status_t, efit_sfl_t, efit_sfl_surface_t, &
              efit_sfl_tangent_t, efit_sfl_diag_t
    public :: contour_t, efit_sfl_contour_at
    public :: efit_sfl_init, efit_sfl_transform, efit_sfl_jvp, efit_sfl_vjp, &
              efit_sfl_raise, efit_sfl_diag_write, efit_sfl_gs_map, efit_sfl_s_of_x, &
              efit_sfl_perturb, &
              efit_sfl_x_of_s, efit_sfl_symflux

    integer, parameter, public :: EFIT_SFL_COVARIANT = 1, EFIT_SFL_JACOBIAN_B2 = 2, &
                                  EFIT_SFL_JACOBIAN_CONST = 3, EFIT_SFL_STRAIGHT_CURRENT = 4, &
                                  EFIT_SFL_PEST = 5
    integer, parameter, public :: EFIT_SFL_OK = 0, EFIT_SFL_ERR_FORCE_BALANCE = 2, &
                                  EFIT_SFL_ERR_NO_PRESSURE = 4, &
                                  EFIT_SFL_ERR_GS = 5, &
                                  EFIT_SFL_ERR_Q = 6, EFIT_SFL_ERR_INPUT = 7
    real(dp), parameter :: pi = 3.14159265358979323846_dp, twopi = 2.0_dp*pi

    type :: efit_sfl_config_t
        integer :: route = EFIT_SFL_COVARIANT
        integer :: n_theta = 1024             ! geometric angles per contour (even)
        integer :: n_q = 128                  ! q nodes for the toroidal flux
        real(dp) :: psi_edge = huge(1.0_dp)   ! file psi gauge and units; default sibry
        real(dp) :: psi_edge_rel = huge(1.0_dp)   ! psi_edge - simag (efit_to_boozer psimax)
        real(dp) :: rho_window(2) = [0.2_dp, 0.95_dp]   ! rho = sqrt(s)
        ! 3 x the maximum over the pass set (Solov'ev, CHEASE circular, FreeGS LSN
        ! and MAST with the standard profile grid at 129-513 points), rounded up;
        ! every FreeGS file written with the shifted profile grid fails
        ! (boozer-bench reports/efit-sfl.md, calibration table). Window medians.
        real(dp) :: tol_gs = 3.0e-3_dp, tol_fb = 2.0e-5_dp, tol_q = 3.0e-3_dp
        logical :: allow_poor_equilibrium = .false.
        character(len=256) :: override_reason = ''
    end type

    type :: efit_sfl_status_t
        integer :: code = EFIT_SFL_OK
        logical :: overridden = .false.
        integer :: isurf = 0
        real(dp) :: s = -1.0_dp, value = 0.0_dp, tol = 0.0_dp
        character(len=:), allocatable :: message
    end type

    type :: efit_sfl_diag_t
        real(dp), allocatable :: x(:), s(:), q(:), qpsi(:), eps_q(:), eps_gs(:), &
                                 eps_fb(:), curr_i(:), curr_g(:), pp_file(:), pp_fb(:)
        logical, allocatable :: in_window(:)
        real(dp) :: eps_gs_window = 0, eps_q_window = 0, eps_fb_median = 0, eps_gs_pres_fpol = 0, &
                    eps_gs_median = 0, &
                    eps_pp = 0, eps_ff = 0, eps_axis = 0, gs_factor = 1
        integer :: exp_bpol = 0
        logical :: override = .false.
        type(efit_sfl_status_t) :: status
    end type

    !> One surface of the output: spectra in the fromefit.bc convention
    !> (Gaussian units), cos/sin pairs of R, Z, v = -nu/(2 pi) and |B| over
    !> theta*, plus flux functions. Also the type of tangents and adjoints.
    type :: efit_sfl_surface_t
        real(dp) :: s = 0, x = 0, psi = 0, iota = 0, curr_i = 0, curr_g = 0, &
                    sqrtg00 = 0, pp_file = 0, pp_fb = 0, eps_fb = 0, d_family = -1, &
                    eps_theta = 0
        real(dp), allocatable :: rmnc(:), rmns(:), zmnc(:), zmns(:), vmnc(:), &
                                 vmns(:), bmnc(:), bmns(:)
    end type

    type :: efit_sfl_tangent_t
        real(dp), allocatable :: psirz(:, :), fpol(:)
    end type

    !> Points of one contour with the primal data the derivatives need.
    type :: contour_t
        real(dp) :: x = 0, psit = 0, f = 0, fpsi = 0, xf = 0
        real(dp), allocatable :: rho(:), r(:), z(:), praw(:), pzraw(:), prho(:), &
                                 prr(:), prz(:), pzz(:), out(:, :), jac(:, :, :), &
                                 ca(:), sa(:)
    end type

    type :: efit_sfl_t
        type(efit_sfl_config_t) :: cfg
        type(geqdsk_t) :: g
        integer :: nw = 0, nh = 0, nt = 0, nq = 0
        real(dp) :: r1 = 0, z1 = 0, hr = 0, hz = 0, hx = 0
        real(dp), allocatable :: cpsi(:, :, :, :), wr(:, :, :), wz(:, :, :)
        real(dp), allocatable :: cf(:, :), wf(:, :, :), cpp(:, :), cffp(:, :), &
                                 cpres(:, :), cqpsi(:, :), wq(:, :, :)
        real(dp) :: r0 = 0, z0 = 0, psi_ax = 0, hess(2, 2) = 0, psi_edge = 0
        real(dp) :: sig_psi = 1, sig_f = 1, sig_theta = 1, psitor_max = 0, phi1 = 0
        real(dp) :: rminor = 0
        real(dp), allocatable :: aop(:, :)
        type(contour_t), allocatable :: qn(:)
        real(dp), allocatable :: xq(:), qq(:), iq(:), cq(:, :), ciq(:, :)
        type(efit_sfl_diag_t) :: diag
    end type

contains

    include '../generated/efit_sfl_leaves.inc'

    ! ---------------------------------------------------------------- setup

    subroutine efit_sfl_init(self, gfile, cfg, status)
        type(efit_sfl_t), intent(inout) :: self
        character(*), intent(in) :: gfile
        type(efit_sfl_config_t), intent(in) :: cfg
        type(efit_sfl_status_t), intent(out), optional :: status
        type(efit_sfl_status_t) :: st
        real(dp) :: c

        self%cfg = cfg
        call geqdsk_read(self%g, gfile)
        call setup(self, st)
        if (st%code == EFIT_SFL_OK) then
            c = gs_fit_factor(self)
            ! psi in Wb (COCOS 11-18) instead of Wb/rad: Delta* psi = (2 pi)^2 x rhs.
            if (abs(log(abs(c)/(twopi**2))) < abs(log(abs(c)))) then
                self%g%psirz = self%g%psirz/twopi
                self%g%simag = self%g%simag/twopi
                self%g%sibry = self%g%sibry/twopi
                self%g%pprime = self%g%pprime*twopi
                self%g%ffprim = self%g%ffprim*twopi
                if (self%cfg%psi_edge < huge(1.0_dp)) &
                    self%cfg%psi_edge = self%cfg%psi_edge/twopi
                if (self%cfg%psi_edge_rel < huge(1.0_dp)) &
                    self%cfg%psi_edge_rel = self%cfg%psi_edge_rel/twopi
                call setup(self, st)
                self%diag%exp_bpol = 1
            end if
            self%diag%gs_factor = gs_fit_factor(self)
        end if
        if (st%code == EFIT_SFL_OK) call run_checks(self, st)
        self%diag%status = st
        call finish_status(self, st, status)
    end subroutine efit_sfl_init

    !> State with the primaries psirz + h dpsirz and fpol + h dfpol (finite
    !> differences in memory; checks are not rerun).
    subroutine efit_sfl_perturb(self, tan, h, pert)
        type(efit_sfl_t), intent(in) :: self
        type(efit_sfl_tangent_t), intent(in) :: tan
        real(dp), intent(in) :: h
        type(efit_sfl_t), intent(out) :: pert
        type(efit_sfl_status_t) :: st

        pert = self
        pert%g%psirz = self%g%psirz + h*tan%psirz
        pert%g%fpol = self%g%fpol + h*tan%fpol
        call setup(pert, st)
        call efit_sfl_raise(st)
    end subroutine efit_sfl_perturb

    subroutine setup(self, st)
        type(efit_sfl_t), intent(inout) :: self
        type(efit_sfl_status_t), intent(inout) :: st

        st = efit_sfl_status_t()
        self%nw = self%g%nw
        self%nh = self%g%nh
        self%r1 = self%g%R_eqd(1)
        self%z1 = self%g%Z_eqd(1)
        self%hr = self%g%R_eqd(2) - self%g%R_eqd(1)
        self%hz = self%g%Z_eqd(2) - self%g%Z_eqd(1)
        self%hx = 1.0_dp/real(self%nw - 1, dp)
        call spline_rows(self%nw, self%hr, self%wr)
        call spline_rows(self%nh, self%hz, self%wz)
        call spline_rows(self%nw, self%hx, self%wf)
        call psi_coefficients(self, self%g%psirz, self%cpsi)
        call profile(self%g%fpol, self%cf)
        call profile(self%g%pprime, self%cpp)
        call profile(self%g%ffprim, self%cffp)
        call profile(self%g%pres, self%cpres)
        call profile(self%g%qpsi, self%cqpsi)

        self%nt = self%cfg%n_theta
        if (mod(self%nt, 2) /= 0 .or. self%nt < 16) then
            call set_status(st, EFIT_SFL_ERR_INPUT, 0, -1.0_dp, real(self%nt, dp), 16.0_dp, &
                            'n_theta must be even and >= 16')
            return
        end if
        call find_axis(self, st)
        if (st%code /= EFIT_SFL_OK) return
        self%psi_edge = self%g%sibry
        if (self%cfg%psi_edge < huge(1.0_dp)) self%psi_edge = self%cfg%psi_edge
        if (self%cfg%psi_edge_rel < huge(1.0_dp)) &
            self%psi_edge = self%g%simag + self%cfg%psi_edge_rel
        self%sig_psi = sign(1.0_dp, self%psi_edge - self%psi_ax)
        self%sig_f = sign(1.0_dp, self%g%fpol(1))
        ! Legacy sigma = sign(B_Z B_phi) at the outboard midplane.
        self%sig_theta = self%sig_psi*self%sig_f
        if (allocated(self%aop)) deallocate (self%aop)
        call antiderivative_matrix(self%nt, self%aop)
        call q_nodes(self, st)
    contains
        subroutine profile(y, c)
            real(dp), intent(in) :: y(:)
            real(dp), allocatable, intent(out) :: c(:, :)
            allocate (c(0:5, self%nw))
            c = 0.0_dp
            c(0, :) = y
            call spl_reg(5, self%nw, self%hx, c)
        end subroutine profile
    end subroutine setup

    !> Rows w(i, a, k) = coefficient k at node i of the quintic spline of the
    !> unit vector e_a, so that the spline is linear in the data with these rows.
    subroutine spline_rows(n, h, w)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), allocatable, intent(out) :: w(:, :, :)
        real(dp), allocatable :: c(:, :)
        integer :: a, k

        allocate (w(n, n, 0:5), c(0:5, n))
        do a = 1, n
            c = 0.0_dp
            c(0, a) = 1.0_dp
            call spl_reg(5, n, h, c)
            do k = 0, 5
                w(:, a, k) = c(k, :)
            end do
        end do
    end subroutine spline_rows

    !> Tensor quintic spline coefficients c(k, l, i, j) of psi(R, Z).
    subroutine psi_coefficients(self, psi, c)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: psi(:, :)
        real(dp), allocatable, intent(out) :: c(:, :, :, :)
        real(dp), allocatable :: cr(:, :), cz(:, :), tmp(:, :, :)
        integer :: i, j, k

        allocate (c(0:5, 0:5, self%nw, self%nh), tmp(0:5, self%nw, self%nh))
        allocate (cr(0:5, self%nw), cz(0:5, self%nh))
        do j = 1, self%nh
            cr = 0.0_dp
            cr(0, :) = psi(:, j)
            call spl_reg(5, self%nw, self%hr, cr)
            tmp(:, :, j) = cr
        end do
        do i = 1, self%nw
            do k = 0, 5
                cz = 0.0_dp
                cz(0, :) = tmp(k, i, :)
                call spl_reg(5, self%nh, self%hz, cz)
                c(k, :, i, :) = cz
            end do
        end do
    end subroutine psi_coefficients

    !> psi and its first and second derivatives from coefficients c at (R, Z).
    pure subroutine psi_eval(self, c, r, z, v)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: c(0:, 0:, :, :), r, z
        real(dp), intent(out) :: v(6)        ! psi, R, Z, RR, RZ, ZZ
        real(dp) :: t, u, a0(0:5), a1(0:5), a2(0:5)
        integer :: i, j, k, l

        call locate(self, r, z, i, j, t, u)
        a0 = 0.0_dp
        a1 = 0.0_dp
        a2 = 0.0_dp
        do k = 0, 5
            do l = 5, 0, -1
                a0(k) = c(k, l, i, j) + u*a0(k)
            end do
            do l = 5, 1, -1
                a1(k) = real(l, dp)*c(k, l, i, j) + u*a1(k)
            end do
            do l = 5, 2, -1
                a2(k) = real(l*(l - 1), dp)*c(k, l, i, j) + u*a2(k)
            end do
        end do
        v = 0.0_dp
        do k = 5, 0, -1
            v(1) = a0(k) + t*v(1)
            v(3) = a1(k) + t*v(3)
            v(6) = a2(k) + t*v(6)
        end do
        do k = 5, 1, -1
            v(2) = real(k, dp)*a0(k) + t*v(2)
            v(5) = real(k, dp)*a1(k) + t*v(5)
        end do
        do k = 5, 2, -1
            v(4) = real(k*(k - 1), dp)*a0(k) + t*v(4)
        end do
    end subroutine psi_eval

    !> Transpose of psi_eval for the value and the first (vb(1:3)) and, if
    !> size(vb) = 6, the second derivatives.
    pure subroutine psi_eval_t(self, r, z, vb, cb)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: r, z, vb(:)
        real(dp), intent(inout) :: cb(:, :, 0:, 0:)
        real(dp) :: t, u, tk(0:5), uk(0:5), dtk(0:5), duk(0:5), ddtk(0:5), dduk(0:5)
        integer :: i, j, k, l

        call locate(self, r, z, i, j, t, u)
        tk(0) = 1.0_dp
        uk(0) = 1.0_dp
        dtk(0) = 0.0_dp
        duk(0) = 0.0_dp
        ddtk(0:1) = 0.0_dp
        dduk(0:1) = 0.0_dp
        do k = 1, 5
            tk(k) = tk(k - 1)*t
            uk(k) = uk(k - 1)*u
            dtk(k) = real(k, dp)*tk(k - 1)
            duk(k) = real(k, dp)*uk(k - 1)
        end do
        do k = 2, 5
            ddtk(k) = real(k*(k - 1), dp)*tk(k - 2)
            dduk(k) = real(k*(k - 1), dp)*uk(k - 2)
        end do
        do l = 0, 5
            do k = 0, 5
                cb(i, j, k, l) = cb(i, j, k, l) + vb(1)*tk(k)*uk(l) + vb(2)*dtk(k)*uk(l) &
                                 + vb(3)*tk(k)*duk(l)
                if (size(vb) == 6) cb(i, j, k, l) = cb(i, j, k, l) + vb(4)*ddtk(k)*uk(l) &
                    + vb(5)*dtk(k)*duk(l) + vb(6)*tk(k)*dduk(l)
            end do
        end do
    end subroutine psi_eval_t

    !> Third derivatives (RRR, RRZ, RZZ, ZZZ) of psi from coefficients c.
    pure subroutine psi_eval3(self, c, r, z, v)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: c(0:, 0:, :, :), r, z
        real(dp), intent(out) :: v(4)
        real(dp) :: t, u, a(0:5, 0:3)
        integer :: i, j, k, l, d

        call locate(self, r, z, i, j, t, u)
        a = 0.0_dp
        do k = 0, 5
            do d = 0, 3
                do l = 5, d, -1
                    a(k, d) = real(fall(l, d), dp)*c(k, l, i, j) + u*a(k, d)
                end do
            end do
        end do
        v = 0.0_dp
        do k = 5, 0, -1
            if (k >= 3) v(1) = real(fall(k, 3), dp)*a(k, 0) + t*v(1)
            if (k >= 2) v(2) = real(fall(k, 2), dp)*a(k, 1) + t*v(2)
            if (k >= 1) v(3) = real(k, dp)*a(k, 2) + t*v(3)
            v(4) = a(k, 3) + t*v(4)
        end do
    contains
        pure integer function fall(n, m)
            integer, intent(in) :: n, m
            integer :: p
            fall = 1
            do p = 0, m - 1
                fall = fall*(n - p)
            end do
        end function fall
    end subroutine psi_eval3

    pure subroutine locate(self, r, z, i, j, t, u)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        integer, intent(out) :: i, j
        real(dp), intent(out) :: t, u

        i = max(1, min(self%nw - 1, floor((r - self%r1)/self%hr) + 1))
        j = max(1, min(self%nh - 1, floor((z - self%z1)/self%hz) + 1))
        t = r - (self%r1 + real(i - 1, dp)*self%hr)
        u = z - (self%z1 + real(j - 1, dp)*self%hz)
    end subroutine locate

    pure logical function inside(self, r, z)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        inside = r >= self%r1 .and. r <= self%r1 + real(self%nw - 1, dp)*self%hr .and. &
                 z >= self%z1 .and. z <= self%z1 + real(self%nh - 1, dp)*self%hz
    end function inside

    !> Value and x-derivative of a profile spline at x_file (clamped to [0, 1],
    !> constant beyond as in field_divB0).
    pure subroutine prof_eval(self, c, x, v, dv)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: c(0:, :), x
        real(dp), intent(out) :: v, dv
        real(dp) :: t, xc
        integer :: i, k

        xc = max(0.0_dp, min(1.0_dp, x))
        i = max(1, min(self%nw - 1, floor(xc/self%hx) + 1))
        t = xc - real(i - 1, dp)*self%hx
        v = 0.0_dp
        dv = 0.0_dp
        do k = 5, 0, -1
            v = c(k, i) + t*v
        end do
        do k = 5, 1, -1
            dv = real(k, dp)*c(k, i) + t*dv
        end do
        if (x /= xc) dv = 0.0_dp
    end subroutine prof_eval

    pure real(dp) function xfile(self, psi)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: psi
        xfile = (psi - self%g%simag)/(self%g%sibry - self%g%simag)
    end function xfile

    subroutine find_axis(self, st)
        type(efit_sfl_t), intent(inout) :: self
        type(efit_sfl_status_t), intent(inout) :: st
        real(dp) :: v(6), det, dr, dz
        integer :: it

        self%r0 = self%g%rmaxis
        self%z0 = self%g%zmaxis
        do it = 1, 100
            call psi_eval(self, self%cpsi, self%r0, self%z0, v)
            det = v(4)*v(6) - v(5)**2
            dr = -(v(6)*v(2) - v(5)*v(3))/det
            dz = -(-v(5)*v(2) + v(4)*v(3))/det
            self%r0 = self%r0 + dr
            self%z0 = self%z0 + dz
            if (abs(dr) + abs(dz) < 1.0e-12_dp*(self%hr + self%hz)) exit
        end do
        call psi_eval(self, self%cpsi, self%r0, self%z0, v)
        self%psi_ax = v(1)
        self%hess = reshape([v(4), v(5), v(5), v(6)], [2, 2])
        if (v(4)*v(6) - v(5)**2 <= 0.0_dp .or. .not. inside(self, self%r0, self%z0)) &
            call set_status(st, EFIT_SFL_ERR_INPUT, 0, 0.0_dp, 0.0_dp, 0.0_dp, &
                            'no magnetic axis (extremum of psi) near rmaxis, zmaxis')
        self%diag%eps_axis = abs(xfile(self, self%psi_ax))
    end subroutine find_axis

    ! ------------------------------------------------------------- contours

    !> Contour psi = psi_t on rays at uniform geometric angles theta_g (oriented
    !> with the poloidal motion of the field line), Newton on the psi spline.
    !> The trapezoidal rule in theta_g converges geometrically; near an X-point
    !> n_theta must grow (efit_sfl_surface_t%eps_theta reports the half-grid
    !> difference of the spectra).
    subroutine contour(self, x, ct, st)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        type(contour_t), intent(out) :: ct
        type(efit_sfl_status_t), intent(inout) :: st
        real(dp) :: f, fp, th(self%nt)
        integer :: j, n

        n = self%nt
        allocate (ct%rho(n), ct%r(n), ct%z(n), ct%praw(n), ct%pzraw(n), ct%prho(n), &
                  ct%prr(n), ct%prz(n), ct%pzz(n), ct%out(5, n), ct%jac(5, 5, n), &
                  ct%ca(n), ct%sa(n))
        ct%x = x
        ct%psit = self%psi_ax + x*(self%psi_edge - self%psi_ax)
        ct%xf = xfile(self, ct%psit)
        call prof_eval(self, self%cf, ct%xf, f, fp)
        ct%f = abs(f)
        ct%fpsi = fp/(self%g%sibry - self%g%simag)
        do j = 1, n
            th(j) = twopi*real(j - 1, dp)/real(n, dp)
        end do
        call rays(self, ct, th, st)
    end subroutine contour

    subroutine rays(self, ct, th, st)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(inout) :: ct
        real(dp), intent(in) :: th(:)
        type(efit_sfl_status_t), intent(inout) :: st
        real(dp) :: v(6), target, rho, lo, hi, fl, fp, hc, e(2), step
        integer :: j, it

        target = abs(ct%psit - self%psi_ax)
        step = 0.5_dp*min(self%hr, self%hz)
        do j = 1, self%nt
            e = [cos(self%sig_theta*th(j)), sin(self%sig_theta*th(j))]
            hc = self%sig_psi*dot_product(e, matmul(self%hess, e))
            rho = sqrt(2.0_dp*target/max(hc, tiny(1.0_dp)))
            lo = 0.0_dp
            hi = -1.0_dp
            do it = 1, 200
                if (.not. inside(self, self%r0 + rho*e(1), self%z0 + rho*e(2))) then
                    hi = rho
                    rho = 0.5_dp*(lo + hi)
                    cycle
                end if
                call psi_eval(self, self%cpsi, self%r0 + rho*e(1), self%z0 + rho*e(2), v)
                fl = self%sig_psi*(v(1) - self%psi_ax) - target
                fp = self%sig_psi*(e(1)*v(2) + e(2)*v(3))
                if (fl < 0.0_dp) then
                    lo = max(lo, rho)
                else
                    hi = merge(rho, min(hi, rho), hi < 0.0_dp)
                end if
                if (abs(fl) <= 1.0e-14_dp*max(target, abs(self%psi_edge - self%psi_ax)) &
                    .and. fp > 0.0_dp) exit
                if (fp > 0.0_dp) then
                    rho = rho - fl/fp
                else
                    rho = -1.0_dp
                end if
                if (hi < 0.0_dp) then
                    if (rho <= lo) rho = lo + step
                else if (rho <= lo .or. rho >= hi) then
                    rho = 0.5_dp*(lo + hi)
                end if
            end do
            call psi_eval(self, self%cpsi, self%r0 + rho*e(1), self%z0 + rho*e(2), v)
            ct%ca(j) = e(1)
            ct%sa(j) = e(2)
            ct%rho(j) = rho
            ct%r(j) = self%r0 + rho*e(1)
            ct%z(j) = self%z0 + rho*e(2)
            ct%praw(j) = v(2)
            ct%pzraw(j) = v(3)
            ct%prho(j) = e(1)*v(2) + e(2)*v(3)
            ct%prr(j) = v(4)
            ct%prz(j) = v(5)
            ct%pzz(j) = v(6)
            if (self%sig_psi*ct%prho(j) <= 0.0_dp .or. it > 200) then
                call set_status(st, EFIT_SFL_ERR_INPUT, 0, ct%x, 0.0_dp, 0.0_dp, &
                    'flux surface is not star-shaped about the axis or leaves the grid; '// &
                    'lower psi_edge (psimax) below the separatrix')
                return
            end if
            call kernel(rho, e(1), e(2), ct%r(j), self%sig_psi*v(2), self%sig_psi*v(3), &
                        ct%f, ct%out(:, j), ct%jac(:, :, j))
        end do
    end subroutine rays

    !> Contour at x for diagnostics and tests.
    subroutine efit_sfl_contour_at(self, x, ct)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        type(contour_t), intent(out) :: ct
        type(efit_sfl_status_t) :: st
        call contour(self, x, ct, st)
        call efit_sfl_raise(st)
    end subroutine efit_sfl_contour_at

    !> K7 with the Jacobian columns (rho, R, psi_R, psi_Z, f).
    pure subroutine kernel(rho, c, s, r, pr, pz, f, y, jac)
        real(dp), intent(in) :: rho, c, s, r, pr, pz, f
        real(dp), intent(out) :: y(5), jac(5, 5)
        real(dp) :: dc(5), ds(5)

        call efit_sfl_contour(rho, c, s, r, pr, pz, f, &
            y(1), jac(1, 1), dc(1), ds(1), jac(1, 2), jac(1, 3), jac(1, 4), jac(1, 5), &
            y(2), jac(2, 1), dc(2), ds(2), jac(2, 2), jac(2, 3), jac(2, 4), jac(2, 5), &
            y(3), jac(3, 1), dc(3), ds(3), jac(3, 2), jac(3, 3), jac(3, 4), jac(3, 5), &
            y(4), jac(4, 1), dc(4), ds(4), jac(4, 2), jac(4, 3), jac(4, 4), jac(4, 5), &
            y(5), jac(5, 1), dc(5), ds(5), jac(5, 2), jac(5, 3), jac(5, 4), jac(5, 5))
    end subroutine kernel

    !> Periodic antiderivative on the uniform grid theta_j = 2 pi (j-1)/n,
    !> trigonometric interpolant without the mean and Nyquist modes, gauge
    !> A(theta_1) = 0: A(j, l) = (2/n) sum_k (sin(k(theta_j - theta_l)) + sin(k theta_l))/k.
    subroutine antiderivative_matrix(n, a)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: a(:, :)
        real(dp) :: sk(0:n - 1)
        integer :: m, k, j, l

        do m = 0, n - 1
            sk(m) = 0.0_dp
            do k = 1, n/2 - 1
                sk(m) = sk(m) + sin(twopi*real(k*m, dp)/real(n, dp))/real(k, dp)
            end do
        end do
        allocate (a(n, n))
        do l = 1, n
            do j = 1, n
                a(j, l) = 2.0_dp/real(n, dp)*(sk(modulo(j - l, n)) + sk(l - 1))
            end do
        end do
    end subroutine antiderivative_matrix

    ! -------------------------------------------------- q nodes and s labels

    subroutine q_nodes(self, st)
        type(efit_sfl_t), intent(inout) :: self
        type(efit_sfl_status_t), intent(inout) :: st
        integer :: j, nq

        nq = self%cfg%n_q
        self%nq = nq
        if (allocated(self%qn)) deallocate (self%qn, self%xq, self%qq, self%iq)
        allocate (self%qn(nq), self%xq(nq), self%qq(nq), self%iq(nq))
        do j = 1, nq
            self%xq(j) = xmap((real(j, dp) - 0.5_dp)/real(nq, dp))
            call contour(self, self%xq(j), self%qn(j), st)
            if (st%code /= EFIT_SFL_OK) return
            self%qq(j) = sum(self%qn(j)%out(1, :))/real(self%nt, dp)
            self%iq(j) = self%sig_f*sum(self%qn(j)%out(2, :))/real(self%nt, dp)
        end do
        call spline_rows(nq, 1.0_dp/real(nq, dp), self%wq)
        ! Phi integrand in u: q dx/du (nodes clustered towards psi_edge)
        call spline1(self%qq*[(dxdu((real(j, dp) - 0.5_dp)/real(nq, dp)), j=1, nq)], self%cq)
        call spline1(self%iq, self%ciq)
        self%phi1 = qint(self, self%cq, 1.0_dp)
        self%psitor_max = (self%psi_edge - self%psi_ax)*self%phi1
        block
            type(contour_t) :: edge
            call contour(self, 1.0_dp, edge, st)
            if (st%code /= EFIT_SFL_OK) return
            self%rminor = sqrt(sum(edge%rho**2)/real(self%nt, dp))
        end block
    contains
        subroutine spline1(y, c)
            real(dp), intent(in) :: y(:)
            real(dp), allocatable, intent(out) :: c(:, :)
            allocate (c(0:5, nq))
            c = 0.0_dp
            c(0, :) = y
            call spl_reg(5, nq, 1.0_dp/real(nq, dp), c)
        end subroutine spline1
    end subroutine q_nodes

    pure subroutine qseg(self, x, i, t)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(out) :: i
        real(dp), intent(out) :: t
        real(dp) :: h

        real(dp) :: u
        h = 1.0_dp/real(self%nq, dp)
        u = umap(x)
        i = max(1, min(self%nq - 1, floor(u/h - 0.5_dp) + 1))
        t = u - (real(i, dp) - 0.5_dp)*h
    end subroutine qseg

    !> Radial node variable u in [0, 1]: x = (1 + kappa) u - kappa u^2, kappa = 0.9,
    !> so that the q nodes are ten times denser in x at psi_edge than at the axis
    !> (q grows logarithmically towards a separatrix).
    pure real(dp) function xmap(u)
        real(dp), intent(in) :: u
        xmap = 1.9_dp*u - 0.9_dp*u**2
    end function xmap

    pure real(dp) function dxdu(u)
        real(dp), intent(in) :: u
        dxdu = 1.9_dp - 1.8_dp*u
    end function dxdu

    pure real(dp) function umap(x)
        real(dp), intent(in) :: x
        umap = (1.9_dp - sqrt(max(1.9_dp**2 - 3.6_dp*x, 0.0_dp)))/1.8_dp
    end function umap

    !> q(x) and dq/dx from the Phi-integrand spline, I(x) and dI/dx.
    pure subroutine qfun(self, x, q, dqdx)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp), intent(out) :: q, dqdx
        real(dp) :: g, dg, u
        u = umap(x)
        call qval(self, self%cq, x, g, dg)
        q = g/dxdu(u)
        dqdx = (dg/dxdu(u) + 1.8_dp*g/dxdu(u)**2)/dxdu(u)
    end subroutine qfun

    !> Integral from 0 to x of the node spline with coefficients c (end
    !> segments extrapolated to 0 and 1).
    pure real(dp) function qint(self, c, x)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: c(0:, :), x
        real(dp) :: h, t
        integer :: i, m

        h = 1.0_dp/real(self%nq, dp)
        call qseg(self, x, i, t)
        qint = pint(c(:, 1), -0.5_dp*h, 0.0_dp)
        do m = 1, i - 1
            qint = qint + pint(c(:, m), 0.0_dp, h)
        end do
        qint = qint + pint(c(:, i), 0.0_dp, t)
    end function qint

    pure real(dp) function pint(c, a, b)
        real(dp), intent(in) :: c(0:5), a, b
        integer :: k
        pint = 0.0_dp
        do k = 0, 5
            pint = pint + c(k)*(b**(k + 1) - a**(k + 1))/real(k + 1, dp)
        end do
    end function pint

    pure subroutine qval(self, c, x, v, dv)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: c(0:, :), x
        real(dp), intent(out) :: v, dv
        real(dp) :: t
        integer :: i, k

        call qseg(self, x, i, t)
        v = 0.0_dp
        dv = 0.0_dp
        do k = 5, 0, -1
            v = c(k, i) + t*v
        end do
        do k = 5, 1, -1
            dv = real(k, dp)*c(k, i) + t*dv
        end do
    end subroutine qval

    !> Weights of qint(x) and of the spline value at x on the node data.
    pure subroutine qweights(self, x, wint, wval)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp), intent(out) :: wint(:), wval(:)
        real(dp) :: c(0:5, self%nq), dv
        integer :: a, k

        do a = 1, self%nq
            do k = 0, 5
                c(k, :) = self%wq(:, a, k)
            end do
            wint(a) = qint(self, c, x)
            call qval(self, c, x, wval(a), dv)
        end do
    end subroutine qweights

    real(dp) function efit_sfl_s_of_x(self, x)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        efit_sfl_s_of_x = qint(self, self%cq, x)/self%phi1
    end function efit_sfl_s_of_x

    real(dp) function efit_sfl_x_of_s(self, s)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: s
        real(dp) :: x, q, dq
        integer :: it

        x = s
        do it = 1, 100
            call qfun(self, x, q, dq)
            x = x - (qint(self, self%cq, x)/self%phi1 - s)*self%phi1/q
            x = max(0.0_dp, min(1.0_dp, x))
            if (abs(qint(self, self%cq, x)/self%phi1 - s) < 1.0e-15_dp) exit
        end do
        efit_sfl_x_of_s = x
    end function efit_sfl_x_of_s

    ! ------------------------------------------------------------- spectra

    !> Surface map: route integrand g, nu = A g, theta* - theta_g = A(h + g)/q,
    !> Jacobian (h + g)/q, spectra of R, Z, v = -nu/(2 pi), |B|.
    subroutine surface_map(self, route, mpol, ct, sv, g, u, jb)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: route, mpol
        type(contour_t), intent(in) :: ct
        type(efit_sfl_surface_t), intent(inout) :: sv
        real(dp), intent(out) :: g(:), u(:), jb(:)
        real(dp) :: q, mb, nu(self%nt), fld(self%nt, 4)

        q = sum(ct%out(1, :))/real(self%nt, dp)
        mb = sum(ct%out(2, :))/real(self%nt, dp)
        call route_integrand(self, route, ct, g)
        nu = matmul(self%aop, g)
        u = matmul(self%aop, ct%out(1, :) + g)/q
        jb = (ct%out(1, :) + g)/q
        fld(:, 1) = ct%r
        fld(:, 2) = ct%z
        fld(:, 3) = -nu/twopi
        fld(:, 4) = ct%out(3, :)
        call spectra(mpol, th_grid(self%nt), fld, jb, u, sv)
        sv%iota = 1.0_dp/q
        sv%curr_i = self%sig_f*mb
        sv%curr_g = self%sig_f*ct%f
    end subroutine surface_map

    pure subroutine route_integrand(self, route, ct, g)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: route
        type(contour_t), intent(in) :: ct
        real(dp), intent(out) :: g(:)
        real(dp) :: q, r, a(size(g))

        q = sum(ct%out(1, :))/real(self%nt, dp)
        select case (route)
        case (EFIT_SFL_COVARIANT)
            r = sum(ct%out(2, :))/real(self%nt, dp)/q
            g = (ct%out(2, :) - ct%out(1, :)*r)/(ct%f + r)
        case (EFIT_SFL_JACOBIAN_B2)
            g = q*ct%out(5, :)/(sum(ct%out(5, :))/real(self%nt, dp)) - ct%out(1, :)
        case (EFIT_SFL_JACOBIAN_CONST)
            g = q*ct%out(4, :)/(sum(ct%out(4, :))/real(self%nt, dp)) - ct%out(1, :)
        case (EFIT_SFL_STRAIGHT_CURRENT)
            a = -ct%out(4, :)*(ct%prr - ct%praw/ct%r + ct%pzz)/ct%r**2
            r = sum(a)/real(self%nt, dp)
            g = (q*a - ct%out(1, :)*r)/(r - q*self%sig_f*ct%fpsi)
        case default
            g = 0.0_dp
        end select
    end subroutine route_integrand

    !> NO_PRESSURE check of straight_current: the divisor per ds, |G' + iota I'|,
    !> against |K| = |G + iota I| = f + <bg>/q, and the file p'.
    subroutine straight_current_check(self, ct, st)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        type(efit_sfl_status_t), intent(inout) :: st
        real(dp) :: q, ma, den, ck, pp, pmax, dummy

        q = sum(ct%out(1, :))/real(self%nt, dp)
        ma = -sum(ct%out(4, :)*(ct%prr - ct%praw/ct%r + ct%pzz)/ct%r**2)/real(self%nt, dp)
        den = (ma/q - self%sig_f*ct%fpsi)*abs(self%psi_edge - self%psi_ax)*self%phi1/q
        ck = ct%f + sum(ct%out(2, :))/real(self%nt, dp)/q
        call prof_eval(self, self%cpp, ct%xf, pp, dummy)
        pmax = maxval(abs(self%g%pprime))
        if (pmax == 0.0_dp .or. abs(pp) <= 1.0e-8_dp*pmax .or. abs(den) < 1.0e-8_dp*ck) &
            call set_status(st, EFIT_SFL_ERR_NO_PRESSURE, 0, -1.0_dp, abs(den)/ck, 1.0e-8_dp, &
                'route straight_current (Hamada II) needs p'' /= 0, but p''(psi) = '//fmt(pp)// &
                ' (max |pprime| '//fmt(pmax)//') and |G'' + iota I''|/|K| = '//fmt(abs(den)/ck)// &
                '. eta = -(beta - <beta>)/(G'' + iota I'') is undefined without pressure '// &
                'gradient; use jacobian_const (Hamada I, route 3) or covariant (route 1)')
    end subroutine straight_current_check

    !> Quadrature b_m = w_m sum_j f_j J_j exp(i m theta*_j) over the points
    !> 1, 1 + stride, ... of the uniform grid th.
    subroutine spectra(mpol, th, fld, jb, u, sv, stride)
        integer, intent(in) :: mpol
        real(dp), intent(in) :: th(:), fld(:, :), jb(:), u(:)
        type(efit_sfl_surface_t), intent(inout) :: sv
        integer, intent(in), optional :: stride
        complex(dp), allocatable :: e(:), em(:)
        complex(dp) :: amp(4)
        real(dp) :: w
        integer :: j, m, st, n

        st = 1
        if (present(stride)) st = stride
        n = size(th)/st
        allocate (e(n), em(n))
        call alloc_surface(sv, mpol)
        e = exp(cmplx(0.0_dp, th(1::st) + u(1::st), dp))
        em = cmplx(jb(1::st), 0.0_dp, dp)
        do m = 0, mpol
            if (m > 0) em = em*e
            w = merge(1.0_dp, 2.0_dp, m == 0)/real(n, dp)
            do j = 1, 4
                amp(j) = w*sum(fld(1::st, j)*em)
            end do
            sv%rmnc(m) = real(amp(1))
            sv%rmns(m) = aimag(amp(1))
            sv%zmnc(m) = real(amp(2))
            sv%zmns(m) = aimag(amp(2))
            sv%vmnc(m) = real(amp(3))
            sv%vmns(m) = aimag(amp(3))
            sv%bmnc(m) = real(amp(4))
            sv%bmns(m) = aimag(amp(4))
        end do
    end subroutine spectra

    pure function th_grid(n) result(th)
        integer, intent(in) :: n
        real(dp) :: th(n)
        integer :: j
        do j = 1, n
            th(j) = twopi*real(j - 1, dp)/real(n, dp)
        end do
    end function th_grid

    subroutine alloc_surface(sv, mpol)
        type(efit_sfl_surface_t), intent(inout) :: sv
        integer, intent(in) :: mpol
        if (allocated(sv%rmnc)) then
            if (size(sv%rmnc) == mpol + 1) return
            deallocate (sv%rmnc, sv%rmns, sv%zmnc, sv%zmns, sv%vmnc, sv%vmns, sv%bmnc, sv%bmns)
        end if
        allocate (sv%rmnc(0:mpol), sv%rmns(0:mpol), sv%zmnc(0:mpol), sv%zmns(0:mpol), &
                  sv%vmnc(0:mpol), sv%vmns(0:mpol), sv%bmnc(0:mpol), sv%bmns(0:mpol))
        sv%rmnc = 0; sv%rmns = 0; sv%zmnc = 0; sv%zmns = 0
        sv%vmnc = 0; sv%vmns = 0; sv%bmnc = 0; sv%bmns = 0
    end subroutine alloc_surface

    !> Spectra on the surfaces s_list with the configured route. An error
    !> (no surface, or straight_current without p') is returned in status if
    !> present, raised otherwise.
    subroutine efit_sfl_transform(self, s_list, mpol, out, status)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: s_list(:)
        integer, intent(in) :: mpol
        type(efit_sfl_surface_t), intent(out) :: out(:)
        type(efit_sfl_status_t), intent(out), optional :: status
        type(efit_sfl_status_t) :: st
        type(contour_t) :: ct
        type(efit_sfl_surface_t) :: alt
        real(dp) :: g(self%nt), u(self%nt), jb(self%nt), b00
        integer :: is, other

        do is = 1, size(s_list)
            out(is)%s = s_list(is)
            out(is)%x = efit_sfl_x_of_s(self, s_list(is))
            call contour(self, out(is)%x, ct, st)
            if (st%code == EFIT_SFL_OK .and. self%cfg%route == EFIT_SFL_STRAIGHT_CURRENT) &
                call straight_current_check(self, ct, st)
            if (st%code /= EFIT_SFL_OK) then
                st%s = s_list(is)
                st%isurf = is
                if (present(status)) then
                    status = st
                    return
                end if
                call efit_sfl_raise(st)
            end if
            out(is)%psi = ct%psit
            call surface_map(self, self%cfg%route, mpol, ct, out(is), g, u, jb)
            call flux_functions(self, ct, out(is))
            other = merge(EFIT_SFL_JACOBIAN_B2, EFIT_SFL_COVARIANT, &
                          self%cfg%route == EFIT_SFL_COVARIANT)
            ! d_H: Hamada II against Hamada I measures the force residual
            if (self%cfg%route == EFIT_SFL_STRAIGHT_CURRENT) other = EFIT_SFL_JACOBIAN_CONST
            if (self%cfg%route == EFIT_SFL_COVARIANT .or. &
                self%cfg%route == EFIT_SFL_JACOBIAN_B2 .or. &
                self%cfg%route == EFIT_SFL_STRAIGHT_CURRENT) then
                call surface_map(self, other, mpol, ct, alt, g, u, jb)
                b00 = abs(out(is)%bmnc(0))
                out(is)%d_family = max(maxval(abs(out(is)%bmnc - alt%bmnc)), &
                                       maxval(abs(out(is)%bmns - alt%bmns)))/b00
            end if
        end do
        call resolution_check(self, s_list, mpol, out)
        if (present(status)) status = st
    end subroutine efit_sfl_transform

    !> Surface x over the symmetry-flux (PEST) angle theta, the angle of the
    !> legacy magdata grids: spectra of R, Z, v = -nu/(2 pi) with nu of the
    !> configured route, and |B|, by quadrature over theta_g; flux functions
    !> as in efit_sfl_transform. The contour is returned for further integrals.
    subroutine efit_sfl_symflux(self, x, mpol, sv, ct, st)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: mpol
        type(efit_sfl_surface_t), intent(inout) :: sv
        type(contour_t), intent(out) :: ct
        type(efit_sfl_status_t), intent(inout) :: st
        real(dp) :: g(self%nt), u(self%nt), jb(self%nt), fld(self%nt, 4)

        call contour(self, x, ct, st)
        if (st%code /= EFIT_SFL_OK) return
        call surface_map(self, EFIT_SFL_PEST, mpol, ct, sv, g, u, jb)
        call route_integrand(self, self%cfg%route, ct, g)
        fld(:, 1) = ct%r
        fld(:, 2) = ct%z
        fld(:, 3) = -matmul(self%aop, g)/twopi
        fld(:, 4) = ct%out(3, :)
        call spectra(mpol, th_grid(self%nt), fld, jb, u, sv)
        sv%x = x
        sv%psi = ct%psit
    end subroutine efit_sfl_symflux

    !> Recompute the outermost surface with 2 n_theta and record the change of
    !> b_mn relative to b_00 in eps_theta of that surface (warning above 1e-9).
    subroutine resolution_check(self, s_list, mpol, out)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: s_list(:)
        integer, intent(in) :: mpol
        type(efit_sfl_surface_t), intent(inout) :: out(:)
        type(efit_sfl_t), allocatable :: fine
        type(efit_sfl_surface_t) :: sv
        type(efit_sfl_status_t) :: st
        type(contour_t) :: ct
        integer :: is
        real(dp), allocatable :: g(:), u(:), jb(:)

        is = maxloc(s_list, dim=1)
        allocate (fine)
        fine = self
        fine%nt = 2*self%nt
        deallocate (fine%aop)
        call antiderivative_matrix(fine%nt, fine%aop)
        allocate (g(fine%nt), u(fine%nt), jb(fine%nt))
        call contour(fine, out(is)%x, ct, st)
        if (st%code /= EFIT_SFL_OK) return
        call surface_map(fine, self%cfg%route, mpol, ct, sv, g, u, jb)
        out(is)%eps_theta = max(maxval(abs(sv%bmnc - out(is)%bmnc)), &
                                maxval(abs(sv%bmns - out(is)%bmns)), &
                                maxval(abs(sv%rmnc - out(is)%rmnc))*abs(out(is)%bmnc(0)/out(is)%rmnc(0)), &
                                maxval(abs(sv%rmns - out(is)%rmns))*abs(out(is)%bmnc(0)/out(is)%rmnc(0)))/ &
                            abs(out(is)%bmnc(0))
        if (out(is)%eps_theta > 1.0e-9_dp) write (error_unit, '(a, es9.2, a, f7.4, a)') &
            'efit_sfl: warning: doubling n_theta changes the spectra by ', out(is)%eps_theta, &
            ' (relative) at s = ', out(is)%s, '; raise n_theta'
    end subroutine resolution_check

    !> sqrt g(0,0) (Jacobian in psi_tor, mean <jc>/q, sign of the legacy
    !> header), file and force-balance p'(s) and the averaged force balance
    !> eps_fb = |4 pi J p' + G' + iota I'|/(|4 pi J p'| + |G'| + |iota I'|), d/ds,
    !> with p' and G' = F F'/F from the pprime and ffprim columns.
    subroutine flux_functions(self, ct, sv)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        type(efit_sfl_surface_t), intent(inout) :: sv
        real(dp) :: q, dq, di, dummy, pp, ffp, dpsids, jac, dg_ds, di_ds, a

        call qfun(self, ct%x, q, dq)
        call qval(self, self%ciq, ct%x, dummy, di)
        di = di/dxdu(umap(ct%x))
        jac = sum(ct%out(4, :))/real(self%nt, dp)/(sum(ct%out(1, :))/real(self%nt, dp))
        sv%sqrtg00 = self%sig_f*jac
        dpsids = (self%psi_edge - self%psi_ax)*self%phi1/q
        call prof_eval(self, self%cpp, ct%xf, pp, dummy)
        sv%pp_file = pp*dpsids
        call prof_eval(self, self%cffp, ct%xf, ffp, dummy)
        dg_ds = ffp/(self%sig_f*ct%f)*dpsids
        di_ds = di*self%phi1/q
        sv%pp_fb = -(dg_ds + sv%iota*di_ds)/(2.0_dp*twopi*sv%sqrtg00)
        a = 2.0_dp*twopi*sv%sqrtg00*sv%pp_file
        sv%eps_fb = abs(a + dg_ds + sv%iota*di_ds)/max(abs(a) + abs(dg_ds) + &
                                                       abs(sv%iota*di_ds), tiny(1.0_dp))
    end subroutine flux_functions

    ! -------------------------------------------------------------- checks

    !> Least-squares factor c in Delta* psi = c x gsrhs over the window nodes.
    real(dp) function gs_fit_factor(self)
        type(efit_sfl_t), intent(in) :: self
        real(dp) :: num, den, d, rhs, w
        integer :: j, k

        num = 0.0_dp
        den = 0.0_dp
        do j = 1, self%nq
            if (.not. in_window(self, j)) cycle
            do k = 1, self%nt
                call gs_point(self, self%qn(j), k, d, rhs, w)
                num = num + w*d*rhs
                den = den + w*rhs**2
            end do
        end do
        gs_fit_factor = num/den
    end function gs_fit_factor

    logical function in_window(self, j)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: j
        real(dp) :: rho
        rho = sqrt(qint(self, self%cq, self%xq(j))/self%phi1)
        in_window = rho >= self%cfg%rho_window(1) .and. rho <= self%cfg%rho_window(2)
    end function in_window

    subroutine gs_point(self, ct, k, dstar, rhs, w, primaries)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        integer, intent(in) :: k
        real(dp), intent(out) :: dstar, rhs, w
        logical, intent(in), optional :: primaries
        real(dp) :: pp, ffp, f, y(6), jac(6, 9)

        call gs_profiles(self, ct%xf, present(primaries), pp, ffp, f)
        call field_kernel(ct%r(k), ct%praw(k), ct%pzraw(k), ct%prr(k), ct%prz(k), &
                          ct%pzz(k), f, 2.0_dp*twopi*pp, ffp, y, jac)
        dstar = y(5)
        rhs = y(6)
        w = abs(ct%rho(k)/ct%prho(k))   ! area element per dpsi dtheta_g
    end subroutine gs_point

    !> p'(psi) and F F'(psi) at x_file: the pprime and ffprim columns (the
    !> source terms the equilibrium solver used, default) or the derivatives of
    !> the integrated pres and fpol columns (primaries = .true.). fpol is printed
    !> with 9 digits, so F F' from fpol carries rounding noise ~ n_w 1e-9 F/dF
    !> (FreeGS LSN: 2e-5 at 129, 5e-5 at 257 points, boozer-bench efit-sfl report).
    pure subroutine gs_profiles(self, x, primaries, pp, ffp, f)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: x
        logical, intent(in) :: primaries
        real(dp), intent(out) :: pp, ffp, f
        real(dp) :: dummy, dpdx, dfdx

        call prof_eval(self, self%cf, x, f, dfdx)
        if (.not. primaries) then
            call prof_eval(self, self%cpp, x, pp, dummy)
            call prof_eval(self, self%cffp, x, ffp, dummy)
        else
            call prof_eval(self, self%cpres, x, dummy, dpdx)
            pp = dpdx/(self%g%sibry - self%g%simag)
            ffp = f*dfdx/(self%g%sibry - self%g%simag)
        end if
    end subroutine gs_profiles

    pure subroutine field_kernel(r, pr, pz, prr, prz, pzz, f, mpp, ffp, y, jac)
        real(dp), intent(in) :: r, pr, pz, prr, prz, pzz, f, mpp, ffp
        real(dp), intent(out) :: y(6), jac(6, 9)
        call efit_sfl_field(r, pr, pz, prr, prz, pzz, f, mpp, ffp, &
            y(1), jac(1, 1), jac(1, 2), jac(1, 3), jac(1, 4), jac(1, 5), jac(1, 6), jac(1, 7), &
            jac(1, 8), jac(1, 9), &
            y(2), jac(2, 1), jac(2, 2), jac(2, 3), jac(2, 4), jac(2, 5), jac(2, 6), jac(2, 7), &
            jac(2, 8), jac(2, 9), &
            y(3), jac(3, 1), jac(3, 2), jac(3, 3), jac(3, 4), jac(3, 5), jac(3, 6), jac(3, 7), &
            jac(3, 8), jac(3, 9), &
            y(4), jac(4, 1), jac(4, 2), jac(4, 3), jac(4, 4), jac(4, 5), jac(4, 6), jac(4, 7), &
            jac(4, 8), jac(4, 9), &
            y(5), jac(5, 1), jac(5, 2), jac(5, 3), jac(5, 4), jac(5, 5), jac(5, 6), jac(5, 7), &
            jac(5, 8), jac(5, 9), &
            y(6), jac(6, 1), jac(6, 2), jac(6, 3), jac(6, 4), jac(6, 5), jac(6, 6), jac(6, 7), &
            jac(6, 8), jac(6, 9))
    end subroutine field_kernel

    subroutine run_checks(self, st)
        type(efit_sfl_t), intent(inout) :: self
        type(efit_sfl_status_t), intent(inout) :: st
        type(efit_sfl_diag_t) :: d
        type(efit_sfl_surface_t) :: sv
        real(dp) :: num, den, nj, dj, dstar, rhs, w, qp, dummy, x, pres_x, pp, f, fx, ffp, &
                    numc, denc
        real(dp), allocatable :: fb(:)
        integer :: j, k, nq, jworst

        nq = self%nq
        d%eps_axis = self%diag%eps_axis
        d%gs_factor = self%diag%gs_factor
        d%exp_bpol = self%diag%exp_bpol
        allocate (d%x(nq), d%s(nq), d%q(nq), d%qpsi(nq), d%eps_q(nq), d%eps_gs(nq), &
                  d%eps_fb(nq), d%curr_i(nq), d%curr_g(nq), d%pp_file(nq), d%pp_fb(nq), &
                  d%in_window(nq))
        num = 0.0_dp
        den = 0.0_dp
        numc = 0.0_dp
        denc = 0.0_dp
        d%eps_pp = 0.0_dp
        d%eps_ff = 0.0_dp
        do j = 1, nq
            d%x(j) = self%xq(j)
            d%s(j) = qint(self, self%cq, self%xq(j))/self%phi1
            d%in_window(j) = in_window(self, j)
            d%q(j) = self%qq(j)
            call prof_eval(self, self%cqpsi, self%qn(j)%xf, qp, dummy)
            d%qpsi(j) = abs(qp)
            d%eps_q(j) = abs(self%qq(j)/abs(qp) - 1.0_dp)
            nj = 0.0_dp
            dj = 0.0_dp
            do k = 1, self%nt
                call gs_point(self, self%qn(j), k, dstar, rhs, w)
                nj = nj + w*(dstar - rhs)**2
                dj = dj + w*rhs**2
                if (in_window(self, j)) then
                    call gs_point(self, self%qn(j), k, dstar, rhs, w, primaries=.true.)
                    numc = numc + w*(dstar - rhs)**2
                    denc = denc + w*rhs**2
                end if
            end do
            d%eps_gs(j) = sqrt(nj/dj)
            if (d%in_window(j)) then
                num = num + nj
                den = den + dj
            end if
            sv%x = self%xq(j)
            sv%iota = 1.0_dp/self%qq(j)
            call flux_functions(self, self%qn(j), sv)
            d%curr_i(j) = self%iq(j)
            d%curr_g(j) = self%sig_f*self%qn(j)%f
            d%pp_file(j) = sv%pp_file
            d%pp_fb(j) = sv%pp_fb
            d%eps_fb(j) = sv%eps_fb
            ! Column consistency: pprime vs d pres/dpsi, ffprim vs F dF/dpsi.
            if (d%in_window(j)) then
                x = self%qn(j)%xf
                call prof_eval(self, self%cpres, x, dummy, pres_x)
                call prof_eval(self, self%cpp, x, pp, dummy)
                call prof_eval(self, self%cf, x, f, fx)
                call prof_eval(self, self%cffp, x, ffp, dummy)
                d%eps_pp = max(d%eps_pp, abs(pres_x/(self%g%sibry - self%g%simag) - pp)/ &
                               max(maxval(abs(self%g%pprime)), tiny(1.0_dp)))
                d%eps_ff = max(d%eps_ff, abs(f*fx/(self%g%sibry - self%g%simag) - ffp)/ &
                               max(maxval(abs(self%g%ffprim)), tiny(1.0_dp)))
            end if
        end do
        d%eps_gs_window = sqrt(num/max(den, tiny(1.0_dp)))
        d%eps_gs_pres_fpol = sqrt(numc/max(denc, tiny(1.0_dp)))
        d%eps_q_window = maxval(d%eps_q, mask=d%in_window)
        allocate (fb(count(d%in_window)))
        fb = pack(d%eps_fb, d%in_window)
        d%eps_fb_median = median(fb)
        fb = pack(d%eps_gs, d%in_window)
        d%eps_gs_median = median(fb)
        self%diag = d

        if (d%eps_gs_median > self%cfg%tol_gs .or. self%diag%gs_factor < 0.0_dp) then
            jworst = maxloc(d%eps_gs, mask=d%in_window, dim=1)
            call set_status(st, EFIT_SFL_ERR_GS, jworst, d%s(jworst), d%eps_gs_median, &
                self%cfg%tol_gs, 'Grad-Shafranov residual of psi with the file pprime and '// &
                'ffprim, median over the window eps_gs = '//fmt(d%eps_gs_median)//' > tol '// &
                fmt(self%cfg%tol_gs)//' (worst surface eps_gs = '//fmt(d%eps_gs(jworst))// &
                ', fitted factor '//fmt(self%diag%gs_factor)//'). Regenerate the G-EQDSK '// &
                'with a converged solver, a finer grid or profile columns on the psi grid')
        else if (d%eps_fb_median > self%cfg%tol_fb) then
            jworst = maxloc(d%eps_fb, mask=d%in_window, dim=1)
            call set_status(st, EFIT_SFL_ERR_FORCE_BALANCE, jworst, d%s(jworst), &
                d%eps_fb_median, self%cfg%tol_fb, 'averaged force balance 4 pi J p'' + G'' '// &
                '+ iota I'' = 0 with the file pprime, median eps_fb = '// &
                fmt(d%eps_fb_median)//' > tol '//fmt(self%cfg%tol_fb)//' (worst surface '// &
                fmt(d%eps_fb(jworst))//'). Regenerate the G-EQDSK with a converged solver')
        else if (d%eps_q_window > self%cfg%tol_q) then
            jworst = maxloc(d%eps_q, mask=d%in_window, dim=1)
            call set_status(st, EFIT_SFL_ERR_Q, jworst, d%s(jworst), d%eps_q_window, &
                self%cfg%tol_q, 'qpsi column disagrees with the contour q of psi and F, '// &
                'q_contour = '//fmt(d%q(jworst))//', qpsi = '//fmt(d%qpsi(jworst))// &
                '. The output uses the contour q; check how qpsi was produced')
        end if
    end subroutine run_checks

    real(dp) function median(a)
        real(dp), intent(in) :: a(:)
        real(dp) :: b(size(a)), t
        integer :: i, j
        if (size(a) == 0) then
            median = 0.0_dp
            return
        end if
        b = a
        do i = 2, size(b)
            t = b(i)
            j = i - 1
            do while (j >= 1)
                if (b(j) <= t) exit
                b(j + 1) = b(j)
                j = j - 1
            end do
            b(j + 1) = t
        end do
        median = 0.5_dp*(b((size(b) + 1)/2) + b(size(b)/2 + 1))
    end function median

    function fmt(v) result(s)
        real(dp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write (buf, '(es10.2)') v
        buf = adjustl(buf)
        allocate (character(len=len_trim(buf)) :: s)
        s = trim(buf)
    end function fmt

    subroutine set_status(st, code, isurf, s, value, tol, cause)
        type(efit_sfl_status_t), intent(inout) :: st
        integer, intent(in) :: code, isurf
        real(dp), intent(in) :: s, value, tol
        character(*), intent(in) :: cause
        character(len=64) :: where

        st%code = code
        st%isurf = isurf
        st%s = s
        st%value = value
        st%tol = tol
        if (s >= 0.0_dp) then
            write (where, '("s = ", f7.4, " (rho ", f6.3, ")")') s, sqrt(max(s, 0.0_dp))
        else
            where = 'file level'
        end if
        if (code == EFIT_SFL_ERR_INPUT .or. code == EFIT_SFL_ERR_NO_PRESSURE) then
            call put('efit_sfl: invalid input at '//trim(where)//': '//cause//'.')
        else
            call put('efit_sfl: inconsistent equilibrium at '//trim(where)//': '//cause// &
                     '. Override: allow_poor_equilibrium (logged and written to the output).')
        end if
    contains
        subroutine put(text)
            character(*), intent(in) :: text
            if (allocated(st%message)) deallocate (st%message)
            allocate (character(len=len(text)) :: st%message)
            st%message = text
        end subroutine put
    end subroutine set_status

    subroutine finish_status(self, st, status)
        type(efit_sfl_t), intent(inout) :: self
        type(efit_sfl_status_t), intent(inout) :: st
        type(efit_sfl_status_t), intent(out), optional :: status

        if (st%code /= EFIT_SFL_OK .and. st%code /= EFIT_SFL_ERR_INPUT .and. &
            self%cfg%allow_poor_equilibrium) then
            st%overridden = .true.
            self%diag%override = .true.
            write (error_unit, '(a)') 'efit_sfl: OVERRIDE ('// &
                trim(self%cfg%override_reason)//') '//st%message
        end if
        self%diag%status = st
        if (present(status)) then
            status = st
        else if (st%code /= EFIT_SFL_OK .and. .not. st%overridden) then
            call efit_sfl_raise(st)
        end if
    end subroutine finish_status

    subroutine efit_sfl_raise(st)
        type(efit_sfl_status_t), intent(in) :: st
        if (st%code == EFIT_SFL_OK .or. st%overridden) return
        if (allocated(st%message)) then
            write (error_unit, '(a)') st%message
        else
            write (error_unit, '(a, i0)') 'efit_sfl: error code ', st%code
        end if
        error stop 1
    end subroutine efit_sfl_raise

    !> Per-node diagnostics table, comment lines start with '#'.
    subroutine efit_sfl_diag_write(self, unit)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: unit
        integer :: j
        associate (d => self%diag)
            write (unit, '(a)') '# efit_sfl diagnostics (Gaussian units, p'' = dp/ds)'
            write (unit, '(a, es12.4, a, es12.4, a, es12.4, a, es12.4, a, es12.4)') &
                '# eps_gs ', d%eps_gs_window, ' eps_q ', d%eps_q_window, ' eps_fb_median ', &
                d%eps_fb_median, ' eps_axis ', d%eps_axis, ' eps_gs_median ', d%eps_gs_median
            write (unit, '(a, es12.4, a, es12.4, a, es12.4, a, i2, a, es12.4)') &
                '# eps_pprime_pres ', d%eps_pp, ' eps_ffprim_fpol ', d%eps_ff, ' gs_factor ', &
                d%gs_factor, ' exp_bpol ', d%exp_bpol, ' eps_gs_pres_fpol ', d%eps_gs_pres_fpol
            write (unit, '(a, 2f6.3, a, es10.2, a, es10.2, a, l2, a, i2)') '# rho_window', &
                self%cfg%rho_window, ' tol_gs', self%cfg%tol_gs, ' tol_q', self%cfg%tol_q, &
                ' override', d%override, ' status', d%status%code
            write (unit, '(a)') '#  x  s  q_contour  qpsi  eps_q  eps_gs  eps_fb  I  G  '// &
                'pp_file  pp_fb  in_window'
            do j = 1, size(d%x)
                write (unit, '(11es24.15, l3)') d%x(j), d%s(j), d%q(j), d%qpsi(j), &
                    d%eps_q(j), d%eps_gs(j), d%eps_fb(j), d%curr_i(j), d%curr_g(j), &
                    d%pp_file(j), d%pp_fb(j), d%in_window(j)
            end do
        end associate
    end subroutine efit_sfl_diag_write

    !> Grad-Shafranov residual (Delta* psi - gsrhs)/max|gsrhs| on the file grid;
    !> mask = 1 inside psi_edge (rays from the axis), else 0.
    subroutine efit_sfl_gs_map(self, res, mask)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(out) :: res(:, :), mask(:, :)
        real(dp) :: v(6), x, pp, ffp, f, y(6), jac(6, 9), scale
        integer :: i, j

        scale = 0.0_dp
        do i = 1, self%nw
            do j = 1, self%nh
                call psi_eval(self, self%cpsi, self%g%R_eqd(i), self%g%Z_eqd(j), v)
                x = (v(1) - self%psi_ax)/(self%psi_edge - self%psi_ax)
                mask(i, j) = merge(1.0_dp, 0.0_dp, x >= 0.0_dp .and. x <= 1.0_dp .and. &
                                   enclosed(self, self%g%R_eqd(i), self%g%Z_eqd(j)))
                call gs_profiles(self, xfile(self, v(1)), .false., pp, ffp, f)
                call field_kernel(self%g%R_eqd(i), v(2), v(3), v(4), v(5), v(6), f, &
                                  2.0_dp*twopi*pp, ffp, y, jac)
                res(i, j) = y(5) - y(6)
                if (mask(i, j) > 0.0_dp) scale = max(scale, abs(y(6)))
            end do
        end do
        res = res/scale
    end subroutine efit_sfl_gs_map

    !> True if (R, Z) lies inside the psi_edge contour: psi monotonic on the
    !> segment from the axis.
    logical function enclosed(self, r, z)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: r, z
        real(dp) :: v(6), prev, cur
        integer :: k
        integer, parameter :: nk = 32
        enclosed = .true.
        prev = 0.0_dp
        do k = 1, nk
            call psi_eval(self, self%cpsi, self%r0 + (r - self%r0)*real(k, dp)/nk, &
                          self%z0 + (z - self%z0)*real(k, dp)/nk, v)
            cur = self%sig_psi*(v(1) - self%psi_ax)
            if (cur < prev .or. cur > abs(self%psi_edge - self%psi_ax)) then
                enclosed = .false.
                return
            end if
            prev = cur
        end do
    end function enclosed

    ! ------------------------------------------------------- derivatives

    !> Tangent of efit_sfl_transform for a perturbation of psirz and fpol.
    subroutine efit_sfl_jvp(self, s_list, mpol, tan, out_d)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: s_list(:)
        integer, intent(in) :: mpol
        type(efit_sfl_tangent_t), intent(in) :: tan
        type(efit_sfl_surface_t), intent(out) :: out_d(:)
        real(dp), allocatable :: dc(:, :, :, :), dcf(:, :)
        real(dp) :: dpsi_ax, dr0(2), dq(self%nq), wint(self%nq), wval(self%nq), &
                    wint1(self%nq), dphi1, dx, dpsit, q, dqx, s
        type(efit_sfl_status_t) :: st
        type(contour_t) :: ct
        integer :: j, is

        call psi_coefficients(self, tan%psirz, dc)
        allocate (dcf(0:5, self%nw))
        dcf = 0.0_dp
        dcf(0, :) = tan%fpol
        call spl_reg(5, self%nw, self%hx, dcf)
        call axis_jvp(self, dc, dpsi_ax, dr0)
        do j = 1, self%nq
            dpsit = dpsi_ax*(1.0_dp - self%xq(j))
            call contour_mean_jvp(self, self%qn(j), dc, dcf, dpsit, dr0, dq(j))
            dq(j) = dq(j)*dxdu((real(j, dp) - 0.5_dp)/real(self%nq, dp))
        end do
        call qweights(self, 1.0_dp, wint1, wval)
        dphi1 = dot_product(wint1, dq)
        do is = 1, size(s_list)
            s = s_list(is)
            out_d(is)%x = efit_sfl_x_of_s(self, s)
            call qweights(self, out_d(is)%x, wint, wval)
            call qfun(self, out_d(is)%x, q, dqx)
            dx = -(dot_product(wint, dq) - s*dphi1)/q
            call contour(self, out_d(is)%x, ct, st)
            dpsit = dpsi_ax*(1.0_dp - ct%x) + (self%psi_edge - self%psi_ax)*dx
            call surface_jvp(self, mpol, ct, dc, dcf, dpsit, dr0, out_d(is))
            out_d(is)%s = 0.0_dp
            out_d(is)%x = dx
            out_d(is)%psi = dpsit
        end do
    end subroutine efit_sfl_jvp

    subroutine axis_jvp(self, dc, dpsi_ax, dr0)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: dc(0:, 0:, :, :)
        real(dp), intent(out) :: dpsi_ax, dr0(2)
        real(dp) :: v(6), det

        call psi_eval(self, dc, self%r0, self%z0, v)
        dpsi_ax = v(1)
        det = self%hess(1, 1)*self%hess(2, 2) - self%hess(1, 2)**2
        dr0(1) = -(self%hess(2, 2)*v(2) - self%hess(1, 2)*v(3))/det
        dr0(2) = -(-self%hess(1, 2)*v(2) + self%hess(1, 1)*v(3))/det
    end subroutine axis_jvp

    !> Pointwise tangents (rho, R, Z, kernel outputs) of a contour.
    subroutine contour_jvp(self, ct, dc, dcf, dpsit, dr0, dout, drr, dzz, df)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: dc(0:, 0:, :, :), dcf(0:, :), dpsit, dr0(2)
        real(dp), intent(out) :: dout(:, :), drr(:), dzz(:), df
        real(dp) :: v(6), drho, dpr, dpz, din(5), dfd, dummy
        integer :: k

        call prof_eval(self, dcf, ct%xf, dfd, dummy)
        df = self%sig_f*(dfd + ct%fpsi*dpsit)
        do k = 1, self%nt
            call psi_eval(self, dc, ct%r(k), ct%z(k), v)
            drho = (dpsit - v(1) - ct%praw(k)*dr0(1) - ct%pzraw(k)*dr0(2))/ct%prho(k)
            drr(k) = dr0(1) + ct%ca(k)*drho
            dzz(k) = dr0(2) + ct%sa(k)*drho
            dpr = self%sig_psi*(v(2) + ct%prr(k)*drr(k) + ct%prz(k)*dzz(k))
            dpz = self%sig_psi*(v(3) + ct%prz(k)*drr(k) + ct%pzz(k)*dzz(k))
            din = [drho, drr(k), dpr, dpz, df]
            dout(:, k) = matmul(ct%jac(:, :, k), din)
        end do
    end subroutine contour_jvp

    subroutine contour_mean_jvp(self, ct, dc, dcf, dpsit, dr0, dq)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: dc(0:, 0:, :, :), dcf(0:, :), dpsit, dr0(2)
        real(dp), intent(out) :: dq
        real(dp) :: dout(5, self%nt), drr(self%nt), dzz(self%nt), df
        call contour_jvp(self, ct, dc, dcf, dpsit, dr0, dout, drr, dzz, df)
        dq = sum(dout(1, :))/real(self%nt, dp)
    end subroutine contour_mean_jvp

    subroutine surface_jvp(self, mpol, ct, dc, dcf, dpsit, dr0, sd)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: mpol
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: dc(0:, 0:, :, :), dcf(0:, :), dpsit, dr0(2)
        type(efit_sfl_surface_t), intent(inout) :: sd
        real(dp) :: dout(5, self%nt), drr(self%nt), dzz(self%nt), df
        real(dp) :: g(self%nt), dg(self%nt), a(self%nt), da(self%nt), u(self%nt), &
                    du(self%nt), jb(self%nt), djb(self%nt), nu(self%nt), dnu(self%nt)
        real(dp) :: q, dqq, mb, dmb, fld(self%nt, 4), dfld(self%nt, 4)
        complex(dp) :: e(self%nt), em(self%nt), amp(4)
        real(dp) :: w
        integer :: j, m

        call contour_jvp(self, ct, dc, dcf, dpsit, dr0, dout, drr, dzz, df)
        associate (h => ct%out(1, :), dh => dout(1, :))
            q = sum(h)/real(self%nt, dp)
            dqq = sum(dh)/real(self%nt, dp)
            mb = sum(ct%out(2, :))/real(self%nt, dp)
            dmb = sum(dout(2, :))/real(self%nt, dp)
            call route_integrand(self, self%cfg%route, ct, g)
            if (self%cfg%route == EFIT_SFL_STRAIGHT_CURRENT) then
                call straight_current_jvp(self, ct, dc, dcf, dpsit, dout, drr, dzz, dg)
            else
                call route_jvp(self, self%cfg%route, ct, dout, df, dg)
            end if
            a = h + g
            da = dh + dg
            nu = matmul(self%aop, g)
            dnu = matmul(self%aop, dg)
            u = matmul(self%aop, a)/q
            du = matmul(self%aop, da)/q - u*dqq/q
            jb = a/q
            djb = da/q - jb*dqq/q
        end associate
        fld(:, 1) = ct%r
        fld(:, 2) = ct%z
        fld(:, 3) = -nu/twopi
        fld(:, 4) = ct%out(3, :)
        dfld(:, 1) = drr
        dfld(:, 2) = dzz
        dfld(:, 3) = -dnu/twopi
        dfld(:, 4) = dout(3, :)
        call alloc_surface(sd, mpol)
        do j = 1, self%nt
            e(j) = exp(cmplx(0.0_dp, twopi*real(j - 1, dp)/real(self%nt, dp) + u(j), dp))
        end do
        em = cmplx(1.0_dp, 0.0_dp, dp)
        do m = 0, mpol
            if (m > 0) em = em*e
            w = merge(1.0_dp, 2.0_dp, m == 0)/real(self%nt, dp)
            do j = 1, 4
                amp(j) = w*sum((dfld(:, j)*jb + fld(:, j)*djb + &
                                fld(:, j)*jb*cmplx(0.0_dp, real(m, dp)*du, dp))*em)
            end do
            sd%rmnc(m) = real(amp(1)); sd%rmns(m) = aimag(amp(1))
            sd%zmnc(m) = real(amp(2)); sd%zmns(m) = aimag(amp(2))
            sd%vmnc(m) = real(amp(3)); sd%vmns(m) = aimag(amp(3))
            sd%bmnc(m) = real(amp(4)); sd%bmns(m) = aimag(amp(4))
        end do
        sd%iota = -dqq/q**2
        sd%curr_i = self%sig_f*dmb
        sd%curr_g = self%sig_f*df
    end subroutine surface_jvp

    pure subroutine route_jvp(self, route, ct, dout, df, dg)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: route
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: dout(:, :), df
        real(dp), intent(out) :: dg(:)
        real(dp) :: q, dq, mb, dmb, r, dr, den, dden, mj, dmj
        real(dp) :: num(self%nt), dnum(self%nt)
        integer :: k

        q = sum(ct%out(1, :))/real(self%nt, dp)
        dq = sum(dout(1, :))/real(self%nt, dp)
        select case (route)
        case (EFIT_SFL_COVARIANT)
            mb = sum(ct%out(2, :))/real(self%nt, dp)
            dmb = sum(dout(2, :))/real(self%nt, dp)
            r = mb/q
            dr = dmb/q - mb*dq/q**2
            num = ct%out(2, :) - ct%out(1, :)*r
            dnum = dout(2, :) - dout(1, :)*r - ct%out(1, :)*dr
            den = ct%f + r
            dden = df + dr
            dg = dnum/den - num*dden/den**2
        case (EFIT_SFL_JACOBIAN_B2, EFIT_SFL_JACOBIAN_CONST)
            k = merge(5, 4, route == EFIT_SFL_JACOBIAN_B2)
            mj = sum(ct%out(k, :))/real(self%nt, dp)
            dmj = sum(dout(k, :))/real(self%nt, dp)
            dg = dq*ct%out(k, :)/mj + q*dout(k, :)/mj - q*ct%out(k, :)*dmj/mj**2 - dout(1, :)
        case default
            dg = 0.0_dp
        end select
    end subroutine route_jvp

    !> Adjoint of efit_sfl_jvp: tan_b = J^T out_b.
    subroutine efit_sfl_vjp(self, s_list, mpol, out_b, tan_b)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: s_list(:)
        integer, intent(in) :: mpol
        type(efit_sfl_surface_t), intent(in) :: out_b(:)
        type(efit_sfl_tangent_t), intent(out) :: tan_b
        real(dp), allocatable :: cb(:, :, :, :), cfb(:)
        real(dp) :: psi_ax_b, r0_b(2), qb(self%nq), wint(self%nq), wval(self%nq), &
                    wint1(self%nq), xb, psitb, q, dqx, s, outb(5, self%nt)
        type(efit_sfl_status_t) :: st
        type(contour_t) :: ct
        integer :: j, is

        allocate (cb(self%nw, self%nh, 0:5, 0:5), cfb(self%nw))
        cb = 0.0_dp
        cfb = 0.0_dp
        psi_ax_b = 0.0_dp
        r0_b = 0.0_dp
        qb = 0.0_dp
        call qweights(self, 1.0_dp, wint1, wval)
        do is = 1, size(s_list)
            s = s_list(is)
            call contour(self, efit_sfl_x_of_s(self, s), ct, st)
            psitb = out_b(is)%psi
            call surface_vjp(self, mpol, ct, out_b(is), cb, cfb, psitb, r0_b)
            psi_ax_b = psi_ax_b + psitb*(1.0_dp - ct%x)
            xb = out_b(is)%x + psitb*(self%psi_edge - self%psi_ax)
            call qweights(self, ct%x, wint, wval)
            call qfun(self, ct%x, q, dqx)
            qb = qb - xb*(wint - s*wint1)/q
        end do
        do j = 1, self%nq
            outb = 0.0_dp
            outb(1, :) = qb(j)*dxdu((real(j, dp) - 0.5_dp)/real(self%nq, dp))/real(self%nt, dp)
            psitb = 0.0_dp
            call contour_vjp(self, self%qn(j), outb, [(0.0_dp, is=1, self%nt)], &
                             [(0.0_dp, is=1, self%nt)], 0.0_dp, cb, cfb, psitb, r0_b)
            psi_ax_b = psi_ax_b + psitb*(1.0_dp - self%xq(j))
        end do
        call axis_vjp(self, psi_ax_b, r0_b, cb)
        allocate (tan_b%psirz(self%nw, self%nh), tan_b%fpol(self%nw))
        call psi_coefficients_t(self, cb, tan_b%psirz)
        tan_b%fpol = cfb
    end subroutine efit_sfl_vjp

    subroutine axis_vjp(self, psi_ax_b, r0_b, cb)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: psi_ax_b, r0_b(2)
        real(dp), intent(inout) :: cb(:, :, 0:, 0:)
        real(dp) :: det, gb(2)

        det = self%hess(1, 1)*self%hess(2, 2) - self%hess(1, 2)**2
        gb(1) = -(self%hess(2, 2)*r0_b(1) - self%hess(1, 2)*r0_b(2))/det
        gb(2) = -(-self%hess(1, 2)*r0_b(1) + self%hess(1, 1)*r0_b(2))/det
        call psi_eval_t(self, self%r0, self%z0, [psi_ax_b, gb(1), gb(2)], cb)
    end subroutine axis_vjp

    !> data_b(a, b) = sum_{k,l} wr_k^T cb_kl wz_l.
    subroutine psi_coefficients_t(self, cb, psib)
        type(efit_sfl_t), intent(in) :: self
        real(dp), intent(in) :: cb(:, :, 0:, 0:)
        real(dp), intent(out) :: psib(:, :)
        integer :: k, l

        psib = 0.0_dp
        do l = 0, 5
            do k = 0, 5
                psib = psib + matmul(transpose(self%wr(:, :, k)), &
                                     matmul(cb(:, :, k, l), self%wz(:, :, l)))
            end do
        end do
    end subroutine psi_coefficients_t

    !> Transpose of contour_jvp: bars on kernel outputs, R, Z and f to the
    !> spline data, the surface psi_t, the axis position and fpol.
    subroutine contour_vjp(self, ct, outb, rb, zb, fb, cb, cfb, psitb, r0_b)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: outb(:, :), rb(:), zb(:), fb
        real(dp), intent(inout) :: cb(:, :, 0:, 0:), cfb(:), psitb, r0_b(2)
        real(dp) :: inb(5), rt, zt, gr, rhob, ftot, t, xc
        integer :: k, i, kk

        ftot = fb
        do k = 1, self%nt
            inb = matmul(transpose(ct%jac(:, :, k)), outb(:, k))
            ftot = ftot + inb(5)
            rt = rb(k) + inb(2) + self%sig_psi*(inb(3)*ct%prr(k) + inb(4)*ct%prz(k))
            zt = zb(k) + self%sig_psi*(inb(3)*ct%prz(k) + inb(4)*ct%pzz(k))
            rhob = inb(1) + ct%ca(k)*rt + ct%sa(k)*zt
            r0_b = r0_b + [rt, zt]
            gr = rhob/ct%prho(k)
            psitb = psitb + gr
            r0_b = r0_b - gr*[ct%praw(k), ct%pzraw(k)]
            call psi_eval_t(self, ct%r(k), ct%z(k), &
                            [-gr, self%sig_psi*inb(3), self%sig_psi*inb(4)], cb)
        end do
        ! df = sig_f (dF_data(x_f) + F'(psi_t) dpsi_t)
        psitb = psitb + ftot*self%sig_f*ct%fpsi
        xc = max(0.0_dp, min(1.0_dp, ct%xf))
        i = max(1, min(self%nw - 1, floor(xc/self%hx) + 1))
        t = xc - real(i - 1, dp)*self%hx
        do kk = 0, 5
            cfb = cfb + self%sig_f*ftot*t**kk*self%wf(i, :, kk)
        end do
    end subroutine contour_vjp

    subroutine surface_vjp(self, mpol, ct, sb, cb, cfb, psitb, r0_b)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: mpol
        type(contour_t), intent(in) :: ct
        type(efit_sfl_surface_t), intent(in) :: sb
        real(dp), intent(inout) :: cb(:, :, 0:, 0:), cfb(:), psitb, r0_b(2)
        real(dp) :: g(self%nt), a(self%nt), u(self%nt), jb(self%nt), nu(self%nt)
        real(dp) :: fld(self%nt, 4), fldb(self%nt, 4), jbb(self%nt), ub(self%nt), &
                    nub(self%nt), ab(self%nt), gb(self%nt), outb(5, self%nt), q, qb, fb
        complex(dp) :: e(self%nt), em(self%nt), beta(4)
        real(dp) :: w
        integer :: j, m

        q = sum(ct%out(1, :))/real(self%nt, dp)
        call route_integrand(self, self%cfg%route, ct, g)
        a = ct%out(1, :) + g
        nu = matmul(self%aop, g)
        u = matmul(self%aop, a)/q
        jb = a/q
        fld(:, 1) = ct%r
        fld(:, 2) = ct%z
        fld(:, 3) = -nu/twopi
        fld(:, 4) = ct%out(3, :)
        do j = 1, self%nt
            e(j) = exp(cmplx(0.0_dp, twopi*real(j - 1, dp)/real(self%nt, dp) + u(j), dp))
        end do
        fldb = 0.0_dp
        jbb = 0.0_dp
        ub = 0.0_dp
        em = cmplx(1.0_dp, 0.0_dp, dp)
        do m = 0, mpol
            if (m > 0) em = em*e
            w = merge(1.0_dp, 2.0_dp, m == 0)/real(self%nt, dp)
            beta = w*cmplx([sb%rmnc(m), sb%zmnc(m), sb%vmnc(m), sb%bmnc(m)], &
                           -[sb%rmns(m), sb%zmns(m), sb%vmns(m), sb%bmns(m)], dp)
            do j = 1, 4
                fldb(:, j) = fldb(:, j) + real(beta(j)*em)*jb
                jbb = jbb + real(beta(j)*em)*fld(:, j)
                ub = ub + real(beta(j)*em*cmplx(0.0_dp, real(m, dp), dp))*fld(:, j)*jb
            end do
        end do
        nub = -fldb(:, 3)/twopi
        ab = matmul(transpose(self%aop), ub)/q + jbb/q
        qb = -sum(ub*u)/q - sum(jbb*a)/q**2 - sb%iota/q**2
        gb = matmul(transpose(self%aop), nub) + ab
        outb = 0.0_dp
        outb(1, :) = ab
        outb(3, :) = fldb(:, 4)
        fb = self%sig_f*sb%curr_g
        outb(2, :) = self%sig_f*sb%curr_i/real(self%nt, dp)
        if (self%cfg%route == EFIT_SFL_STRAIGHT_CURRENT) then
            call straight_current_vjp(self, ct, gb, outb, qb, fldb(:, 1), fldb(:, 2), cb, cfb, psitb)
        else
            call route_vjp(self, self%cfg%route, ct, gb, outb, fb, qb)
        end if
        outb(1, :) = outb(1, :) + qb/real(self%nt, dp)
        call contour_vjp(self, ct, outb, fldb(:, 1), fldb(:, 2), fb, cb, cfb, psitb, r0_b)
    end subroutine surface_vjp

    !> Transpose of route_jvp; the q bar is returned in qb (to be spread
    !> over h by the caller).
    pure subroutine route_vjp(self, route, ct, gb, outb, fb, qb)
        type(efit_sfl_t), intent(in) :: self
        integer, intent(in) :: route
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: gb(:)
        real(dp), intent(inout) :: outb(:, :), fb, qb
        real(dp) :: q, mb, r, den, rb, denb, mj, mjb
        real(dp) :: num(self%nt), numb(self%nt)
        integer :: k, n

        n = self%nt
        q = sum(ct%out(1, :))/real(n, dp)
        select case (route)
        case (EFIT_SFL_COVARIANT)
            mb = sum(ct%out(2, :))/real(n, dp)
            r = mb/q
            num = ct%out(2, :) - ct%out(1, :)*r
            den = ct%f + r
            numb = gb/den
            denb = -sum(gb*num)/den**2
            outb(2, :) = outb(2, :) + numb
            outb(1, :) = outb(1, :) - numb*r
            rb = -sum(numb*ct%out(1, :)) + denb
            fb = fb + denb
            outb(2, :) = outb(2, :) + rb/q/real(n, dp)
            qb = qb - rb*mb/q**2
        case (EFIT_SFL_JACOBIAN_B2, EFIT_SFL_JACOBIAN_CONST)
            k = merge(5, 4, route == EFIT_SFL_JACOBIAN_B2)
            mj = sum(ct%out(k, :))/real(n, dp)
            qb = qb + sum(gb*ct%out(k, :))/mj
            outb(k, :) = outb(k, :) + gb*q/mj
            mjb = -sum(gb*q*ct%out(k, :))/mj**2
            outb(k, :) = outb(k, :) + mjb/real(n, dp)
            outb(1, :) = outb(1, :) - gb
        case default
        end select
    end subroutine route_vjp

    !> Pointwise pieces of straight_current: Delta* psi (ds), its derivatives
    !> along R and Z at fixed coefficients (third derivatives of the spline),
    !> a = -jc ds/R^2, and fs = sig_f F'(psi) with its x_file derivative factor.
    pure subroutine straight_current_parts(self, ct, ds, dsr, dsz, a, fs, fxx, i, t)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(out) :: ds(:), dsr(:), dsz(:), a(:), fs, fxx, t
        integer, intent(out) :: i
        real(dp) :: v3(4), xc
        integer :: k

        do k = 1, self%nt
            ds(k) = ct%prr(k) - ct%praw(k)/ct%r(k) + ct%pzz(k)
            call psi_eval3(self, self%cpsi, ct%r(k), ct%z(k), v3)
            dsr(k) = v3(1) - ct%prr(k)/ct%r(k) + ct%praw(k)/ct%r(k)**2 + v3(3)
            dsz(k) = v3(2) - ct%prz(k)/ct%r(k) + v3(4)
        end do
        a = -ct%out(4, :)*ds/ct%r**2
        fs = self%sig_f*ct%fpsi
        ! d^2 F/dx_file^2 (zero beyond [0, 1], where F is held constant)
        xc = max(0.0_dp, min(1.0_dp, ct%xf))
        i = max(1, min(self%nw - 1, floor(xc/self%hx) + 1))
        t = xc - real(i - 1, dp)*self%hx
        fxx = 0.0_dp
        do k = 5, 2, -1
            fxx = real(k*(k - 1), dp)*self%cf(k, i) + t*fxx
        end do
        if (ct%xf /= xc) fxx = 0.0_dp
    end subroutine straight_current_parts

    !> Tangent of the straight_current integrand g = (q a - h <a>)/(<a> - q fs).
    pure subroutine straight_current_jvp(self, ct, dc, dcf, dpsit, dout, drr, dzz, dg)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: dc(0:, 0:, :, :), dcf(0:, :), dpsit, dout(:, :), drr(:), dzz(:)
        real(dp), intent(out) :: dg(:)
        real(dp), dimension(self%nt) :: ds, dsr, dsz, a, da, dds
        real(dp) :: v(6), q, dq, ma, dma, fs, dfs, fxx, t, den, dden, dfx, dummy, hpsi
        integer :: k, i

        call straight_current_parts(self, ct, ds, dsr, dsz, a, fs, fxx, i, t)
        do k = 1, self%nt
            call psi_eval(self, dc, ct%r(k), ct%z(k), v)
            dds(k) = v(4) - v(2)/ct%r(k) + v(6) + dsr(k)*drr(k) + dsz(k)*dzz(k)
        end do
        da = -(dout(4, :)*ds + ct%out(4, :)*dds)/ct%r**2 - 2.0_dp*a*drr/ct%r
        hpsi = self%g%sibry - self%g%simag
        call prof_eval(self, dcf, ct%xf, dummy, dfx)
        dfs = self%sig_f*(dfx + fxx*dpsit/hpsi)/hpsi
        q = sum(ct%out(1, :))/real(self%nt, dp)
        dq = sum(dout(1, :))/real(self%nt, dp)
        ma = sum(a)/real(self%nt, dp)
        dma = sum(da)/real(self%nt, dp)
        den = ma - q*fs
        dden = dma - dq*fs - q*dfs
        dg = (dq*a + q*da - dout(1, :)*ma - ct%out(1, :)*dma)/den &
             - (q*a - ct%out(1, :)*ma)*dden/den**2
    end subroutine straight_current_jvp

    !> Transpose of straight_current_jvp: bars to the kernel outputs (outb),
    !> q (qb), the point positions (rb, zb), the spline coefficients, fpol and
    !> psi_t.
    subroutine straight_current_vjp(self, ct, gb, outb, qb, rb, zb, cb, cfb, psitb)
        type(efit_sfl_t), intent(in) :: self
        type(contour_t), intent(in) :: ct
        real(dp), intent(in) :: gb(:)
        real(dp), intent(inout) :: outb(:, :), qb, rb(:), zb(:), cb(:, :, 0:, 0:), cfb(:), psitb
        real(dp), dimension(self%nt) :: ds, dsr, dsz, a, g, nb, ab, dsb
        real(dp) :: q, ma, fs, fxx, t, den, denb, mab, fsb, hpsi
        integer :: k, i, kk

        call straight_current_parts(self, ct, ds, dsr, dsz, a, fs, fxx, i, t)
        q = sum(ct%out(1, :))/real(self%nt, dp)
        ma = sum(a)/real(self%nt, dp)
        den = ma - q*fs
        g = (q*a - ct%out(1, :)*ma)/den
        nb = gb/den
        denb = -sum(gb*g)/den
        qb = qb + sum(nb*a) - denb*fs
        ab = q*nb
        outb(1, :) = outb(1, :) - nb*ma
        mab = -sum(nb*ct%out(1, :)) + denb
        fsb = -denb*q
        ab = ab + mab/real(self%nt, dp)
        outb(4, :) = outb(4, :) - ab*ds/ct%r**2
        dsb = -ab*ct%out(4, :)/ct%r**2
        rb = rb - 2.0_dp*ab*a/ct%r + dsb*dsr
        zb = zb + dsb*dsz
        do k = 1, self%nt
            call psi_eval_t(self, ct%r(k), ct%z(k), &
                            [0.0_dp, -dsb(k)/ct%r(k), 0.0_dp, dsb(k), 0.0_dp, dsb(k)], cb)
        end do
        hpsi = self%g%sibry - self%g%simag
        psitb = psitb + fsb*self%sig_f*fxx/hpsi**2
        if (ct%xf == max(0.0_dp, min(1.0_dp, ct%xf))) then
            do kk = 1, 5
                cfb = cfb + fsb*self%sig_f/hpsi*real(kk, dp)*t**(kk - 1)*self%wf(i, :, kk)
            end do
        end if
    end subroutine straight_current_vjp

end module efit_sfl
