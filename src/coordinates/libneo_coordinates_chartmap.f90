module libneo_coordinates_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates_base, only: coordinate_system_t, &
                                       UNKNOWN, CYL, VMEC, BOOZER, &
                                       RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM, &
                                       chartmap_from_cyl_ok, chartmap_from_cyl_err_max_iter, &
                                       chartmap_from_cyl_err_singular, &
                                       chartmap_from_cyl_err_out_of_bounds, &
                                       chartmap_from_cyl_err_invalid
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
    use math_constants, only: TWOPI
    use netcdf, only: NF90_CHAR, NF90_GLOBAL, NF90_NOERR, nf90_get_att, nf90_get_var, &
                      nf90_inq_varid, nf90_inquire_attribute
    use fortnum_status, only: fortnum_status_t
    use fortnum_multiroot, only: multiroot_hybrid
    implicit none
    private

    public :: chartmap_coordinate_system_t
    public :: make_chartmap_coordinate_system

    ! Geometric status of the Cartesian inverse cs%invert_cart. These describe the
    ! geometry of the located root ONLY; they carry no confinement-loss meaning. A
    ! caller that needs a loss decision (e.g. SIMPLE's full-orbit pusher) keys it on
    ! its own guiding-centre test, not on this status.
    !   LOCATED      converged interior root, residual at chart tolerance.
    !   CLAMPED_EDGE solver pinned rho to 1 and could not reduce the residual below
    !                tolerance, but the residual is still a small fraction of a
    !                radial cell: the target sits essentially on the last surface.
    !   OUTSIDE      best root has rho clamped to 1 with a residual that is a large
    !                fraction of a radial cell: the target is past the last surface,
    !                where the forward map cannot represent it.
    !   NO_ROOT      singular Jacobian or no convergence from any seed.
    integer, parameter, public :: CHARTMAP_LOCATED = 0
    integer, parameter, public :: CHARTMAP_CLAMPED_EDGE = 1
    integer, parameter, public :: CHARTMAP_OUTSIDE = 2
    integer, parameter, public :: CHARTMAP_NO_ROOT = 3

    ! A located root has residual below an absolute tolerance OR below this fraction
    ! of a radial cell |dx/drho| (scale-correct near the axis and on reactor-size
    ! charts where the absolute floor is meaningless). A clamped-edge root within the
    ! same fraction is CLAMPED_EDGE; a larger residual at rho=1 is OUTSIDE.
    real(dp), parameter :: chartmap_invert_edge_frac = 0.05_dp
    real(dp), parameter :: chartmap_invert_accept_tol = 1.0e-6_dp
    ! Warm guesses with rho this close to 1 take the bracketed 1-D edge solve, which
    ! cannot overshoot rho=1 (bracket [rho_lo, 1]); interior guesses take the 3-D
    ! pseudo-Cartesian multiroot solve.
    real(dp), parameter :: chartmap_invert_edge_band = 0.1_dp

    type, extends(coordinate_system_t) :: chartmap_coordinate_system_t
        type(BatchSplineData3D) :: spl_cart
        type(BatchSplineData3D) :: spl_rz
        logical :: has_spl_rz = .false.
        integer :: nrho = 0
        integer :: ntheta = 0
        integer :: nzeta = 0
        integer :: num_field_periods = 1
        integer :: zeta_convention = UNKNOWN
        integer :: rho_convention = UNKNOWN
        real(dp) :: tol_newton = 1.0e-12_dp
    contains
        procedure :: evaluate_cart => chartmap_evaluate_cart
        procedure :: evaluate_cyl => chartmap_evaluate_cyl
        procedure :: covariant_basis => chartmap_covariant_basis
        procedure :: metric_tensor => chartmap_metric_tensor
        procedure :: christoffel => chartmap_christoffel
        procedure :: from_cyl => chartmap_from_cyl
        procedure :: from_cart => chartmap_from_cart
        procedure :: invert_cart => chartmap_invert_cart_warm
    end type chartmap_coordinate_system_t

    ! Read-only context threaded through the fortnum callbacks: the chart and the
    ! Cartesian target. The solvers never inspect it (they pass it back to the
    ! callback unchanged); it is the active parameter for the implicit-rule
    ! sensitivity and the only state the callbacks need.
    type :: chartmap_invert_ctx_t
        class(chartmap_coordinate_system_t), pointer :: self => null()
        real(dp) :: x_target(3) = 0.0_dp
        real(dp) :: zeta = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: rho_lo = 0.0_dp
        real(dp) :: zeta_period = TWOPI
    end type chartmap_invert_ctx_t

contains

    logical function chartmap_trace_once_enabled()
        character(len=32) :: value
        integer :: len_value, status
        logical, save :: used = .false.

        chartmap_trace_once_enabled = .false.
        if (used) return

        call get_environment_variable('LIBNEO_CHARTMAP_TRACE_ONCE', value, &
                                      length=len_value, status=status)
        if (status /= 0 .or. len_value <= 0) return
        if (trim(value(1:len_value)) == '0') return

        used = .true.
        chartmap_trace_once_enabled = .true.
    end function chartmap_trace_once_enabled

    subroutine chartmap_read_num_field_periods(ncid, num_field_periods)
        integer, intent(in) :: ncid
        integer, intent(out) :: num_field_periods

        integer :: var_num_field_periods
        integer :: var_nfp

        num_field_periods = 1
        if (nf90_inq_varid(ncid, 'num_field_periods', var_num_field_periods) == &
            NF90_NOERR) then
            if (nf90_get_var(ncid, var_num_field_periods, num_field_periods) /= &
                NF90_NOERR) num_field_periods = 1
        else if (nf90_inq_varid(ncid, 'nfp', var_nfp) == NF90_NOERR) then
            if (nf90_get_var(ncid, var_nfp, num_field_periods) /= NF90_NOERR) &
                num_field_periods = 1
        end if
        if (num_field_periods < 1) num_field_periods = 1
    end subroutine chartmap_read_num_field_periods

    subroutine chartmap_read_zeta_convention(ncid, zeta_convention)
        integer, intent(in) :: ncid
        integer, intent(out) :: zeta_convention

        integer :: status
        integer :: att_type
        integer :: att_len
        character(len=:), allocatable :: value

        zeta_convention = UNKNOWN

        status = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'zeta_convention', &
                                        xtype=att_type, len=att_len)
        if (status /= NF90_NOERR) return
        if (att_type /= NF90_CHAR .or. att_len < 1) return

        allocate (character(len=att_len) :: value)
        status = nf90_get_att(ncid, NF90_GLOBAL, 'zeta_convention', value)
        if (status /= NF90_NOERR) return

        select case (trim(value))
        case ('cyl')
            zeta_convention = CYL
        case ('vmec')
            zeta_convention = VMEC
        case ('boozer')
            zeta_convention = BOOZER
        case default
            zeta_convention = UNKNOWN
        end select
    end subroutine chartmap_read_zeta_convention

    subroutine chartmap_read_rho_convention(ncid, rho_convention)
        integer, intent(in) :: ncid
        integer, intent(out) :: rho_convention

        integer :: status
        integer :: att_type
        integer :: att_len
        character(len=:), allocatable :: value

        rho_convention = UNKNOWN

        status = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'rho_convention', &
                                        xtype=att_type, len=att_len)
        if (status /= NF90_NOERR) return
        if (att_type /= NF90_CHAR .or. att_len < 1) return

        allocate (character(len=att_len) :: value)
        status = nf90_get_att(ncid, NF90_GLOBAL, 'rho_convention', value)
        if (status /= NF90_NOERR) return

        select case (trim(value))
        case ('rho_tor', 'vmec', 'vmec_extended', 'vmec_to_wall')
            rho_convention = RHO_TOR
        case ('rho_pol')
            rho_convention = RHO_POL
        case ('psi_tor_norm', 's')
            rho_convention = PSI_TOR_NORM
        case ('psi_pol_norm')
            rho_convention = PSI_POL_NORM
        case default
            rho_convention = UNKNOWN
        end select
    end subroutine chartmap_read_rho_convention

    subroutine chartmap_extend_theta(nrho, ntheta, nzeta, theta, x, y, z, theta_spl, &
                                     x2, y2, z2)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(in) :: theta(ntheta)
        real(dp), intent(in) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
                                z(nrho, ntheta, nzeta)
        real(dp), allocatable, intent(out) :: theta_spl(:)
        real(dp), allocatable, intent(out) :: x2(:, :, :), y2(:, :, :), z2(:, :, :)

        allocate (theta_spl(ntheta + 1))
        theta_spl(1:ntheta) = theta
        theta_spl(ntheta + 1) = theta(1) + TWOPI

        allocate (x2(nrho, ntheta + 1, nzeta))
        allocate (y2(nrho, ntheta + 1, nzeta))
        allocate (z2(nrho, ntheta + 1, nzeta))
        x2(:, 1:ntheta, :) = x
        y2(:, 1:ntheta, :) = y
        z2(:, 1:ntheta, :) = z
        x2(:, ntheta + 1, :) = x(:, 1, :)
        y2(:, ntheta + 1, :) = y(:, 1, :)
        z2(:, ntheta + 1, :) = z(:, 1, :)
    end subroutine chartmap_extend_theta

    subroutine chartmap_extend_zeta(nrho, ntheta, nzeta, zeta, zeta_period, x, y, z, &
                                    zeta_spl, &
                                    x2, y2, z2)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(in) :: zeta(nzeta)
        real(dp), intent(in) :: zeta_period
        real(dp), intent(in) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
                                z(nrho, ntheta, nzeta)
        real(dp), allocatable, intent(out) :: zeta_spl(:)
        real(dp), allocatable, intent(out) :: x2(:, :, :), y2(:, :, :), z2(:, :, :)

        allocate (zeta_spl(nzeta + 1))
        zeta_spl(1:nzeta) = zeta
        zeta_spl(nzeta + 1) = zeta(1) + zeta_period

        ! Periodic closure: the closing slice one field period ahead equals the first.
        ! This is correct because the cart-spline is built in the CO-ROTATING frame
        ! (chartmap_corotate), where the stored quantity is x_tilde = Rot_z(-zeta) x,
        ! genuinely periodic in zeta. The eval rotates back by Rot_z(+zeta).
        allocate (x2(nrho, ntheta, nzeta + 1))
        allocate (y2(nrho, ntheta, nzeta + 1))
        allocate (z2(nrho, ntheta, nzeta + 1))
        x2(:, :, 1:nzeta) = x
        y2(:, :, 1:nzeta) = y
        z2(:, :, 1:nzeta) = z
        x2(:, :, nzeta + 1) = x(:, :, 1)
        y2(:, :, nzeta + 1) = y(:, :, 1)
        z2(:, :, nzeta + 1) = z(:, :, 1)
    end subroutine chartmap_extend_zeta

    ! Rotate the stored Cartesian (x,y) of each zeta slice into the co-rotating frame
    ! x_tilde = Rot_z(-zeta) x, in place. The full Cartesian map is periodic in zeta
    ! only up to the field-period rotation, x(zeta+2pi/nfp) = Rot_z(2pi/nfp) x(zeta), so
    ! it cannot be splined periodically. x_tilde is genuinely periodic in zeta (the
    ! rotation that winds with zeta is divided out), so a periodic spline of x_tilde is
    ! exact and continuous across field-period seams. chartmap_eval_cart rotates back.
    subroutine chartmap_corotate(nrho, ntheta, nzeta, zeta, x, y)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(in) :: zeta(nzeta)
        real(dp), intent(inout) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta)

        real(dp) :: ca, sa, xt(nrho, ntheta), yt(nrho, ntheta)
        integer :: k

        do k = 1, nzeta
            ca = cos(zeta(k))
            sa = sin(zeta(k))
            xt = ca*x(:, :, k) + sa*y(:, :, k)     ! Rot_z(-zeta): [ cos  sin]
            yt = -sa*x(:, :, k) + ca*y(:, :, k)    !               [-sin  cos]
            x(:, :, k) = xt
            y(:, :, k) = yt
        end do
    end subroutine chartmap_corotate

    subroutine chartmap_extend_theta_rz(nrho, ntheta, nzeta, theta, r, zc, theta_spl, &
                                        r2, z2)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(in) :: theta(ntheta)
        real(dp), intent(in) :: r(nrho, ntheta, nzeta), zc(nrho, ntheta, nzeta)
        real(dp), allocatable, intent(out) :: theta_spl(:)
        real(dp), allocatable, intent(out) :: r2(:, :, :), z2(:, :, :)

        allocate (theta_spl(ntheta + 1))
        theta_spl(1:ntheta) = theta
        theta_spl(ntheta + 1) = theta(1) + TWOPI

        allocate (r2(nrho, ntheta + 1, nzeta))
        allocate (z2(nrho, ntheta + 1, nzeta))
        r2(:, 1:ntheta, :) = r
        z2(:, 1:ntheta, :) = zc
        r2(:, ntheta + 1, :) = r(:, 1, :)
        z2(:, ntheta + 1, :) = zc(:, 1, :)
    end subroutine chartmap_extend_theta_rz

    subroutine chartmap_extend_zeta_rz(nrho, ntheta, nzeta, zeta, zeta_period, r, zc, &
                                       zeta_spl, &
                                       r2, z2)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(in) :: zeta(nzeta)
        real(dp), intent(in) :: zeta_period
        real(dp), intent(in) :: r(nrho, ntheta, nzeta), zc(nrho, ntheta, nzeta)
        real(dp), allocatable, intent(out) :: zeta_spl(:)
        real(dp), allocatable, intent(out) :: r2(:, :, :), z2(:, :, :)

        allocate (zeta_spl(nzeta + 1))
        zeta_spl(1:nzeta) = zeta
        zeta_spl(nzeta + 1) = zeta(1) + zeta_period

        allocate (r2(nrho, ntheta, nzeta + 1))
        allocate (z2(nrho, ntheta, nzeta + 1))
        r2(:, :, 1:nzeta) = r
        z2(:, :, 1:nzeta) = zc
        r2(:, :, nzeta + 1) = r(:, :, 1)
        z2(:, :, nzeta + 1) = zc(:, :, 1)
    end subroutine chartmap_extend_zeta_rz

    subroutine chartmap_eval_cart(self, u, vals)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: vals(3)

        real(dp) :: rz(2)
        real(dp) :: phi
        real(dp) :: phi_mod
        real(dp) :: zeta_period
        real(dp) :: u_eval(3)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: ang, ca, sa, vx, vy

        if (self%has_spl_rz) then
            x_min = self%spl_rz%x_min
            x_max = self%spl_rz%x_min + self%spl_rz%h_step * &
                    real(self%spl_rz%num_points - 1, dp)
        else
            x_min = self%spl_cart%x_min
            x_max = self%spl_cart%x_min + self%spl_cart%h_step * &
                    real(self%spl_cart%num_points - 1, dp)
        end if

        u_eval = u
        u_eval(1) = max(x_min(1), min(u_eval(1), x_max(1)))
        u_eval(2) = modulo(u_eval(2) - x_min(2), TWOPI) + x_min(2)
        zeta_period = TWOPI/real(self%num_field_periods, dp)
        u_eval(3) = modulo(u_eval(3) - x_min(3), zeta_period) + x_min(3)

        if (self%has_spl_rz) then
            phi = u_eval(3)
            phi_mod = modulo(phi - self%spl_rz%x_min(3), zeta_period) + &
                      self%spl_rz%x_min(3)
            call evaluate_batch_splines_3d(self%spl_rz, u_eval, rz)
            vals(1) = rz(1)*cos(phi_mod)
            vals(2) = rz(1)*sin(phi_mod)
            vals(3) = rz(2)
        else
            call evaluate_batch_splines_3d(self%spl_cart, u_eval, vals)
            ! The cart spline holds the co-rotating position x_tilde (periodic in zeta);
            ! rotate it back to the lab frame, x = Rot_z(zeta) x_tilde. Use the full
            ! input zeta u(3), not the wedge-reduced one: x_tilde is periodic, so this is
            ! exact and continuous across field-period seams for any zeta.
            ang = u(3)
            ca = cos(ang); sa = sin(ang)
            vx = vals(1); vy = vals(2)
            vals(1) = ca*vx - sa*vy
            vals(2) = sa*vx + ca*vy
        end if
    end subroutine chartmap_eval_cart

    subroutine chartmap_eval_cart_der(self, u, vals, dvals)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: vals(3)
        real(dp), intent(out) :: dvals(3, 3)

        real(dp) :: rz(2)
        real(dp) :: drz(3, 2)
        real(dp) :: phi
        real(dp) :: phi_mod
        real(dp) :: cph, sph
        real(dp) :: zeta_period
        real(dp) :: u_eval(3)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: ang, ca, sa, rx, ry
        integer :: k

        if (self%has_spl_rz) then
            x_min = self%spl_rz%x_min
            x_max = self%spl_rz%x_min + self%spl_rz%h_step * &
                    real(self%spl_rz%num_points - 1, dp)
        else
            x_min = self%spl_cart%x_min
            x_max = self%spl_cart%x_min + self%spl_cart%h_step * &
                    real(self%spl_cart%num_points - 1, dp)
        end if

        u_eval = u
        u_eval(1) = max(x_min(1), min(u_eval(1), x_max(1)))
        u_eval(2) = modulo(u_eval(2) - x_min(2), TWOPI) + x_min(2)
        zeta_period = TWOPI/real(self%num_field_periods, dp)
        u_eval(3) = modulo(u_eval(3) - x_min(3), zeta_period) + x_min(3)

        if (self%has_spl_rz) then
            phi = u_eval(3)
            phi_mod = modulo(phi - self%spl_rz%x_min(3), zeta_period) + &
                      self%spl_rz%x_min(3)
            cph = cos(phi_mod)
            sph = sin(phi_mod)
            call evaluate_batch_splines_3d_der(self%spl_rz, u_eval, rz, drz)

            vals(1) = rz(1)*cph
            vals(2) = rz(1)*sph
            vals(3) = rz(2)

            dvals(1, 1) = drz(1, 1)*cph
            dvals(2, 1) = drz(1, 1)*sph
            dvals(3, 1) = drz(1, 2)

            dvals(1, 2) = drz(2, 1)*cph
            dvals(2, 2) = drz(2, 1)*sph
            dvals(3, 2) = drz(2, 2)

            dvals(1, 3) = drz(3, 1)*cph - rz(1)*sph
            dvals(2, 3) = drz(3, 1)*sph + rz(1)*cph
            dvals(3, 3) = drz(3, 2)
        else
            call evaluate_batch_splines_3d_der(self%spl_cart, u_eval, vals, dvals)
            ! evaluate_batch_splines_3d_der returns dvals(deriv_dim, quantity) =
            ! d x_tilde_quantity / d u_dim. Transpose to the documented covariant-basis
            ! convention dvals(cart_component, coord_index) = d x_tilde_i / d u_k.
            dvals = transpose(dvals)
            ! Co-rotating reconstruction x = Rot_z(zeta) x_tilde (see chartmap_eval_cart).
            ! The rho and theta basis columns rotate; the zeta column also picks up the
            ! derivative of the rotation, Rot_z'(zeta) x_tilde = Rot_z(zeta) (J x_tilde)
            ! with J x_tilde = (-x_tilde_y, x_tilde_x, 0). Add J x_tilde to the zeta
            ! column first (using the unrotated x_tilde), then rotate value and every
            ! column by Rot_z(zeta).
            dvals(1, 3) = dvals(1, 3) - vals(2)
            dvals(2, 3) = dvals(2, 3) + vals(1)
            ang = u(3)
            ca = cos(ang); sa = sin(ang)
            rx = vals(1); ry = vals(2)
            vals(1) = ca*rx - sa*ry
            vals(2) = sa*rx + ca*ry
            do k = 1, 3
                rx = dvals(1, k); ry = dvals(2, k)
                dvals(1, k) = ca*rx - sa*ry
                dvals(2, k) = sa*rx + ca*ry
            end do
        end if
    end subroutine chartmap_eval_cart_der

    subroutine make_chartmap_coordinate_system(cs, filename)
        use libneo_coordinates_validator, only: validate_chartmap_file
        class(coordinate_system_t), allocatable, intent(out) :: cs
        character(len=*), intent(in) :: filename

        allocate (chartmap_coordinate_system_t :: cs)

        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            call initialize_chartmap(ccs, filename)
        class default
            error stop "make_chartmap_coordinate_system: allocation failed"
        end select
    end subroutine make_chartmap_coordinate_system

    subroutine initialize_chartmap(ccs, filename)
        use libneo_coordinates_validator, only: validate_chartmap_file
        type(chartmap_coordinate_system_t), intent(inout) :: ccs
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: ierr
        integer :: len_rho, len_theta, len_zeta
        integer :: nrho, ntheta, nzeta
        character(len=2048) :: message
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        real(dp), allocatable :: pos_batch(:, :, :, :)
        real(dp), allocatable :: pos_rz(:, :, :, :)
        real(dp), allocatable :: r(:, :, :), zc(:, :, :)
        real(dp), allocatable :: theta_spl(:), zeta_spl(:)
        real(dp), allocatable :: x_th(:, :, :), y_th(:, :, :), z_th(:, :, :)
        real(dp), allocatable :: x_spl(:, :, :), y_spl(:, :, :), z_spl(:, :, :)
        real(dp), allocatable :: r_th(:, :, :), z_th_rz(:, :, :)
        real(dp), allocatable :: r_spl(:, :, :), z_spl_rz(:, :, :)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: x_min_rz(3), x_max_rz(3)
        real(dp) :: zeta_period
        logical :: periodic_cart(3)
        logical :: periodic_rz(3)
        integer :: order(3)
        integer :: num_field_periods
        integer :: zeta_convention
        integer :: rho_convention

        call validate_chartmap_file(trim(filename), ierr, message)
        if (ierr /= 0) then
            print *, "initialize_chartmap: invalid chartmap file"
            print *, trim(message)
            error stop
        end if

        call nc_open(trim(filename), ncid)

        call nc_inq_dim(ncid, 'rho', len_rho)
        call nc_inq_dim(ncid, 'theta', len_theta)
        call nc_inq_dim(ncid, 'zeta', len_zeta)

        nrho = len_rho
        ntheta = len_theta
        nzeta = len_zeta

        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        call nc_get(ncid, 'rho', rho)
        call nc_get(ncid, 'theta', theta)
        call nc_get(ncid, 'zeta', zeta)

        allocate (x(nrho, ntheta, nzeta))
        allocate (y(nrho, ntheta, nzeta))
        allocate (z(nrho, ntheta, nzeta))
        call nc_get(ncid, 'x', x)
        call nc_get(ncid, 'y', y)
        call nc_get(ncid, 'z', z)
        call chartmap_read_num_field_periods(ncid, num_field_periods)
        call chartmap_read_zeta_convention(ncid, zeta_convention)
        call chartmap_read_rho_convention(ncid, rho_convention)
        call nc_close(ncid)

        if (num_field_periods == 1 .and. nzeta >= 2) then
            block
                real(dp) :: dz, period_est, nfp_est_real
                integer :: nfp_est

                dz = zeta(2) - zeta(1)
                period_est = (zeta(nzeta) - zeta(1)) + dz
                if (period_est > 0.0_dp) then
                    nfp_est_real = TWOPI/period_est
                    nfp_est = max(1, nint(nfp_est_real))
                    if (abs(nfp_est_real - real(nfp_est, dp)) < 0.25_dp) then
                        num_field_periods = nfp_est
                    end if
                end if
            end block
        end if

        order = [3, 3, 3]

        periodic_cart = [.false., .true., .false.]
        periodic_rz = [.false., .true., .false.]

        ccs%num_field_periods = num_field_periods
        ccs%zeta_convention = zeta_convention
        ccs%rho_convention = rho_convention
        zeta_period = TWOPI/real(num_field_periods, dp)

        if (zeta_convention /= CYL .and. zeta_convention /= VMEC .and. &
            zeta_convention /= BOOZER) then
            print *, "initialize_chartmap: unsupported zeta_convention"
            print *, "  must be cyl, vmec, or boozer"
            error stop
        end if

        if (zeta_convention == CYL .or. zeta_convention == VMEC) then
            if (nzeta < 2) then
                print *, "initialize_chartmap: nzeta must be >= 2 for cyl convention"
                error stop
            end if

            allocate (r(nrho, ntheta, nzeta))
            allocate (zc(nrho, ntheta, nzeta))
            r = sqrt(x**2 + y**2)
            zc = z

            call chartmap_extend_theta_rz(nrho, ntheta, nzeta, theta, r, zc, &
                                          theta_spl, &
                                          r_th, z_th_rz)

            periodic_rz(3) = .true.
            call chartmap_extend_zeta_rz(nrho, ntheta + 1, nzeta, zeta, zeta_period, &
                                         r_th, z_th_rz, zeta_spl, r_spl, z_spl_rz)

            x_min_rz = [rho(1), theta_spl(1), zeta_spl(1)]
            x_max_rz = [rho(nrho), theta_spl(size(theta_spl)), zeta_spl(size(zeta_spl))]

            allocate (pos_rz(nrho, size(theta_spl), size(zeta_spl), 2))
            pos_rz(:, :, :, 1) = r_spl
            pos_rz(:, :, :, 2) = z_spl_rz

            call construct_batch_splines_3d(x_min_rz, x_max_rz, pos_rz, order, &
                                            periodic_rz, &
                                            ccs%spl_rz)
            ccs%has_spl_rz = .true.
        else
            if (nzeta < 2) then
                print *, "initialize_chartmap: nzeta must be >= 2"
                error stop
            end if

            call chartmap_extend_theta(nrho, ntheta, nzeta, theta, x, y, z, theta_spl, &
                                       x_th, y_th, z_th)

            ! Spline the position in the co-rotating frame x_tilde = Rot_z(-zeta) x.
            ! The full Cartesian map is periodic in zeta only up to the field-period
            ! rotation; x_tilde divides that rotation out and is genuinely periodic, so
            ! the zeta spline is periodic and exact across field-period seams.
            ! chartmap_eval_cart / _der rotate back by Rot_z(+zeta). theta stays periodic.
            call chartmap_corotate(nrho, ntheta + 1, nzeta, zeta, x_th, y_th)
            periodic_cart(3) = .true.
            call chartmap_extend_zeta(nrho, ntheta + 1, nzeta, zeta, zeta_period, &
                                      x_th, y_th, z_th, zeta_spl, x_spl, y_spl, z_spl)

            x_min = [rho(1), theta_spl(1), zeta_spl(1)]
            x_max = [rho(nrho), theta_spl(size(theta_spl)), zeta_spl(size(zeta_spl))]

            allocate (pos_batch(nrho, size(theta_spl), size(zeta_spl), 3))
            pos_batch(:, :, :, 1) = x_spl
            pos_batch(:, :, :, 2) = y_spl
            pos_batch(:, :, :, 3) = z_spl

            call construct_batch_splines_3d(x_min, x_max, pos_batch, order, &
                                            periodic_cart, &
                                            ccs%spl_cart)
            ccs%has_spl_rz = .false.
        end if

        ccs%nrho = nrho
        ccs%ntheta = ntheta
        ccs%nzeta = nzeta
    end subroutine initialize_chartmap

    subroutine chartmap_evaluate_cart(self, u, x)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        call chartmap_eval_cart(self, u, x)
    end subroutine chartmap_evaluate_cart

    subroutine chartmap_evaluate_cyl(self, u, x)
        use cylindrical_cartesian, only: cart_to_cyl
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xcart(3)
        real(dp) :: rz(2)
        real(dp) :: zeta_period

        if (self%has_spl_rz) then
            zeta_period = TWOPI/real(self%num_field_periods, dp)
            call evaluate_batch_splines_3d(self%spl_rz, u, rz)
            x(1) = rz(1)
            x(2) = modulo(u(3), zeta_period)
            x(3) = rz(2)
        else
            call chartmap_eval_cart(self, u, xcart)
            call cart_to_cyl(xcart, x)
        end if
    end subroutine chartmap_evaluate_cyl

    subroutine chartmap_covariant_basis(self, u, e_cov)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        real(dp) :: vals(3)
        real(dp) :: dvals(3, 3)

        call chartmap_eval_cart_der(self, u, vals, dvals)
        e_cov = dvals
    end subroutine chartmap_covariant_basis

    subroutine chartmap_metric_tensor(self, u, g, ginv, sqrtg)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg

        real(dp) :: e_cov(3, 3)
        real(dp) :: det
        integer :: i, j

        call self%covariant_basis(u, e_cov)

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
            end do
        end do

        det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
              - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
              + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

        sqrtg = sqrt(abs(det))

        ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
        ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
        ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
        ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
        ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
        ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
        ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
        ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
        ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det
    end subroutine chartmap_metric_tensor

    subroutine chartmap_christoffel(self, u, Gamma)
        !> Christoffel symbols of the second kind from central differences of
        !> chartmap_metric_tensor; consistent with the splined chart metric.
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)

        call self%christoffel_fd(u, Gamma)
    end subroutine chartmap_christoffel

    subroutine chartmap_from_cyl(self, xcyl, u, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: x_target(3), rho_theta(2)
        real(dp) :: zeta
        real(dp) :: zeta_period

        x_target(1) = xcyl(1)*cos(xcyl(2))
        x_target(2) = xcyl(1)*sin(xcyl(2))
        x_target(3) = xcyl(3)

        zeta_period = TWOPI/real(self%num_field_periods, dp)
        if (self%zeta_convention == CYL .or. self%zeta_convention == VMEC) then
            zeta = modulo(xcyl(2), zeta_period)
            call newton_slice(self, x_target, zeta, rho_theta, ierr)
            u(1) = min(max(rho_theta(1), 0.0_dp), 1.0_dp)
            u(2) = modulo(rho_theta(2), TWOPI)
            u(3) = zeta
            return
        end if

        call self%from_cart(x_target, u, ierr)
    end subroutine chartmap_from_cyl

    subroutine chartmap_from_cart(self, x, u, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        call chartmap_invert_cart(self, x, u, ierr)
    end subroutine chartmap_from_cart

    ! Warm Cartesian inverse x -> logical u = (rho, theta, phi) from a carried guess,
    ! the single-point inverse the full-orbit pusher needs. The iteration runs in the
    ! pseudo-Cartesian chart w = (X, Y, phi) = (rho cos theta, rho sin theta, phi),
    ! which is regular through the magnetic axis (the polar column dx/dtheta ~ rho
    ! makes det(Jc) -> 0 at rho = 0; the (X, Y) Jacobian stays finite). Two solver
    ! paths, both fortnum:
    !   - warm guess near the edge (rho_guess within an edge band of 1): a bracketed
    !     1-D radial solve on [rho_lo, 1], so the radius cannot overshoot past 1.
    !   - otherwise: a 3-D multiroot_hybrid in w with the analytic Jacobian built
    !     from covariant_basis. The callback rejects any rho > 1 trial (large
    !     residual) so the line search cannot accept a clamped-edge point.
    ! On warm-path failure, fall back to the cold multi-seed chartmap_invert_cart.
    ! geom_status is GEOMETRIC ONLY (located / clamped-edge / outside / no-root); it
    ! never encodes a confinement loss.
    subroutine chartmap_invert_cart_warm(self, x, u_guess, u, geom_status)
        class(chartmap_coordinate_system_t), intent(in), target :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: u_guess(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: geom_status

        real(dp) :: u_try(3), res_try, res_best, rscale, zeta_period
        logical :: converged, have_root

        zeta_period = TWOPI/real(self%num_field_periods, dp)
        res_best = huge(1.0_dp)
        have_root = .false.

        ! Try each solver in turn; keep the lowest-residual root seen so far. A path
        ! that does not converge still contributes its best clamped root (e.g. a
        ! past-edge target pinned at rho=1), so classification can report OUTSIDE
        ! rather than a spurious no-root. Stop early once a path truly converges.
        ! The toroidal coordinate is seeded from the in-wedge geometric angle inside
        ! each solver: the carried Boozer phi goes a full period stale across a
        ! field-period seam, whereas atan2(y, x) is always in-wedge and differs from
        ! logical phi only by the Boozer shift O(0.1 rad). For the cyl/vmec charts
        ! zeta IS the geometric angle, so the seed is exact.
        if (u_guess(1) >= 1.0_dp - chartmap_invert_edge_band) then
            call chartmap_invert_edge_1d(self, x, u_guess, u_try, res_try, converged)
            call keep_best(u_try, res_try, u, res_best, have_root)
            if (converged) go to 100
        end if

        call chartmap_invert_interior_3d(self, x, u_guess, u_try, res_try, converged)
        call keep_best(u_try, res_try, u, res_best, have_root)
        if (converged) go to 100

        call chartmap_invert_cold(self, x, u_try, res_try, converged)
        call keep_best(u_try, res_try, u, res_best, have_root)

100     continue
        if (.not. have_root) then
            u = 0.0_dp
            geom_status = CHARTMAP_NO_ROOT
            return
        end if

        u(1) = min(max(u(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(u(2), TWOPI)
        ! The rz chart is wedge-native (its forward map reduces zeta internally);
        ! the co-rotating cart chart is only 2pi-periodic in zeta, so wrapping into
        ! the wedge would move the returned point by a field-period rotation.
        if (self%has_spl_rz) then
            u(3) = modulo(u(3), zeta_period)
        else
            u(3) = modulo(u(3), TWOPI)
        end if

        rscale = chartmap_radial_scale(self, u)
        geom_status = chartmap_classify_root(u(1), res_best, rscale)
    end subroutine chartmap_invert_cart_warm

    ! Retain the lower-residual root. A finite-residual candidate sets have_root so a
    ! non-converged clamped-edge root still survives for OUTSIDE classification.
    pure subroutine keep_best(u_try, res_try, u_best, res_best, have_root)
        real(dp), intent(in) :: u_try(3), res_try
        real(dp), intent(inout) :: u_best(3), res_best
        logical, intent(inout) :: have_root
        if (res_try == res_try .and. res_try < res_best) then
            u_best = u_try
            res_best = res_try
            have_root = .true.
        end if
    end subroutine keep_best

    ! Bracketed 1-D edge solve along the seeded ray: hold (theta, phi) at the
    ! warm-guess values and solve the scalar radial residual for rho with the bracket
    ! [rho_lo, 1] enforced so the radius cannot overshoot past 1. The driver is
    ! fortnum's multiroot_hybrid at n = 1 (Newton + line search on the single radial
    ! variable); the analytic 1x1 Jacobian is the outward radial sensitivity from one
    ! covariant_basis call (no finite differences). The bracket lives in the
    ! callback: a trial rho > 1 is penalized (large outward residual) so the line
    ! search backtracks, and rho < rho_lo is clamped, so a marker leaving the plasma
    ! stalls cleanly on the clamped edge instead of running away.
    subroutine chartmap_invert_edge_1d(self, x, u_guess, u, res_norm, converged)
        class(chartmap_coordinate_system_t), intent(in), target :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: u_guess(3)
        real(dp), intent(out) :: u(3)
        real(dp), intent(out) :: res_norm
        logical, intent(out) :: converged

        type(chartmap_invert_ctx_t) :: ctx
        type(fortnum_status_t) :: status
        real(dp) :: rho0(1), rho(1), xc(3)

        ctx%self => self
        ctx%x_target = x
        ctx%theta = u_guess(2)
        ctx%zeta = u_guess(3)
        if (self%has_spl_rz) ctx%zeta = atan2(x(2), x(1))
        ctx%rho_lo = max(0.0_dp, u_guess(1) - chartmap_invert_edge_band)

        rho0(1) = min(max(u_guess(1), ctx%rho_lo), 1.0_dp)
        call multiroot_hybrid(chartmap_radial_fdf, 1, rho0, rho, status, &
                              xtol=self%tol_newton, ftol=self%tol_newton, &
                              max_iter=30, ctx=ctx)

        u = [min(max(rho(1), 0.0_dp), 1.0_dp), ctx%theta, ctx%zeta]
        call chartmap_eval_cart(self, u, xc)
        res_norm = sqrt(sum((x - xc)**2))
        ! Converged only when the residual is genuinely small along this fixed ray.
        ! A point whose true (theta, phi) differ from the guess cannot be reached on
        ! the held ray; its residual stays large and the caller falls to the 3-D
        ! solve. The clamped-edge root is still returned for OUTSIDE classification.
        converged = chartmap_root_located(self, u, res_norm)
    end subroutine chartmap_invert_edge_1d

    ! Interior 3-D solve in the pseudo-Cartesian chart w = (X, Y, phi). The callback
    ! returns F = x(w) - x_target and the analytic Jacobian dF/dw built from
    ! covariant_basis; multiroot_hybrid's Newton + line search drives it to a root.
    ! Seed w from the warm guess with phi reseeded to the in-wedge geometric angle.
    subroutine chartmap_invert_interior_3d(self, x, u_guess, u, res_norm, converged)
        class(chartmap_coordinate_system_t), intent(in), target :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: u_guess(3)
        real(dp), intent(out) :: u(3)
        real(dp), intent(out) :: res_norm
        logical, intent(out) :: converged

        type(chartmap_invert_ctx_t) :: ctx
        type(fortnum_status_t) :: status
        real(dp) :: w0(3), w(3), phi_seed, xc(3)

        phi_seed = u_guess(3)
        if (self%has_spl_rz) phi_seed = atan2(x(2), x(1))

        ctx%self => self
        ctx%x_target = x
        ctx%zeta_period = TWOPI/real(self%num_field_periods, dp)

        w0(1) = u_guess(1)*cos(u_guess(2))
        w0(2) = u_guess(1)*sin(u_guess(2))
        w0(3) = phi_seed

        call multiroot_hybrid(chartmap_invert_fdf, 3, w0, w, status, &
                              xtol=self%tol_newton, ftol=self%tol_newton, &
                              max_iter=60, ctx=ctx)

        call chartmap_w_to_u(w, u)
        call chartmap_eval_cart(self, u, xc)
        res_norm = sqrt(sum((x - xc)**2))
        ! Converged when the residual is at the chart tolerance or a small fraction of
        ! a radial cell; a larger residual (a genuine stall or a past-edge target the
        ! line search pinned at rho=1) falls through to the cold multi-seed, while the
        ! best root is still returned for classification.
        converged = chartmap_root_located(self, u, res_norm)
    end subroutine chartmap_invert_interior_3d

    ! Cold multi-seed fallback: reuse the robust zeta-swept chartmap_invert_cart so a
    ! warm guess stale across a seam still locates. Returns the residual for
    ! classification; converged flags a genuine locate.
    subroutine chartmap_invert_cold(self, x, u, res_norm, converged)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u(3)
        real(dp), intent(out) :: res_norm
        logical, intent(out) :: converged

        real(dp) :: xc(3)
        integer :: ierr

        converged = .false.
        res_norm = huge(1.0_dp)
        u = 0.0_dp
        call chartmap_invert_cart(self, x, u, ierr)
        if (ierr /= chartmap_from_cyl_ok) return
        call chartmap_eval_cart(self, u, xc)
        res_norm = sqrt(sum((x - xc)**2))
        converged = chartmap_root_located(self, u, res_norm)
    end subroutine chartmap_invert_cold

    ! A root is located when its Cartesian residual is at the absolute tolerance or a
    ! small fraction of the local radial cell |dx/drho|. The relative test is what
    ! keeps this scale-correct near the axis and on reactor-size charts.
    logical function chartmap_root_located(self, u, res_norm) result(located)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: res_norm
        real(dp) :: rscale
        rscale = chartmap_radial_scale(self, u)
        located = res_norm < chartmap_invert_accept_tol .or. &
                  (rscale > 0.0_dp .and. res_norm < chartmap_invert_edge_frac*rscale)
    end function chartmap_root_located

    ! fortnum multiroot callback: F(w) = x(w) - x_target with the analytic Jacobian
    ! dF/dw = pseudo-Cartesian Jacobian Jw built from covariant_basis. A trial with
    ! rho = sqrt(X^2 + Y^2) > 1 is past the chart edge, where the forward map clamps
    ! rho and would report a spuriously small residual on the flat clamped region; we
    ! penalize it (large outward F) so the line search backtracks instead of
    ! accepting the clamped-edge point. This replicates the over-edge rejection the
    ! full-orbit Newton hand-rolled.
    subroutine chartmap_invert_fdf(w, f, jac, ctx)
        real(dp), intent(in) :: w(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx

        real(dp) :: u(3), e_cov(3, 3), jw(3, 3), xc(3), cth, sth, rho
        real(dp), parameter :: edge_penalty = 1.0e6_dp
        integer :: i

        f = 0.0_dp
        jac = 0.0_dp
        do i = 1, 3
            jac(i, i) = 1.0_dp
        end do
        if (.not. present(ctx)) return

        select type (ctx)
        type is (chartmap_invert_ctx_t)
            call chartmap_w_to_u(w, u)
            rho = u(1)
            if (rho > 1.0_dp) then
                ! Penalize past-edge trials along the outward radial direction so the
                ! Armijo step shrinks; keep the identity Jacobian as a safe descent.
                f = 0.0_dp
                f(1) = edge_penalty*(rho - 1.0_dp)
                return
            end if
            call chartmap_eval_cart(ctx%self, u, xc)
            f = xc - ctx%x_target
            call chartmap_covariant_basis(ctx%self, u, e_cov)
            call chartmap_pseudocart_basis(u, e_cov, jw, cth, sth)
            jac = jw
        end select
    end subroutine chartmap_invert_fdf

    ! fortnum n=1 multiroot callback for the bracketed edge solve: the residual
    ! projected on the outward radial direction along the fixed (theta, phi) ray and
    ! its 1x1 Jacobian. f = (x(rho) - x_target) . rhat, rhat = e_rho/|e_rho|, with
    ! e_rho = dx/drho the covariant radial column (analytic); jac(1,1) = e_rho . rhat
    ! is the dominant radial sensitivity, the off-ray curvature being second order.
    ! The bracket [rho_lo, 1] is enforced here: rho > 1 is penalized (large outward
    ! residual) so the line search backtracks, replicating the upper bracket bound.
    subroutine chartmap_radial_fdf(w, f, jac, ctx)
        real(dp), intent(in) :: w(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx

        real(dp) :: u(3), xc(3), e_cov(3, 3), e_rho(3), rhat(3), e_norm, rho
        real(dp), parameter :: edge_penalty = 1.0e6_dp

        f = 0.0_dp
        jac(1, 1) = 1.0_dp
        if (.not. present(ctx)) return

        select type (ctx)
        type is (chartmap_invert_ctx_t)
            rho = w(1)
            if (rho > 1.0_dp) then
                f(1) = edge_penalty*(rho - 1.0_dp)
                return
            end if
            rho = max(rho, ctx%rho_lo)
            u = [rho, ctx%theta, ctx%zeta]
            call chartmap_eval_cart(ctx%self, u, xc)
            call chartmap_covariant_basis(ctx%self, u, e_cov)
            e_rho = e_cov(:, 1)
            e_norm = sqrt(sum(e_rho**2))
            if (e_norm < 1.0e-30_dp) then
                rhat = 0.0_dp
            else
                rhat = e_rho/e_norm
            end if
            f(1) = dot_product(xc - ctx%x_target, rhat)
            jac(1, 1) = dot_product(e_rho, rhat)
        end select
    end subroutine chartmap_radial_fdf

    ! Pseudo-Cartesian near-axis chart basis. w = (X, Y, phi) = (rho cos theta,
    ! rho sin theta, phi). The polar covariant column e_cov(:,2) = dx/dtheta ~ rho
    ! vanishes at the axis (det -> 0); the (X, Y) columns built by the chain rule of
    ! X, Y stay finite because dx/dtheta ~ rho cancels the 1/rho. Returns the regular
    ! chart Jacobian jw(a,i) = dx_a/dw_i and the trig used to map field components.
    ! Map geometry shared by the inverse solve and the field assembly.
    subroutine chartmap_pseudocart_basis(u, e_cov, jw, cth, sth)
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: e_cov(3, 3)
        real(dp), intent(out) :: jw(3, 3)
        real(dp), intent(out) :: cth, sth
        real(dp) :: rho
        integer :: a

        rho = max(u(1), 1.0e-30_dp)
        cth = cos(u(2))
        sth = sin(u(2))
        do a = 1, 3
            jw(a, 1) = e_cov(a, 1)*cth - e_cov(a, 2)*(sth/rho)   ! dx/dX
            jw(a, 2) = e_cov(a, 1)*sth + e_cov(a, 2)*(cth/rho)   ! dx/dY
            jw(a, 3) = e_cov(a, 3)                               ! dx/dphi
        end do
    end subroutine chartmap_pseudocart_basis

    ! Pseudo-Cartesian w = (X, Y, phi) -> logical u = (rho, theta, phi). rho >= 0 and
    ! the atan2 branch make the axis an ordinary point (no reflect hack).
    pure subroutine chartmap_w_to_u(w, u)
        real(dp), intent(in) :: w(3)
        real(dp), intent(out) :: u(3)
        u(1) = sqrt(w(1)**2 + w(2)**2)
        u(2) = atan2(w(2), w(1))
        u(3) = w(3)
    end subroutine chartmap_w_to_u

    ! Length of one unit-rho radial step |dx/drho| = |e_cov(:,1)|, the chart scale
    ! that sizes a residual: a residual a small fraction of a radial cell means the
    ! target is essentially located; a large fraction is an interior stall.
    real(dp) function chartmap_radial_scale(self, u) result(s)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp) :: e_cov(3, 3)
        call chartmap_covariant_basis(self, u, e_cov)
        s = sqrt(e_cov(1, 1)**2 + e_cov(2, 1)**2 + e_cov(3, 1)**2)
    end function chartmap_radial_scale

    ! Classify a converged root from its radius and residual against the local
    ! radial cell. Interior at tolerance -> LOCATED; an edge-clamped root within a
    ! cell fraction -> CLAMPED_EDGE; a larger residual at rho = 1 -> OUTSIDE.
    pure integer function chartmap_classify_root(rho, res_norm, rscale) result(status)
        real(dp), intent(in) :: rho, res_norm, rscale
        logical :: located, at_edge

        located = res_norm < chartmap_invert_accept_tol .or. &
                  (rscale > 0.0_dp .and. res_norm < chartmap_invert_edge_frac*rscale)
        at_edge = rho >= 1.0_dp - chartmap_invert_edge_frac

        if (located .and. .not. at_edge) then
            status = CHARTMAP_LOCATED
        else if (located .and. at_edge) then
            status = CHARTMAP_CLAMPED_EDGE
        else if (at_edge) then
            status = CHARTMAP_OUTSIDE
        else
            status = CHARTMAP_NO_ROOT
        end if
    end function chartmap_classify_root

    subroutine chartmap_invert_cart(self, x_target, u, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: u_try(3)
        real(dp) :: res_norm, best_res
        real(dp) :: zeta_period
        real(dp) :: phi
        real(dp) :: phi_shift, ca, sa
        real(dp) :: x_wedge(3), x_round(3)
        real(dp) :: accept_tol
        real(dp) :: zeta_seed(5)
        logical :: trace
        integer :: k
        integer :: ierr_seed, ierr_try

        ierr = chartmap_from_cyl_err_max_iter
        best_res = huge(1.0_dp)
        u = 0.0_dp

        zeta_period = TWOPI/real(self%num_field_periods, dp)
        if (.not. self%has_spl_rz) then
            trace = chartmap_trace_once_enabled()
            ! The co-rotating cart chart is field-period symmetric,
            ! x(rho, theta, zeta + k*zeta_period) = Rot_z(k*zeta_period) x(rho,
            ! theta, zeta), but the lab-frame position itself is only 2pi-periodic
            ! in zeta, so a target may sit in any field-period sector. Rotate it
            ! back by the whole field periods of its geometric angle: the rotated
            ! target has its root near the stored wedge, where the seed scan
            ! resolves, and the solved zeta shifts back exactly by the symmetry.
            phi_shift = zeta_period*floor(modulo(atan2(x_target(2), x_target(1)), &
                                                 TWOPI)/zeta_period)
            ca = cos(phi_shift)
            sa = sin(phi_shift)
            x_wedge(1) = ca*x_target(1) + sa*x_target(2)
            x_wedge(2) = -sa*x_target(1) + ca*x_target(2)
            x_wedge(3) = x_target(3)
            call chartmap_initial_guess_cartesian_3d(self, x_wedge, u_try)
            call chartmap_refine_cart_seed(self, x_wedge, zeta_period, u_try, &
                                           res_norm, ierr_seed)
            call chartmap_solve_cart(self, x_wedge, u_try(1), u_try(2), u_try(3), &
                                     ierr_try, trace)
            ! Accept on the final residual, not on the solver exit path: the Newton
            ! step-size exit fires once steps reach tol_step in chart units, which
            ! on a reactor-size chart is a Cartesian residual of order
            ! |dx/du|*tol_step, well above tol_res yet far below the documented
            ! inverse tolerance. A target the chart cannot represent keeps a
            ! residual orders of magnitude above the accept tolerance and stays
            ! rejected. The threshold carries the solver's own convergence floors
            ! (tol_newton and the scale-aware term) so a root the Newton
            ! legitimately declares converged is never rejected here.
            call chartmap_eval_cart(self, u_try, x_round)
            res_norm = sqrt(sum((x_wedge - x_round)**2))
            accept_tol = max(chartmap_invert_accept_tol, self%tol_newton, &
                             100.0_dp*epsilon(1.0_dp)*sqrt(sum(x_wedge**2)))
            if (res_norm < accept_tol) ierr = chartmap_from_cyl_ok
            if (ierr /= chartmap_from_cyl_ok) return
            u = u_try
            u(1) = min(max(u(1), 0.0_dp), 1.0_dp)
            u(2) = modulo(u(2), TWOPI)
            u(3) = modulo(u(3) + phi_shift, TWOPI)
            return
        end if

        phi = atan2(x_target(2), x_target(1))
        zeta_seed = [modulo(phi, zeta_period), 0.0_dp, 0.25_dp*zeta_period, &
                     0.50_dp*zeta_period, 0.75_dp*zeta_period]

        trace = chartmap_trace_once_enabled()

        do k = 1, size(zeta_seed)
            call chartmap_try_cart_seed(self, x_target, zeta_seed(k), u_try, res_norm, &
                                        ierr_try, trace .and. k == 1)
            if (ierr_try /= chartmap_from_cyl_ok) cycle
            if (res_norm < best_res) then
                best_res = res_norm
                u = u_try
                ierr = chartmap_from_cyl_ok
            end if
        end do

        if (ierr /= chartmap_from_cyl_ok) return
        u(1) = min(max(u(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(u(2), TWOPI)
        u(3) = modulo(u(3), zeta_period)
    end subroutine chartmap_invert_cart

    subroutine chartmap_try_cart_seed(self, x_target, zeta_seed, u, res_norm, ierr, &
                                      trace)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: zeta_seed
        real(dp), intent(out) :: u(3)
        real(dp), intent(out) :: res_norm
        integer, intent(out) :: ierr
        logical, intent(in) :: trace

        real(dp) :: rho, theta, zeta
        real(dp) :: x_round(3)

        zeta = zeta_seed
        call chartmap_initial_guess(self, x_target, zeta, rho, theta, ierr)
        if (ierr /= chartmap_from_cyl_ok) then
            res_norm = huge(1.0_dp)
            u = 0.0_dp
            return
        end if

        call chartmap_solve_cart(self, x_target, rho, theta, zeta, ierr, trace)
        u = [rho, theta, zeta]
        if (ierr /= chartmap_from_cyl_ok) then
            res_norm = huge(1.0_dp)
            return
        end if

        call chartmap_eval_cart(self, u, x_round)
        res_norm = sqrt(sum((x_target - x_round)**2))
    end subroutine chartmap_try_cart_seed

    subroutine chartmap_solve_cart(self, x_target, rho, theta, zeta, ierr, trace)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(inout) :: rho
        real(dp), intent(inout) :: theta
        real(dp), intent(inout) :: zeta
        integer, intent(out) :: ierr
        logical, intent(in) :: trace

        real(dp) :: delta(3)
        real(dp) :: xw, yw, xw_new, yw_new, zeta_new
        real(dp) :: res_norm, res_norm_new
        real(dp) :: tol_res, tol_step
        real(dp) :: step_norm
        integer :: iter
        logical :: success

        ! Iterate in the pseudo-Cartesian chart (X, Y) = (rho cos theta,
        ! rho sin theta), not polar (rho, theta). The polar Gauss-Newton is singular
        ! at the magnetic axis: the poloidal column dx/dtheta ~ rho, so the
        ! normal-matrix entry |dx/dtheta|^2 -> 0 and the step stalls (err_max_iter)
        ! for any orbit whose particle reaches rho -> 0. (X, Y) carries the same
        ! information with a Jacobian regular through rho = 0 (the standard near-axis
        ! chart, cf. flux_pseudocartesian and Pfefferle et al. arXiv:1412.5464).
        ! zeta iterates unwrapped: the lab-frame position of the co-rotating cart
        ! chart is only 2pi-periodic in zeta, so wrapping the trial into the stored
        ! wedge would teleport it by a field-period rotation. Callers reduce the
        ! solved zeta to their convention. The residual tolerance carries a
        ! scale-aware floor so reactor-size charts (|x| ~ 1e3 cm) are not asked for
        ! a residual below floating-point resolution.
        ierr = chartmap_from_cyl_ok
        tol_res = max(1.0e-10_dp, self%tol_newton, &
                      100.0_dp*epsilon(1.0_dp)*sqrt(sum(x_target**2)))
        tol_step = max(1.0e-10_dp, self%tol_newton)
        xw = rho*cos(theta)
        yw = rho*sin(theta)
        do iter = 1, 60
            call chartmap_newton_delta_xy(self, x_target, xw, yw, zeta, delta, &
                                          res_norm, ierr)
            if (ierr /= chartmap_from_cyl_ok) exit
            if (res_norm < tol_res) exit
            step_norm = sqrt(sum(delta**2))
            if (step_norm < tol_step) then
                ierr = chartmap_from_cyl_err_max_iter
                exit
            end if
            call chartmap_line_search_xy(self, x_target, xw, yw, zeta, delta, &
                                         res_norm, xw_new, yw_new, zeta_new, &
                                         res_norm_new, success)
            if (.not. success) then
                ierr = chartmap_from_cyl_err_max_iter
                exit
            end if
            if (trace) print *, "chartmap_from_cart iter=", iter, " res=", res_norm, &
                " step=", step_norm
            xw = xw_new
            yw = yw_new
            zeta = zeta_new
            if (iter == 60) ierr = chartmap_from_cyl_err_max_iter
        end do
        rho = sqrt(xw**2 + yw**2)
        theta = modulo(atan2(yw, xw), TWOPI)
    end subroutine chartmap_solve_cart

    ! Gauss-Newton step in the pseudo-Cartesian chart (X, Y, zeta). The Cartesian
    ! partials dx/dX, dx/dY come from the polar partials by the chain rule of
    ! X = rho cos theta, Y = rho sin theta: dx/dX = cos theta dx/drho
    ! - (sin theta / rho) dx/dtheta, dx/dY = sin theta dx/drho
    ! + (cos theta / rho) dx/dtheta. dx/dtheta ~ rho cancels the 1/rho, so both
    ! columns stay finite at the axis.
    subroutine chartmap_newton_delta_xy(self, x_target, xw, yw, zeta, delta, &
                                        res_norm, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: xw, yw, zeta
        real(dp), intent(out) :: delta(3)
        real(dp), intent(out) :: res_norm
        integer, intent(out) :: ierr

        real(dp) :: rho, theta, cth, sth, rsafe
        real(dp) :: residual(3), dx_drho(3), dx_dtheta(3), dx_dzeta(3)
        real(dp) :: dx_dxw(3), dx_dyw(3)
        real(dp) :: jtj(3, 3), jtr(3)

        rho = sqrt(xw**2 + yw**2)
        theta = atan2(yw, xw)
        call chartmap_eval_residual_and_partials_cart(self, x_target, rho, theta, &
                                                      zeta, residual, res_norm, &
                                                      dx_drho, dx_dtheta, dx_dzeta)
        ! Take the gauge direction from theta, not from xw/rho: at an exact axis
        ! seed (xw = yw = 0) the quotient degenerates to cth = sth = 0, which zeroes
        ! both pseudo-Cartesian columns and freezes the iteration on the axis.
        ! atan2(0, 0) = 0 gives a well-defined gauge there. The radial floor caps
        ! the 1/rho amplification of axis-level spline noise in dx/dtheta.
        rsafe = max(rho, 1.0e-8_dp)
        cth = cos(theta)
        sth = sin(theta)
        dx_dxw = cth*dx_drho - (sth/rsafe)*dx_dtheta
        dx_dyw = sth*dx_drho + (cth/rsafe)*dx_dtheta

        jtj(1, 1) = dot_product(dx_dxw, dx_dxw)
        jtj(1, 2) = dot_product(dx_dxw, dx_dyw)
        jtj(1, 3) = dot_product(dx_dxw, dx_dzeta)
        jtj(2, 1) = jtj(1, 2)
        jtj(2, 2) = dot_product(dx_dyw, dx_dyw)
        jtj(2, 3) = dot_product(dx_dyw, dx_dzeta)
        jtj(3, 1) = jtj(1, 3)
        jtj(3, 2) = jtj(2, 3)
        jtj(3, 3) = dot_product(dx_dzeta, dx_dzeta)

        jtr(1) = dot_product(dx_dxw, residual)
        jtr(2) = dot_product(dx_dyw, residual)
        jtr(3) = dot_product(dx_dzeta, residual)

        call chartmap_solve_normal_eq3(jtj, jtr, delta, ierr)
    end subroutine chartmap_newton_delta_xy

    ! Backtracking line search in (X, Y, zeta). rho = sqrt(X^2 + Y^2) >= 0 is
    ! automatic, so only the outer bound needs a guard. zeta moves unwrapped: the
    ! lab-frame position is only 2pi-periodic in zeta, so a wedge wrap would
    ! teleport the trial by a field-period rotation and stall the search.
    subroutine chartmap_line_search_xy(self, x_target, xw, yw, zeta, delta, &
                                       res_norm, xw_new, yw_new, zeta_new, &
                                       res_norm_new, success)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: xw, yw, zeta
        real(dp), intent(in) :: delta(3), res_norm
        real(dp), intent(out) :: xw_new, yw_new, zeta_new, res_norm_new
        logical, intent(out) :: success

        real(dp) :: alpha, rho_t, theta_t, uvec(3), x_trial(3)
        integer :: k

        success = .false.
        xw_new = xw; yw_new = yw; zeta_new = zeta; res_norm_new = res_norm
        alpha = 1.0_dp
        do k = 1, 12
            xw_new = xw + alpha*delta(1)
            yw_new = yw + alpha*delta(2)
            rho_t = sqrt(xw_new**2 + yw_new**2)
            if (rho_t > 1.02_dp) then
                alpha = 0.5_dp*alpha
                cycle
            end if
            theta_t = atan2(yw_new, xw_new)
            zeta_new = zeta + alpha*delta(3)
            uvec = [rho_t, theta_t, zeta_new]
            call chartmap_eval_cart(self, uvec, x_trial)
            res_norm_new = sqrt(sum((x_target - x_trial)**2))
            if (res_norm_new < res_norm) then
                success = .true.
                return
            end if
            alpha = 0.5_dp*alpha
        end do
    end subroutine chartmap_line_search_xy

    subroutine chartmap_eval_residual_and_partials_cart(self, x_target, rho, &
                                                        theta, zeta, &
                                                        residual, res_norm, dx_drho, &
                                                        dx_dtheta, dx_dzeta)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: residual(3)
        real(dp), intent(out) :: res_norm
        real(dp), intent(out) :: dx_drho(3)
        real(dp), intent(out) :: dx_dtheta(3)
        real(dp), intent(out) :: dx_dzeta(3)

        real(dp) :: uvec(3)
        real(dp) :: vals(3)
        real(dp) :: dvals(3, 3)

        uvec = [rho, theta, zeta]
        call chartmap_eval_cart_der(self, uvec, vals, dvals)

        residual = x_target - vals
        res_norm = sqrt(sum(residual**2))
        dx_drho = dvals(:, 1)
        dx_dtheta = dvals(:, 2)
        dx_dzeta = dvals(:, 3)
    end subroutine chartmap_eval_residual_and_partials_cart

    subroutine chartmap_solve_normal_eq3(jtj, jtr, delta, ierr)
        real(dp), intent(in) :: jtj(3, 3)
        real(dp), intent(in) :: jtr(3)
        real(dp), intent(out) :: delta(3)
        integer, intent(out) :: ierr

        real(dp) :: a(3, 3), b(3)
        real(dp) :: lambda
        integer :: k
        logical :: ok

        ierr = chartmap_from_cyl_ok
        delta = 0.0_dp

        lambda = 0.0_dp
        do k = 1, 10
            a = jtj
            a(1, 1) = a(1, 1) + lambda
            a(2, 2) = a(2, 2) + lambda
            a(3, 3) = a(3, 3) + lambda
            b = jtr

            call chartmap_solve_3x3(a, b, delta, ok)
            if (ok) return

            if (lambda == 0.0_dp) then
                lambda = 1.0e-12_dp*max(1.0_dp, jtj(1, 1) + jtj(2, 2) + jtj(3, 3))
            else
                lambda = 10.0_dp*lambda
            end if
        end do

        ierr = chartmap_from_cyl_err_singular
    end subroutine chartmap_solve_normal_eq3

    subroutine chartmap_solve_3x3(a, b, x, ok)
        real(dp), intent(inout) :: a(3, 3)
        real(dp), intent(inout) :: b(3)
        real(dp), intent(out) :: x(3)
        logical, intent(out) :: ok
        real(dp), parameter :: pivtol = 1.0e-18_dp
        real(dp) :: pivot
        real(dp) :: factor
        real(dp) :: tmp_row(3)
        real(dp) :: tmp_b
        integer :: i, k, piv

        ok = .true.; x = 0.0_dp
        do i = 1, 3
            piv = i
            pivot = abs(a(i, i))
            do k = i + 1, 3
                if (abs(a(k, i)) > pivot) then
                    piv = k
                    pivot = abs(a(k, i))
                end if
            end do
            if (pivot < pivtol) then
                ok = .false.
                return
            end if
            if (piv /= i) then
                tmp_row = a(i, :)
                a(i, :) = a(piv, :)
                a(piv, :) = tmp_row
                tmp_b = b(i)
                b(i) = b(piv)
                b(piv) = tmp_b
            end if
            pivot = a(i, i)
            a(i, :) = a(i, :)/pivot
            b(i) = b(i)/pivot
            do k = 1, 3
                if (k == i) cycle
                factor = a(k, i)
                a(k, :) = a(k, :) - factor*a(i, :)
                b(k) = b(k) - factor*b(i)
            end do
        end do
        x = b
    end subroutine chartmap_solve_3x3

    subroutine chartmap_initial_guess(self, x_target, zeta, rho, theta, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: rho
        real(dp), intent(out) :: theta
        integer, intent(out) :: ierr

        real(dp) :: uvec(3)
        real(dp) :: axis_x(3), bound_x(3)
        real(dp) :: vals(3)
        real(dp) :: denom
        real(dp) :: e_r(2)
        real(dp) :: r_rel
        real(dp) :: z_rel

        ierr = chartmap_from_cyl_ok

        if (.not. self%has_spl_rz) then
            call chartmap_initial_guess_cartesian(self, x_target, zeta, rho, theta)
            return
        end if

        uvec = [0.0_dp, 0.0_dp, zeta]
        call chartmap_eval_cart(self, uvec, vals)
        axis_x = vals

        e_r(1) = cos(zeta)
        e_r(2) = sin(zeta)
        r_rel = e_r(1)*(x_target(1) - axis_x(1)) + e_r(2)*(x_target(2) - axis_x(2))
        z_rel = x_target(3) - axis_x(3)
        theta = atan2(z_rel, r_rel)

        uvec = [1.0_dp, theta, zeta]
        call chartmap_eval_cart(self, uvec, vals)
        bound_x = vals

        denom = sqrt(sum((bound_x - axis_x)**2))
        if (denom < 1.0e-14_dp) then
            ierr = chartmap_from_cyl_err_invalid
            rho = 0.0_dp
            return
        end if

        rho = sqrt(sum((x_target - axis_x)**2))/denom
        rho = min(max(rho, 0.0_dp), 1.0_dp)
    end subroutine chartmap_initial_guess

    subroutine chartmap_initial_guess_cartesian(self, x_target, zeta, rho, theta)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: rho
        real(dp), intent(out) :: theta

        real(dp) :: uvec(3)
        real(dp) :: x_try(3)
        real(dp) :: best_dist2, dist2
        integer :: nrho_probe, ntheta_probe
        integer :: ir, it

        nrho_probe = min(max(self%nrho, 5), 9)
        ntheta_probe = min(max(self%ntheta, 16), 32)

        best_dist2 = huge(1.0_dp)
        rho = 0.0_dp
        theta = 0.0_dp

        do it = 1, ntheta_probe
            uvec(2) = TWOPI*real(it - 1, dp)/real(ntheta_probe, dp)
            uvec(3) = zeta
            do ir = 1, nrho_probe
                uvec(1) = real(ir - 1, dp)/real(nrho_probe - 1, dp)
                call chartmap_eval_cart(self, uvec, x_try)
                dist2 = sum((x_target - x_try)**2)
                if (dist2 < best_dist2) then
                    best_dist2 = dist2
                    rho = uvec(1)
                    theta = uvec(2)
                end if
            end do
        end do
    end subroutine chartmap_initial_guess_cartesian

    subroutine chartmap_initial_guess_cartesian_3d(self, x_target, u)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(out) :: u(3)

        real(dp) :: uvec(3)
        real(dp) :: x_try(3)
        real(dp) :: best_dist2, dist2
        real(dp) :: zeta_period
        integer :: nrho_probe, ntheta_probe, nzeta_probe
        integer :: ir, it, iz

        nrho_probe = min(max(self%nrho, 9), 17)
        ntheta_probe = min(max(self%ntheta, 32), 64)
        nzeta_probe = min(max(self%nzeta, 24), 48)
        zeta_period = TWOPI/real(self%num_field_periods, dp)

        best_dist2 = huge(1.0_dp)
        u = 0.0_dp

        do iz = 1, nzeta_probe
            uvec(3) = zeta_period*real(iz - 1, dp)/real(nzeta_probe, dp)
            do it = 1, ntheta_probe
                uvec(2) = TWOPI*real(it - 1, dp)/real(ntheta_probe, dp)
                do ir = 1, nrho_probe
                    uvec(1) = real(ir - 1, dp)/real(nrho_probe - 1, dp)
                    call chartmap_eval_cart(self, uvec, x_try)
                    dist2 = sum((x_target - x_try)**2)
                    if (dist2 < best_dist2) then
                        best_dist2 = dist2
                        u = uvec
                    end if
                end do
            end do
        end do
    end subroutine chartmap_initial_guess_cartesian_3d

    subroutine chartmap_refine_cart_seed(self, x_target, zeta_period, u, res_norm, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: zeta_period
        real(dp), intent(inout) :: u(3)
        real(dp), intent(out) :: res_norm
        integer, intent(out) :: ierr

        real(dp) :: trial(3), best_u(3), x_trial(3)
        real(dp) :: step_rho, step_theta, step_zeta
        real(dp) :: best_dist2, dist2
        integer :: iter, ir, it, iz

        ierr = chartmap_from_cyl_ok
        step_rho = 1.0_dp/16.0_dp
        step_theta = TWOPI/64.0_dp
        step_zeta = zeta_period/48.0_dp

        call chartmap_eval_cart(self, u, x_trial)
        best_dist2 = sum((x_target - x_trial)**2)

        do iter = 1, 64
            best_u = u
            do iz = -1, 1
                trial(3) = modulo(u(3) + real(iz, dp)*step_zeta, zeta_period)
                do it = -1, 1
                    trial(2) = modulo(u(2) + real(it, dp)*step_theta, TWOPI)
                    do ir = -1, 1
                        trial(1) = min(max(u(1) + real(ir, dp)*step_rho, 0.0_dp), 1.0_dp)
                        call chartmap_eval_cart(self, trial, x_trial)
                        dist2 = sum((x_target - x_trial)**2)
                        if (dist2 < best_dist2) then
                            best_dist2 = dist2
                            best_u = trial
                        end if
                    end do
                end do
            end do

            if (any(best_u /= u)) then
                u = best_u
            else
                step_rho = 0.5_dp*step_rho
                step_theta = 0.5_dp*step_theta
                step_zeta = 0.5_dp*step_zeta
                if (max(step_rho, max(step_theta, step_zeta)) < 1.0e-12_dp) exit
            end if
        end do

        call chartmap_eval_cart(self, u, x_trial)
        res_norm = sqrt(sum((x_target - x_trial)**2))
        if (res_norm > 1.0e-4_dp) ierr = chartmap_from_cyl_err_max_iter
    end subroutine chartmap_refine_cart_seed

    subroutine chartmap_newton_delta(self, x_target, rho, theta, zeta, delta, &
                                     res_norm, &
                                     ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: delta(2)
        real(dp), intent(out) :: res_norm
        integer, intent(out) :: ierr

        real(dp) :: residual(3)
        real(dp) :: dx_drho(3), dx_dtheta(3)
        real(dp) :: jtj11, jtj12, jtj22, jtr1, jtr2

        ierr = chartmap_from_cyl_ok
        delta = 0.0_dp

        call chartmap_eval_residual_and_partials(self, x_target, rho, theta, zeta, &
                                                 residual, &
                                                 res_norm, dx_drho, dx_dtheta)

        jtj11 = dot_product(dx_drho, dx_drho)
        jtj12 = dot_product(dx_drho, dx_dtheta)
        jtj22 = dot_product(dx_dtheta, dx_dtheta)
        jtr1 = dot_product(dx_drho, residual)
        jtr2 = dot_product(dx_dtheta, residual)

        call chartmap_solve_normal_eq(jtj11, jtj12, jtj22, jtr1, jtr2, delta, ierr)
    end subroutine chartmap_newton_delta

    subroutine chartmap_solve_normal_eq(jtj11, jtj12, jtj22, jtr1, jtr2, delta, ierr)
        real(dp), intent(in) :: jtj11
        real(dp), intent(in) :: jtj12
        real(dp), intent(in) :: jtj22
        real(dp), intent(in) :: jtr1
        real(dp), intent(in) :: jtr2
        real(dp), intent(out) :: delta(2)
        integer, intent(out) :: ierr

        real(dp) :: a11, a12, a22, det, lambda
        integer :: k

        ierr = chartmap_from_cyl_ok
        delta = 0.0_dp

        lambda = 0.0_dp
        do k = 1, 8
            a11 = jtj11 + lambda
            a12 = jtj12
            a22 = jtj22 + lambda
            det = a11*a22 - a12*a12
            if (abs(det) > 1.0e-18_dp) exit
            if (lambda == 0.0_dp) then
                lambda = 1.0e-12_dp*max(1.0_dp, jtj11 + jtj22)
            else
                lambda = 10.0_dp*lambda
            end if
        end do

        if (abs(det) <= 1.0e-18_dp) then
            ierr = chartmap_from_cyl_err_singular
            return
        end if

        delta(1) = (a22*jtr1 - a12*jtr2)/det
        delta(2) = (a11*jtr2 - a12*jtr1)/det
    end subroutine chartmap_solve_normal_eq

    subroutine chartmap_eval_residual_and_partials(self, x_target, rho, theta, zeta, &
                                                   residual, &
                                                   res_norm, dx_drho, dx_dtheta)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: residual(3)
        real(dp), intent(out) :: res_norm
        real(dp), intent(out) :: dx_drho(3)
        real(dp), intent(out) :: dx_dtheta(3)

        real(dp) :: uvec(3)
        real(dp) :: vals(3)
        real(dp) :: dvals(3, 3)

        uvec = [rho, theta, zeta]
        call chartmap_eval_cart_der(self, uvec, vals, dvals)

        residual = x_target - vals
        res_norm = sqrt(sum(residual**2))
        dx_drho = dvals(:, 1)
        dx_dtheta = dvals(:, 2)
    end subroutine chartmap_eval_residual_and_partials

    subroutine chartmap_line_search(self, x_target, rho, theta, zeta, delta, res_norm, &
                                    rho_new, theta_new, res_norm_new, success)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(in) :: delta(2)
        real(dp), intent(in) :: res_norm
        real(dp), intent(out) :: rho_new
        real(dp), intent(out) :: theta_new
        real(dp), intent(out) :: res_norm_new
        logical, intent(out) :: success

        real(dp) :: alpha
        real(dp) :: uvec(3)
        real(dp) :: x_trial(3)
        real(dp) :: residual(3)
        integer :: k

        success = .false.
        rho_new = rho
        theta_new = theta
        res_norm_new = res_norm

        alpha = 1.0_dp
        do k = 1, 12
            rho_new = rho + alpha*delta(1)
            if (rho_new < -0.02_dp .or. rho_new > 1.02_dp) then
                alpha = 0.5_dp*alpha
                cycle
            end if

            theta_new = modulo(theta + alpha*delta(2), TWOPI)
            uvec = [rho_new, theta_new, zeta]
            call chartmap_eval_cart(self, uvec, x_trial)
            residual = x_target - x_trial
            res_norm_new = sqrt(sum(residual**2))

            if (res_norm_new < res_norm) then
                success = .true.
                return
            end if

            alpha = 0.5_dp*alpha
        end do
    end subroutine chartmap_line_search

    subroutine newton_slice(self, x_target, zeta, rho_theta, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3), zeta
        real(dp), intent(out) :: rho_theta(2)
        integer, intent(out) :: ierr

        real(dp) :: rho, theta
        logical :: trace

        ierr = chartmap_from_cyl_ok
        trace = chartmap_trace_once_enabled()
        if (trace) then
            print *, "chartmap_from_cyl trace:"
            print *, "  zeta=", zeta
            print *, "  x_target=", x_target
        end if

        call chartmap_solve_slice(self, x_target, zeta, rho, theta, ierr, trace)

        rho_theta = [rho, theta]
    end subroutine newton_slice

    subroutine chartmap_solve_slice(self, x_target, zeta, rho, theta, ierr, trace)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: rho
        real(dp), intent(out) :: theta
        integer, intent(out) :: ierr
        logical, intent(in) :: trace

        integer :: iter
        real(dp) :: delta(2)
        real(dp) :: res_norm, res_norm_new
        real(dp) :: rho_new, theta_new
        logical :: ok_step
        real(dp), parameter :: tol_res = 1.0e-10_dp, tol_step = 1.0e-10_dp

        ierr = chartmap_from_cyl_ok
        call chartmap_initial_guess(self, x_target, zeta, rho, theta, ierr)
        if (ierr /= chartmap_from_cyl_ok) return

        do iter = 1, 60
            call chartmap_newton_delta(self, x_target, rho, theta, zeta, delta, &
                                       res_norm, ierr)
            if (ierr /= chartmap_from_cyl_ok) exit
            if (res_norm < tol_res) exit
            if (maxval(abs(delta)) < max(tol_step, self%tol_newton)) exit

            call chartmap_line_search(self, x_target, rho, theta, zeta, delta, &
                                      res_norm, rho_new, &
                                      theta_new, res_norm_new, ok_step)
            if (.not. ok_step) then
                if (maxval(abs(delta)) < max(tol_step, self%tol_newton)) exit
                ierr = chartmap_from_cyl_err_max_iter
                exit
            end if

            rho = rho_new
            theta = theta_new

            if (trace) then
                print *, "  iter=", iter, " rho=", rho, " theta=", theta, &
                    " |res|=", res_norm_new, " |d|=", maxval(abs(delta))
            end if
        end do

        if (ierr == chartmap_from_cyl_ok .and. iter >= 60) then
            ierr = chartmap_from_cyl_err_max_iter
        end if
    end subroutine chartmap_solve_slice

end module libneo_coordinates_chartmap
