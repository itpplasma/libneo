module libneo_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    use libneo_coordinate_conventions
    use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
    use math_constants, only: TWOPI
    use netcdf, only: NF90_BYTE, NF90_CHAR, NF90_DOUBLE, NF90_GLOBAL, NF90_INT, &
                      NF90_INT64, NF90_MAX_VAR_DIMS, NF90_NOERR, NF90_NOWRITE, &
                      NF90_SHORT, &
                      nf90_close, nf90_get_att, nf90_get_var, nf90_inq_dimid, &
                      nf90_inq_varid, nf90_inquire_attribute, &
                      nf90_inquire_dimension, nf90_inquire_variable, &
                      nf90_open, nf90_strerror
    use spline_vmec_sub, only: splint_vmec_data
    use cylindrical_cartesian, only: cyl_to_cart, cart_to_cyl
    use geoflux_coordinates, only: geoflux_to_cyl, assign_geoflux_to_cyl_jacobian
    implicit none

    public :: UNKNOWN
    public :: CYL, VMEC, BOOZER
    public :: RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM

    character(len=*), parameter :: chartmap_netcdf_spec = &
                                   "Chartmap NetCDF conventions:"//new_line('a')// &
                                   "- Dimensions: rho, theta, zeta"//new_line('a')// &
                                   "- Variables (required):"//new_line('a')// &
                                   "  rho(rho), theta(theta)"//new_line('a')// &
                                   "  zeta(zeta)"//new_line('a')// &
                                   "  x(zeta,theta,rho)"//new_line('a')// &
                                   "  y(zeta,theta,rho)"//new_line('a')// &
                                   "  z(zeta,theta,rho)"//new_line('a')// &
                                   "- Optional: num_field_periods (integer >= 1)"// &
                                   new_line('a')// &
                                   "- Required: zeta_convention (global attribute)"// &
                                   new_line('a')// &
                                   "  allowed: cyl, vmec"// &
                                   new_line('a')// &
                                   "- Ranges:"//new_line('a')// &
                                   "  rho in [0,1]"//new_line('a')// &
                                   "  theta in [0,2pi)"//new_line('a')// &
                                   "  zeta in [0,2pi/num_field_periods)"// &
                                   new_line('a')// &
                                   "periodic dims exclude endpoint"//new_line('a')// &
                                   "- Storage order:"//new_line('a')// &
                                   "  file dims (zeta,theta,rho)"//new_line('a')// &
                                   "Fortran reads x(rho,theta,zeta)"//new_line('a')// &
                                   "- Units: x,y,z in cm"

    integer, parameter :: chartmap_from_cyl_ok = 0
    integer, parameter :: chartmap_from_cyl_err_max_iter = 1
    integer, parameter :: chartmap_from_cyl_err_singular = 2
    integer, parameter :: chartmap_from_cyl_err_out_of_bounds = 3
    integer, parameter :: chartmap_from_cyl_err_invalid = 4

    integer, parameter :: refcoords_file_unknown = 0
    integer, parameter :: refcoords_file_chartmap = 1
    integer, parameter :: refcoords_file_vmec_wout = 2

    character(len=*), parameter :: chartmap_from_cyl_ierr_spec = &
                                   "chartmap_from_cyl ierr codes:"//new_line('a')// &
                                   "0 ok"//new_line('a')// &
                                   "1 max iterations or no progress"//new_line('a')// &
                                   "2 singular normal equations"//new_line('a')// &
                                   "3 step out of bounds for rho"//new_line('a')// &
                                   "4 invalid mapping slice"

    type, abstract :: coordinate_system_t
    contains
        procedure(evaluate_cart_if), deferred :: evaluate_cart
        procedure(evaluate_cyl_if), deferred :: evaluate_cyl
        procedure(covariant_basis_if), deferred :: covariant_basis
        procedure(metric_tensor_if), deferred :: metric_tensor
        procedure(from_cyl_if), deferred :: from_cyl
        procedure :: cov_to_cart => coordinate_system_cov_to_cart
        procedure :: ctr_to_cart => coordinate_system_ctr_to_cart
    end type coordinate_system_t

    abstract interface
        subroutine evaluate_cart_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine evaluate_cyl_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine covariant_basis_if(self, u, e_cov)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: e_cov(3, 3)
        end subroutine

        subroutine metric_tensor_if(self, u, g, ginv, sqrtg)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        end subroutine

        subroutine from_cyl_if(self, xcyl, u, ierr)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: xcyl(3)
            real(dp), intent(out) :: u(3)
            integer, intent(out) :: ierr
        end subroutine
    end interface

    type, extends(coordinate_system_t) :: vmec_coordinate_system_t
    contains
        procedure :: evaluate_cart => vmec_evaluate_cart
        procedure :: evaluate_cyl => vmec_evaluate_cyl
        procedure :: covariant_basis => vmec_covariant_basis
        procedure :: metric_tensor => vmec_metric_tensor
        procedure :: from_cyl => vmec_from_cyl
    end type vmec_coordinate_system_t

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
        procedure :: from_cyl => chartmap_from_cyl
        procedure :: from_cart => chartmap_from_cart
    end type chartmap_coordinate_system_t

    type, extends(coordinate_system_t) :: geoflux_coordinate_system_t
    contains
        procedure :: evaluate_cart => geoflux_evaluate_cart
        procedure :: evaluate_cyl => geoflux_evaluate_cyl
        procedure :: covariant_basis => geoflux_covariant_basis
        procedure :: metric_tensor => geoflux_metric_tensor
        procedure :: from_cyl => geoflux_from_cyl
    end type geoflux_coordinate_system_t

    public :: make_vmec_coordinate_system
    public :: make_chartmap_coordinate_system
    public :: make_geoflux_coordinate_system
    public :: validate_chartmap_file
    public :: detect_refcoords_file_type

contains

    subroutine coordinate_system_cov_to_cart(self, u, v_cov, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_cov(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_ctr(3)

        call self%metric_tensor(u, g, ginv, sqrtg)
        call self%covariant_basis(u, e_cov)
        v_ctr = matmul(ginv, v_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_cov_to_cart

    subroutine coordinate_system_ctr_to_cart(self, u, v_ctr, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_ctr(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3)

        call self%covariant_basis(u, e_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_ctr_to_cart

    ! ========================================================================
    ! VMEC coordinate system procedures
    ! ========================================================================

    subroutine make_vmec_coordinate_system(cs)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        allocate(vmec_coordinate_system_t :: cs)
    end subroutine make_vmec_coordinate_system

    subroutine vmec_evaluate_cyl(self, u, x)
        class(vmec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: s, theta, varphi
        real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: dl_ds, dl_dt, dl_dp

        associate(dummy => self)
        end associate

        s = u(1)
        theta = u(2)
        varphi = u(3)

        call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, &
                              dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, &
                              dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

        x(1) = R
        x(2) = varphi
        x(3) = Z
    end subroutine vmec_evaluate_cyl

    subroutine vmec_evaluate_cart(self, u, x)
        class(vmec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xcyl(3)

        call self%evaluate_cyl(u, xcyl)
        call cyl_to_cart(xcyl, x)
    end subroutine vmec_evaluate_cart

    subroutine vmec_covariant_basis(self, u, e_cov)
        class(vmec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3,3)

        real(dp) :: s, theta, varphi
        real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: dl_ds, dl_dt, dl_dp
        real(dp) :: cos_phi, sin_phi

        associate(dummy => self)
        end associate

        s = u(1)
        theta = u(2)
        varphi = u(3)

        call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, &
                              dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, &
                              dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

        cos_phi = cos(varphi)
        sin_phi = sin(varphi)

        e_cov(1, 1) = dR_ds*cos_phi
        e_cov(2, 1) = dR_ds*sin_phi
        e_cov(3, 1) = dZ_ds

        e_cov(1, 2) = dR_dt*cos_phi
        e_cov(2, 2) = dR_dt*sin_phi
        e_cov(3, 2) = dZ_dt

        e_cov(1, 3) = dR_dp*cos_phi - R*sin_phi
        e_cov(2, 3) = dR_dp*sin_phi + R*cos_phi
        e_cov(3, 3) = dZ_dp
    end subroutine vmec_covariant_basis

    subroutine vmec_metric_tensor(self, u, g, ginv, sqrtg)
        class(vmec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg

        real(dp) :: e_cov(3,3)
        real(dp) :: det
        integer :: i, j

        call self%covariant_basis(u, e_cov)

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
            end do
        end do

        det = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
            - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
            + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))

        sqrtg = sqrt(abs(det))

        ginv(1,1) = (g(2,2)*g(3,3) - g(2,3)*g(3,2))/det
        ginv(1,2) = (g(1,3)*g(3,2) - g(1,2)*g(3,3))/det
        ginv(1,3) = (g(1,2)*g(2,3) - g(1,3)*g(2,2))/det
        ginv(2,1) = (g(2,3)*g(3,1) - g(2,1)*g(3,3))/det
        ginv(2,2) = (g(1,1)*g(3,3) - g(1,3)*g(3,1))/det
        ginv(2,3) = (g(1,3)*g(2,1) - g(1,1)*g(2,3))/det
        ginv(3,1) = (g(2,1)*g(3,2) - g(2,2)*g(3,1))/det
        ginv(3,2) = (g(1,2)*g(3,1) - g(1,1)*g(3,2))/det
        ginv(3,3) = (g(1,1)*g(2,2) - g(1,2)*g(2,1))/det
    end subroutine vmec_metric_tensor

    subroutine vmec_from_cyl(self, xcyl, u, ierr)
        class(vmec_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        integer, parameter :: max_iter = 50
        real(dp), parameter :: tol_res = 1.0e-10_dp
        real(dp), parameter :: tol_step = 1.0e-12_dp

        real(dp) :: s, theta, varphi
        real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: R, Z, alam
        real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: dl_ds, dl_dt, dl_dp
        real(dp) :: R_axis, Z_axis, R_bnd, Z_bnd
        real(dp) :: rs, zs, dist_axis, dist_bnd
        real(dp) :: f(2), J(2, 2), det, delta(2)
        real(dp) :: res_norm, res_norm_try, alpha
        real(dp) :: s_try, theta_try
        integer :: iter, k

        associate(dummy => self)
        end associate

        ierr = 0
        varphi = xcyl(2)

        call splint_vmec_data(0.0_dp, 0.0_dp, varphi, A_phi, A_theta, &
                              dA_phi_ds, dA_theta_ds, aiota, R_axis, Z_axis, &
                              alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                              dl_ds, dl_dt, dl_dp)

        rs = xcyl(1) - R_axis
        zs = xcyl(3) - Z_axis
        theta = modulo(atan2(zs, rs), TWOPI)

        call splint_vmec_data(1.0_dp, theta, varphi, A_phi, A_theta, dA_phi_ds, &
                              dA_theta_ds, aiota, R_bnd, Z_bnd, alam, dR_ds, &
                              dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, &
                              dl_dp)

        dist_axis = sqrt(rs**2 + zs**2)
        dist_bnd = sqrt((R_bnd - R_axis)**2 + (Z_bnd - Z_axis)**2)
        if (dist_bnd > 1.0e-12_dp) then
            s = min(1.0_dp, max(0.0_dp, (dist_axis/dist_bnd)**2))
        else
            s = 0.0_dp
        end if

        do iter = 1, max_iter
            call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, &
                                  dA_theta_ds, aiota, R, Z, alam, dR_ds, &
                                  dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, &
                                  dl_dt, dl_dp)

            f(1) = R - xcyl(1)
            f(2) = Z - xcyl(3)
            res_norm = sqrt(f(1)**2 + f(2)**2)
            if (res_norm < tol_res) exit

            J(1, 1) = dR_ds
            J(1, 2) = dR_dt
            J(2, 1) = dZ_ds
            J(2, 2) = dZ_dt
            det = J(1, 1)*J(2, 2) - J(1, 2)*J(2, 1)
            if (abs(det) < 1.0e-14_dp) then
                ierr = 2
                return
            end if

            delta(1) = (-f(1)*J(2, 2) + f(2)*J(1, 2))/det
            delta(2) = (-J(1, 1)*f(2) + J(2, 1)*f(1))/det

            if (maxval(abs(delta)) < tol_step) exit

            alpha = 1.0_dp
            do k = 1, 10
                s_try = min(1.0_dp, max(0.0_dp, s + alpha*delta(1)))
                theta_try = modulo(theta + alpha*delta(2), TWOPI)
                call splint_vmec_data(s_try, theta_try, varphi, A_phi, A_theta, &
                                      dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
                                      dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                      dl_ds, dl_dt, dl_dp)
                res_norm_try = sqrt((R - xcyl(1))**2 + (Z - xcyl(3))**2)
                if (res_norm_try < res_norm) then
                    s = s_try
                    theta = theta_try
                    exit
                end if
                alpha = 0.5_dp*alpha
            end do

            if (k > 10) then
                ierr = 1
                return
            end if
        end do

        if (iter > max_iter) then
            ierr = 1
            return
        end if

        u = [s, theta, varphi]
    end subroutine vmec_from_cyl

    ! ========================================================================
    ! Geoflux coordinate system procedures
    ! ========================================================================

    subroutine make_geoflux_coordinate_system(cs)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        allocate(geoflux_coordinate_system_t :: cs)
    end subroutine make_geoflux_coordinate_system

    subroutine geoflux_evaluate_cyl(self, u, x)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        associate(dummy => self)
        end associate

        call geoflux_to_cyl(u, x)
    end subroutine geoflux_evaluate_cyl

    subroutine geoflux_evaluate_cart(self, u, x)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xcyl(3)

        call self%evaluate_cyl(u, xcyl)
        call cyl_to_cart(xcyl, x)
    end subroutine geoflux_evaluate_cart

    subroutine geoflux_covariant_basis(self, u, e_cov)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3,3)

        real(dp) :: s, theta, phi
        real(dp) :: R, Z
        real(dp) :: jac_cyl(3,3)
        real(dp) :: cos_phi, sin_phi
        real(dp) :: xgeo(3), xcyl(3)

        associate(dummy => self)
        end associate

        s = u(1)
        theta = u(2)
        phi = u(3)

        xgeo = [s, theta, phi]
        call geoflux_to_cyl(xgeo, xcyl)
        R = xcyl(1)
        Z = xcyl(3)

        call assign_geoflux_to_cyl_jacobian(s, theta, phi, R, Z, jac_cyl)

        cos_phi = cos(phi)
        sin_phi = sin(phi)

        e_cov(1, 1) = jac_cyl(1,1) * cos_phi
        e_cov(2, 1) = jac_cyl(1,1) * sin_phi
        e_cov(3, 1) = jac_cyl(3,1)

        e_cov(1, 2) = jac_cyl(1,2) * cos_phi
        e_cov(2, 2) = jac_cyl(1,2) * sin_phi
        e_cov(3, 2) = jac_cyl(3,2)

        e_cov(1, 3) = jac_cyl(1,3) * cos_phi - R * sin_phi
        e_cov(2, 3) = jac_cyl(1,3) * sin_phi + R * cos_phi
        e_cov(3, 3) = jac_cyl(3,3)
    end subroutine geoflux_covariant_basis

    subroutine geoflux_metric_tensor(self, u, g, ginv, sqrtg)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: g(3,3), ginv(3,3), sqrtg

        real(dp) :: e_cov(3,3)
        real(dp) :: det
        integer :: i, j

        call self%covariant_basis(u, e_cov)

        do i = 1, 3
            do j = 1, 3
                g(i, j) = dot_product(e_cov(:, i), e_cov(:, j))
            end do
        end do

        det = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
            - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
            + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))

        sqrtg = sqrt(abs(det))

        ginv(1,1) = (g(2,2)*g(3,3) - g(2,3)*g(3,2))/det
        ginv(1,2) = (g(1,3)*g(3,2) - g(1,2)*g(3,3))/det
        ginv(1,3) = (g(1,2)*g(2,3) - g(1,3)*g(2,2))/det
        ginv(2,1) = (g(2,3)*g(3,1) - g(2,1)*g(3,3))/det
        ginv(2,2) = (g(1,1)*g(3,3) - g(1,3)*g(3,1))/det
        ginv(2,3) = (g(1,3)*g(2,1) - g(1,1)*g(2,3))/det
        ginv(3,1) = (g(2,1)*g(3,2) - g(2,2)*g(3,1))/det
        ginv(3,2) = (g(1,2)*g(3,1) - g(1,1)*g(3,2))/det
        ginv(3,3) = (g(1,1)*g(2,2) - g(1,2)*g(2,1))/det
    end subroutine geoflux_metric_tensor

    subroutine geoflux_from_cyl(self, xcyl, u, ierr)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        associate(dummy => self)
        end associate
        associate(dummy2 => xcyl)
        end associate
        u = 0.0_dp
        ierr = 0

        error stop "geoflux_from_cyl: not implemented"
    end subroutine geoflux_from_cyl

    ! ========================================================================
    ! File detection procedure
    ! ========================================================================

    subroutine detect_refcoords_file_type(filename, file_type, ierr, message)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: file_type
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: ncid
        integer :: status
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_x, var_y, var_z
        integer :: var_rmnc
        integer :: ierr_val
        character(len=2048) :: message_val

        ierr = 0
        message = ""
        file_type = refcoords_file_unknown

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            ierr = 1
            message = "open failed: "//trim(nf90_strerror(status))
            return
        end if

        if (nf90_inq_dimid(ncid, "rho", dim_rho) == NF90_NOERR .and. &
            nf90_inq_dimid(ncid, "theta", dim_theta) == NF90_NOERR .and. &
            nf90_inq_dimid(ncid, "zeta", dim_zeta) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "x", var_x) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "y", var_y) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "z", var_z) == NF90_NOERR) then
            status = nf90_close(ncid)
            file_type = refcoords_file_chartmap

            call validate_chartmap_file(filename, ierr_val, message_val)
            if (ierr_val /= 0) then
                ierr = ierr_val
                message = trim(message_val)
            end if
            return
        end if

        if (nf90_inq_varid(ncid, "rmnc", var_rmnc) == NF90_NOERR) then
            file_type = refcoords_file_vmec_wout
            status = nf90_close(ncid)
            return
        end if

        status = nf90_close(ncid)
    end subroutine detect_refcoords_file_type

    ! ========================================================================
    ! Chartmap validator procedures
    ! ========================================================================

    subroutine validate_chartmap_file(filename, ierr, message)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer, parameter :: ok = 0
        integer, parameter :: err_open = 1
        integer, parameter :: err_missing_dim = 2
        integer, parameter :: err_invalid_dim = 3
        integer, parameter :: err_missing_var = 4
        integer, parameter :: err_bad_var = 5
        integer, parameter :: err_missing_units = 6
        integer, parameter :: err_bad_units = 7
        integer, parameter :: err_bad_grid = 8

        integer :: status
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: len_rho, len_theta, len_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer :: num_field_periods
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp) :: period
        character(len=16) :: zeta_convention

        ierr = ok
        message = ""

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            ierr = err_open
            message = "open failed: "//trim(nf90_strerror(status))
            return
        end if

        do
            call validator_require_dim(ncid, "rho", dim_rho, len_rho, ierr, message)
            if (ierr /= ok) exit
            call validator_require_dim(ncid, "theta", dim_theta, len_theta, ierr, &
                                       message)
            if (ierr /= ok) exit
            call validator_require_dim(ncid, "zeta", dim_zeta, len_zeta, ierr, message)
            if (ierr /= ok) exit

            if (len_rho < 2 .or. len_theta < 2 .or. len_zeta < 2) then
                ierr = err_invalid_dim
                message = "invalid dimension length(s)"
                exit
            end if

            call validator_require_1d_real64(ncid, "rho", dim_rho, var_rho, ierr, &
                                             message)
            if (ierr /= ok) exit
            call validator_require_1d_real64(ncid, "theta", dim_theta, var_theta, &
                                             ierr, message)
            if (ierr /= ok) exit
            call validator_require_1d_real64(ncid, "zeta", dim_zeta, var_zeta, ierr, &
                                             message)
            if (ierr /= ok) exit

            call validator_require_3d_real64(ncid, "x", dim_rho, dim_theta, dim_zeta, &
                                             var_x, ierr, message)
            if (ierr /= ok) exit
            call validator_require_3d_real64(ncid, "y", dim_rho, dim_theta, dim_zeta, &
                                             var_y, ierr, message)
            if (ierr /= ok) exit
            call validator_require_3d_real64(ncid, "z", dim_rho, dim_theta, dim_zeta, &
                                             var_z, ierr, message)
            if (ierr /= ok) exit

            call validator_require_units_cm(ncid, "x", var_x, ierr, message)
            if (ierr /= ok) exit
            call validator_require_units_cm(ncid, "y", var_y, ierr, message)
            if (ierr /= ok) exit
            call validator_require_units_cm(ncid, "z", var_z, ierr, message)
            if (ierr /= ok) exit

            call validator_check_optional_zeta_convention(ncid, zeta_convention, &
                                                          ierr, message)
            if (ierr /= ok) exit

            call validator_read_optional_num_field_periods(ncid, num_field_periods, &
                                                           ierr, message)
            if (ierr /= ok) exit

            allocate (rho(len_rho), theta(len_theta), zeta(len_zeta))
            call validator_read_grid(ncid, var_rho, "rho", rho, ierr, message)
            if (ierr /= ok) exit
            call validator_read_grid(ncid, var_theta, "theta", theta, ierr, message)
            if (ierr /= ok) exit
            call validator_read_grid(ncid, var_zeta, "zeta", zeta, ierr, message)
            if (ierr /= ok) exit

            call validator_check_rho_grid(rho, ierr, message)
            if (ierr /= ok) exit

            call validator_check_periodic_grid("theta", theta, TWOPI, ierr, message)
            if (ierr /= ok) exit

            period = TWOPI/real(num_field_periods, dp)
            call validator_check_periodic_grid("zeta", zeta, period, ierr, message)
            if (ierr /= ok) exit

            if (trim(zeta_convention) == "cyl" .or. trim(zeta_convention) == &
                "vmec") then
                call validator_check_cyl_phi_contract(ncid, var_x, var_y, len_rho, &
                                                      len_theta, zeta, &
                                                      period, ierr, message)
                if (ierr /= ok) exit
            end if

            exit
        end do

        status = nf90_close(ncid)
    end subroutine validate_chartmap_file

    subroutine validator_check_cyl_phi_contract(ncid, var_x, var_y, len_rho, &
                                                len_theta, zeta, period, ierr, &
                                                message)
        integer, intent(in) :: ncid
        integer, intent(in) :: var_x, var_y
        integer, intent(in) :: len_rho, len_theta
        real(dp), intent(in) :: zeta(:)
        real(dp), intent(in) :: period
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer, parameter :: err_bad_grid = 8
        real(dp), parameter :: tol_phi = 1.0e-10_dp
        integer :: ir, it, iz
        integer :: i
        integer, dimension(2) :: ir_list
        integer, dimension(3) :: it_list
        integer, dimension(4) :: iz_list
        real(dp) :: xbuf(1, 1, 1), ybuf(1, 1, 1)
        real(dp) :: phi, dphi, r_xy
        integer :: status

        ierr = 0
        message = ""

        if (size(zeta) < 2) then
            ierr = err_bad_grid
            message = "zeta_convention=cyl requires at least 2 zeta points"
            return
        end if

        ir_list = [max(1, 1 + len_rho/2), len_rho]
        it_list = [1, 1 + len_theta/3, 1 + (2*len_theta)/3]
        iz_list = [1, 2, 1 + size(zeta)/2, size(zeta)]

        do ir = 1, size(ir_list)
            do i = 1, size(iz_list)
                iz = max(1, min(size(zeta), iz_list(i)))
                do it = 1, size(it_list)
                    status = nf90_get_var(ncid, var_x, xbuf, &
                                          start=[ir_list(ir), it_list(it), iz], &
                                          count=[1, 1, 1])
                    if (status /= NF90_NOERR) then
                        ierr = err_bad_grid
                        message = "could not read x for cyl contract check"
                        return
                    end if
                    status = nf90_get_var(ncid, var_y, ybuf, &
                                          start=[ir_list(ir), it_list(it), iz], &
                                          count=[1, 1, 1])
                    if (status /= NF90_NOERR) then
                        ierr = err_bad_grid
                        message = "could not read y for cyl contract check"
                        return
                    end if

                    r_xy = sqrt(xbuf(1, 1, 1)**2 + ybuf(1, 1, 1)**2)
                    if (r_xy < 1.0e-12_dp) cycle

                    phi = atan2(ybuf(1, 1, 1), xbuf(1, 1, 1))
                    dphi = abs(modulo(phi - zeta(iz) + 0.5_dp*period, period) - &
                               0.5_dp*period)

                    if (dphi > tol_phi) then
                        ierr = err_bad_grid
                        write (message, '(a,1x,es12.4)') &
                            "cyl requires atan2(y,x)=zeta; dphi=", dphi
                        return
                    end if
                end do
            end do
        end do
    end subroutine validator_check_cyl_phi_contract

    subroutine validator_check_optional_zeta_convention(ncid, value_out, ierr, &
                                                        message)
        integer, intent(in) :: ncid
        character(len=*), intent(out) :: value_out
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status
        integer :: att_type
        integer :: att_len
        character(len=:), allocatable :: value

        ierr = 0
        message = ""
        value_out = ""

        status = nf90_inquire_attribute(ncid, NF90_GLOBAL, "zeta_convention", &
                                        xtype=att_type, len=att_len)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "missing zeta_convention global attribute"
            return
        end if
        if (att_type /= NF90_CHAR .or. att_len < 1) then
            ierr = 5
            message = "zeta_convention must be a global string attribute"
            return
        end if

        allocate (character(len=att_len) :: value)
        status = nf90_get_att(ncid, NF90_GLOBAL, "zeta_convention", value)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read zeta_convention"
            return
        end if

        select case (trim(value))
        case ("cyl", "vmec")
            value_out = trim(value)
        case default
            ierr = 5
            message = "unsupported zeta_convention (must be cyl or vmec)"
        end select
    end subroutine validator_check_optional_zeta_convention

    subroutine validator_require_dim(ncid, name, dimid, dimlen, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimid
        integer, intent(out) :: dimlen
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status

        ierr = 0
        message = ""

        status = nf90_inq_dimid(ncid, trim(name), dimid)
        if (status /= NF90_NOERR) then
            ierr = 2
            message = "missing dimension "//trim(name)
            return
        end if

        status = nf90_inquire_dimension(ncid, dimid, len=dimlen)
        if (status /= NF90_NOERR) then
            ierr = 3
            message = "invalid dimension "//trim(name)
            return
        end if
    end subroutine validator_require_dim

    subroutine validator_require_1d_real64(ncid, name, dimid, varid, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: dimid
        integer, intent(out) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, xtype, ndims
        integer :: dimids(NF90_MAX_VAR_DIMS)

        ierr = 0
        message = ""

        status = nf90_inq_varid(ncid, trim(name), varid)
        if (status /= NF90_NOERR) then
            ierr = 4
            message = "missing variable "//trim(name)
            return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, &
                                       dimids=dimids)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable "//trim(name)
            return
        end if

        if (xtype /= NF90_DOUBLE .or. ndims /= 1 .or. dimids(1) /= dimid) then
            ierr = 5
            message = "variable "//trim(name)//" must be float64 with correct dim"
            return
        end if
    end subroutine validator_require_1d_real64

    subroutine validator_require_3d_real64(ncid, name, dim_rho, dim_theta, dim_zeta, &
                                           varid, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim_rho, dim_theta, dim_zeta
        integer, intent(out) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, xtype, ndims
        integer :: dimids(NF90_MAX_VAR_DIMS)

        ierr = 0
        message = ""

        status = nf90_inq_varid(ncid, trim(name), varid)
        if (status /= NF90_NOERR) then
            ierr = 4
            message = "missing variable "//trim(name)
            return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, &
                                       dimids=dimids)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable "//trim(name)
            return
        end if

        if (xtype /= NF90_DOUBLE .or. ndims /= 3) then
            ierr = 5
            message = "variable "//trim(name)//" must be float64 with 3 dimensions"
            return
        end if

        if (dimids(1) /= dim_rho .or. dimids(2) /= dim_theta .or. dimids(3) /= &
            dim_zeta) then
            ierr = 5
            message = "variable "//trim(name)// &
                      " has wrong dimension order; expected (rho,theta,zeta)"
            return
        end if
    end subroutine validator_require_3d_real64

    subroutine validator_require_units_cm(ncid, name, varid, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, att_len, att_type
        character(len=:), allocatable :: units

        ierr = 0
        message = ""

        status = nf90_inquire_attribute(ncid, varid, "units", xtype=att_type, &
                                        len=att_len)
        if (status /= NF90_NOERR) then
            ierr = 6
            message = "missing units attribute for "//trim(name)
            return
        end if

        if (att_type /= NF90_CHAR .or. att_len < 1) then
            ierr = 7
            message = "invalid units attribute for "//trim(name)
            return
        end if

        allocate (character(len=att_len) :: units)
        status = nf90_get_att(ncid, varid, "units", units)
        if (status /= NF90_NOERR) then
            ierr = 7
            message = "could not read units attribute for "//trim(name)
            return
        end if

        if (trim(units) /= "cm") then
            ierr = 7
            message = "units for "//trim(name)//" must be cm"
            return
        end if
    end subroutine validator_require_units_cm

    subroutine validator_read_optional_num_field_periods(ncid, num_field_periods, &
                                                         ierr, message)
        integer, intent(in) :: ncid
        integer, intent(out) :: num_field_periods
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, varid, xtype, ndims
        character(len=*), parameter :: var_primary = "num_field_periods"
        character(len=*), parameter :: var_legacy = "nfp"

        ierr = 0
        message = ""
        num_field_periods = 1

        status = nf90_inq_varid(ncid, var_primary, varid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_varid(ncid, var_legacy, varid)
            if (status /= NF90_NOERR) return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable num_field_periods"
            return
        end if

        if (ndims /= 0 .or. .not. validator_is_integer_xtype(xtype)) then
            ierr = 5
            message = "num_field_periods must be a scalar integer variable"
            return
        end if

        status = nf90_get_var(ncid, varid, num_field_periods)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read num_field_periods"
            return
        end if

        if (num_field_periods < 1) then
            ierr = 5
            message = "num_field_periods must be >= 1"
            return
        end if
    end subroutine validator_read_optional_num_field_periods

    subroutine validator_read_grid(ncid, varid, name, x, ierr, message)
        integer, intent(in) :: ncid
        integer, intent(in) :: varid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status

        ierr = 0
        message = ""

        status = nf90_get_var(ncid, varid, x)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read "//trim(name)
        end if
    end subroutine validator_read_grid

    subroutine validator_check_rho_grid(rho, ierr, message)
        real(dp), intent(in) :: rho(:)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp), parameter :: rho_min_max = 0.01_dp

        ierr = 0
        message = ""

        if (.not. validator_is_strictly_increasing(rho)) then
            ierr = 8
            message = "rho grid must be strictly increasing"
            return
        end if

        if (rho(1) < 0.0_dp .or. rho(1) > rho_min_max) then
            ierr = 8
            message = "rho(1) must be in [0, 0.01] (axis or small positive for singularity avoidance)"
            return
        end if

        if (abs(rho(size(rho)) - 1.0_dp) > tol) then
            ierr = 8
            message = "rho must end at 1.0"
            return
        end if
    end subroutine validator_check_rho_grid

    subroutine validator_check_periodic_grid(name, x, period, ierr, message)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: period
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: step_expected

        ierr = 0
        message = ""

        if (size(x) == 1) then
            if (abs(x(1) - 0.0_dp) > tol) then
                ierr = 8
                message = trim(name)//" must start at 0"
            end if
            return
        end if

        if (.not. validator_is_strictly_increasing(x)) then
            ierr = 8
            message = trim(name)//" grid must be strictly increasing"
            return
        end if

        step_expected = period/real(size(x), dp)
        if (abs(x(1) - 0.0_dp) > tol) then
            ierr = 8
            message = trim(name)//" must start at 0"
            return
        end if

        if (.not. validator_is_uniform_step(x, step_expected, tol)) then
            ierr = 8
            message = trim(name)//" must be uniform with endpoint excluded"
            return
        end if

        if (x(size(x)) > (period - step_expected + tol)) then
            ierr = 8
            message = trim(name)//" must exclude the period endpoint"
            return
        end if
    end subroutine validator_check_periodic_grid

    logical function validator_is_integer_xtype(xtype)
        integer, intent(in) :: xtype
        validator_is_integer_xtype = (xtype == NF90_BYTE) .or. &
                                     (xtype == NF90_SHORT) .or. &
                                     (xtype == NF90_INT) .or. (xtype == NF90_INT64)
    end function validator_is_integer_xtype

    logical function validator_is_strictly_increasing(x)
        real(dp), intent(in) :: x(:)
        integer :: i

        validator_is_strictly_increasing = .true.
        do i = 2, size(x)
            if (x(i) <= x(i - 1)) then
                validator_is_strictly_increasing = .false.
                return
            end if
        end do
    end function validator_is_strictly_increasing

    logical function validator_is_uniform_step(x, step, tol)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: step
        real(dp), intent(in) :: tol

        real(dp) :: dx_max

        if (size(x) < 2) then
            validator_is_uniform_step = .true.
            return
        end if

        dx_max = maxval(abs((x(2:) - x(:size(x) - 1)) - step))
        validator_is_uniform_step = (dx_max <= 10.0_dp*tol)
    end function validator_is_uniform_step

    ! ========================================================================
    ! Chartmap coordinate system procedures - include from file
    ! ========================================================================

#include "libneo_coordinates_chartmap_impl.inc"

end module libneo_coordinates
