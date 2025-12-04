submodule (libneo_coordinates) libneo_coordinates_babin
    use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
    use math_constants, only: TWOPI
    use interpolate, only: SplineData3D, construct_splines_3d, &
        evaluate_splines_3d, evaluate_splines_3d_der, destroy_splines_3d
    implicit none

contains

    module subroutine make_babin_coordinate_system(cs, filename)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        character(len=*), intent(in) :: filename

        allocate(babin_coordinate_system_t :: cs)

        select type (bcs => cs)
        type is (babin_coordinate_system_t)
            call initialize_babin(bcs, filename)
        class default
            error stop "make_babin_coordinate_system: allocation failed"
        end select
    end subroutine make_babin_coordinate_system

    subroutine initialize_babin(bcs, filename)
        type(babin_coordinate_system_t), intent(inout) :: bcs
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: len_r, len_theta, len_z
        integer :: nrho, ntheta, nzeta
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: rgrid(:, :, :)
        real(dp), allocatable :: zgrid(:, :, :)
        real(dp) :: x_min(3), x_max(3)
        logical :: periodic(3)
        integer :: order(3)

        call nc_open(trim(filename), ncid)

        call nc_inq_dim(ncid, 'rho',  len_r)
        call nc_inq_dim(ncid, 'theta', len_theta)
        call nc_inq_dim(ncid, 'zeta', len_z)

        nrho   = len_r
        ntheta = len_theta
        nzeta  = len_z

        allocate(rho(nrho), theta(ntheta), zeta(nzeta))
        call nc_get(ncid, 'rho',   rho)
        call nc_get(ncid, 'theta', theta)
        call nc_get(ncid, 'zeta',  zeta)

        allocate(rgrid(nrho, ntheta, nzeta))
        allocate(zgrid(nrho, ntheta, nzeta))
        call nc_get(ncid, 'R', rgrid)
        call nc_get(ncid, 'Z', zgrid)
        call nc_close(ncid)

        order = [3, 3, 3]
        x_min = [rho(1), theta(1), zeta(1)]
        x_max = [rho(nrho), theta(ntheta), zeta(nzeta)]

        periodic(1) = .false.
        periodic(2) = .true.
        periodic(3) = (nzeta > 1)

        call construct_splines_3d(x_min, x_max, rgrid, order, periodic, &
            bcs%spl_r)
        call construct_splines_3d(x_min, x_max, zgrid, order, periodic, &
            bcs%spl_z)

        bcs%nrho = nrho
        bcs%ntheta = ntheta
        bcs%nzeta = nzeta
    end subroutine initialize_babin

    subroutine babin_evaluate_point(self, u, x)
        class(babin_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: R, Z

        call evaluate_splines_3d(self%spl_r, u, R)
        call evaluate_splines_3d(self%spl_z, u, Z)

        x(1) = R
        x(2) = u(3)
        x(3) = Z
    end subroutine babin_evaluate_point

    subroutine babin_covariant_basis(self, u, e_cov)
        class(babin_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        real(dp) :: R, Z
        real(dp) :: dR(3), dZ(3)
        real(dp) :: cos_phi, sin_phi

        call evaluate_splines_3d_der(self%spl_r, u, R, dR)
        call evaluate_splines_3d_der(self%spl_z, u, Z, dZ)

        cos_phi = cos(u(3))
        sin_phi = sin(u(3))

        e_cov(1, 1) = dR(1) * cos_phi
        e_cov(2, 1) = dR(1) * sin_phi
        e_cov(3, 1) = dZ(1)

        e_cov(1, 2) = dR(2) * cos_phi
        e_cov(2, 2) = dR(2) * sin_phi
        e_cov(3, 2) = dZ(2)

        e_cov(1, 3) = dR(3) * cos_phi - R * sin_phi
        e_cov(2, 3) = dR(3) * sin_phi + R * cos_phi
        e_cov(3, 3) = dZ(3)
    end subroutine babin_covariant_basis

    subroutine babin_metric_tensor(self, u, g, ginv, sqrtg)
        class(babin_coordinate_system_t), intent(in) :: self
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

        det = g(1, 1) * (g(2, 2) * g(3, 3) - g(2, 3) * g(3, 2)) &
            - g(1, 2) * (g(2, 1) * g(3, 3) - g(2, 3) * g(3, 1)) &
            + g(1, 3) * (g(2, 1) * g(3, 2) - g(2, 2) * g(3, 1))

        sqrtg = sqrt(abs(det))

        ginv(1, 1) = (g(2, 2) * g(3, 3) - g(2, 3) * g(3, 2)) / det
        ginv(1, 2) = (g(1, 3) * g(3, 2) - g(1, 2) * g(3, 3)) / det
        ginv(1, 3) = (g(1, 2) * g(2, 3) - g(1, 3) * g(2, 2)) / det
        ginv(2, 1) = (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)) / det
        ginv(2, 2) = (g(1, 1) * g(3, 3) - g(1, 3) * g(3, 1)) / det
        ginv(2, 3) = (g(1, 3) * g(2, 1) - g(1, 1) * g(2, 3)) / det
        ginv(3, 1) = (g(2, 1) * g(3, 2) - g(2, 2) * g(3, 1)) / det
        ginv(3, 2) = (g(1, 2) * g(3, 1) - g(1, 1) * g(3, 2)) / det
        ginv(3, 3) = (g(1, 1) * g(2, 2) - g(1, 2) * g(2, 1)) / det
    end subroutine babin_metric_tensor

    subroutine babin_from_cyl(self, xcyl, u, ierr)
        class(babin_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: rho_theta(2)

        call newton_slice(self, xcyl(1), xcyl(3), xcyl(2), rho_theta, ierr)
        u(1) = min(max(rho_theta(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(rho_theta(2), TWOPI)
        u(3) = modulo(xcyl(2), TWOPI)
    end subroutine babin_from_cyl

    subroutine newton_slice(self, R_target, Z_target, zeta, rho_theta, ierr)
        class(babin_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: R_target, Z_target, zeta
        real(dp), intent(out) :: rho_theta(2)
        integer, intent(out) :: ierr

        integer :: iter
        real(dp) :: uvec(3)
        real(dp) :: R_val, Z_val
        real(dp) :: dR(3), dZ(3)
        real(dp) :: det
        real(dp) :: delta_rho, delta_theta
        real(dp) :: axis_R, axis_Z
        real(dp) :: bound_R, bound_Z
        real(dp) :: rho, theta

        ierr = 0
        uvec = [0.0_dp, 0.0_dp, zeta]
        call evaluate_splines_3d(self%spl_r, uvec, axis_R)
        call evaluate_splines_3d(self%spl_z, uvec, axis_Z)
        theta = atan2(Z_target - axis_Z, R_target - axis_R)
        uvec = [1.0_dp, theta, zeta]
        call evaluate_splines_3d(self%spl_r, uvec, bound_R)
        call evaluate_splines_3d(self%spl_z, uvec, bound_Z)
        rho = sqrt((R_target - axis_R) ** 2 + (Z_target - axis_Z) ** 2) &
            / max(1.0e-12_dp, sqrt((bound_R - axis_R) ** 2 + (bound_Z - axis_Z) ** 2))
        rho = min(max(rho, 0.0_dp), 1.0_dp)

        do iter = 1, 30
            uvec = [rho, theta, zeta]
            call evaluate_splines_3d_der(self%spl_r, uvec, R_val, dR)
            call evaluate_splines_3d_der(self%spl_z, uvec, Z_val, dZ)

            det = dR(1) * dZ(2) - dR(2) * dZ(1)
            if (abs(det) < 1.0e-18_dp) then
                ierr = 2
                exit
            end if

            delta_rho = (dZ(2) * (R_target - R_val) - dR(2) &
                * (Z_target - Z_val)) / det
            delta_theta = (dR(1) * (Z_target - Z_val) - dZ(1) &
                * (R_target - R_val)) / det

            rho = rho + delta_rho
            theta = modulo(theta + delta_theta, TWOPI)

            if (abs(delta_rho) < self%tol_newton .and. abs(delta_theta) &
                < self%tol_newton) exit

            if (rho < -0.05_dp .or. rho > 1.05_dp) then
                ierr = 3
                exit
            end if
        end do

        if (ierr == 0 .and. (abs(delta_rho) >= self%tol_newton .or. &
            abs(delta_theta) >= self%tol_newton)) ierr = 4

        rho_theta = [rho, theta]
    end subroutine newton_slice

end submodule libneo_coordinates_babin
