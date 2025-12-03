submodule (libneo_coordinates) libneo_coordinates_gframe
    use gframe_boundary, only: gframe_boundary_t
    use interpolate, only: SplineData3D, construct_splines_3d, &
        evaluate_splines_3d, evaluate_splines_3d_der, destroy_splines_3d
    implicit none

contains

    module subroutine make_gframe_coordinate_system(cs, boundary, nrho)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        type(gframe_boundary_t), intent(in) :: boundary
        integer, intent(in), optional :: nrho

        type(gframe_coordinate_system_t), pointer :: gcs
        integer :: nrho_eff

        if (present(nrho)) then
            nrho_eff = nrho
        else
            nrho_eff = 33
        end if
        allocate(gframe_coordinate_system_t :: cs)

        select type (gcs => cs)
        type is (gframe_coordinate_system_t)
            call initialize_gframe(gcs, boundary, nrho_eff)
        class default
            error stop "make_gframe_coordinate_system: allocation failed"
        end select
    end subroutine make_gframe_coordinate_system

    subroutine initialize_gframe(gcs, boundary, nrho)
        type(gframe_coordinate_system_t), intent(inout) :: gcs
        type(gframe_boundary_t), intent(in) :: boundary
        integer, intent(in) :: nrho

        real(dp), allocatable :: rgrid(:, :, :)
        real(dp), allocatable :: zgrid(:, :, :)
        real(dp), allocatable :: rcos(:, :), rsin(:, :)
        real(dp), allocatable :: zcos(:, :), zsin(:, :)
        real(dp) :: x_min(3), x_max(3)
        logical :: periodic(3)
        integer :: order(3)
        logical :: axisymmetric

        order = [3, 3, 3]
        periodic = [.false., .true., .true.]
        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max = [1.0_dp, two_pi(), two_pi()]

        call build_harmonic_volume(boundary, nrho, rgrid, zgrid, rcos, rsin, &
            zcos, zsin, axisymmetric)

        call construct_splines_3d(x_min, x_max, rgrid, order, periodic, &
            gcs%spl_r)
        call construct_splines_3d(x_min, x_max, zgrid, order, periodic, &
            gcs%spl_z)

        gcs%nrho = size(rgrid, 1)
        gcs%ntheta = size(rgrid, 2)
        gcs%nzeta = size(rgrid, 3)
        gcs%axisymmetric = axisymmetric
        gcs%mmax = size(rcos, 1) - 1

        call move_alloc(rcos, gcs%rcos)
        call move_alloc(rsin, gcs%rsin)
        call move_alloc(zcos, gcs%zcos)
        call move_alloc(zsin, gcs%zsin)

        call destroy_volume(rgrid, zgrid)
    end subroutine initialize_gframe

    subroutine build_harmonic_volume(boundary, nrho, rgrid, zgrid, rcos, rsin, &
        zcos, zsin, axisymmetric)
        type(gframe_boundary_t), intent(in) :: boundary
        integer, intent(in) :: nrho
        real(dp), allocatable, intent(out) :: rgrid(:, :, :)
        real(dp), allocatable, intent(out) :: zgrid(:, :, :)
        real(dp), allocatable, intent(out) :: rcos(:, :)
        real(dp), allocatable, intent(out) :: rsin(:, :)
        real(dp), allocatable, intent(out) :: zcos(:, :)
        real(dp), allocatable, intent(out) :: zsin(:, :)
        logical, intent(out) :: axisymmetric

        integer :: mmax
        integer :: jt, jz, ir, m
        real(dp) :: theta_val, rho_val
        real(dp) :: r_acc, z_acc
        real(dp) :: rb_first, zb_first

        mmax = (boundary%ntheta - 1)/2
        allocate(rcos(mmax + 1, boundary%nzeta))
        allocate(rsin(mmax + 1, boundary%nzeta))
        allocate(zcos(mmax + 1, boundary%nzeta))
        allocate(zsin(mmax + 1, boundary%nzeta))
        call fourier_project(boundary, mmax, rcos, rsin, zcos, zsin)

        allocate(rgrid(nrho, boundary%ntheta, boundary%nzeta))
        allocate(zgrid(nrho, boundary%ntheta, boundary%nzeta))

        axisymmetric = .true.
        do jz = 2, boundary%nzeta
            rb_first = maxval(abs(boundary%Rb(:, jz) - boundary%Rb(:, 1)))
            zb_first = maxval(abs(boundary%Zb(:, jz) - boundary%Zb(:, 1)))
            if (rb_first > 1.0e-12_dp .or. zb_first > 1.0e-12_dp) then
                axisymmetric = .false.
                exit
            end if
        end do

        do jz = 1, boundary%nzeta
            do jt = 1, boundary%ntheta
                theta_val = boundary%theta(jt)
                do ir = 1, nrho
                    rho_val = real(ir - 1, dp) / real(nrho - 1, dp)
                    call reconstruct_point(theta_val, rho_val, jz, mmax, rcos, &
                        rsin, zcos, zsin, r_acc, z_acc)
                    rgrid(ir, jt, jz) = r_acc
                    zgrid(ir, jt, jz) = z_acc
                end do
            end do
        end do
    end subroutine build_harmonic_volume

    subroutine reconstruct_point(theta, rho, iz, mmax, rcos, rsin, zcos, zsin, &
        rout, zout)
        real(dp), intent(in) :: theta, rho
        integer, intent(in) :: iz, mmax
        real(dp), intent(in) :: rcos(:, :), rsin(:, :)
        real(dp), intent(in) :: zcos(:, :), zsin(:, :)
        real(dp), intent(out) :: rout, zout

        integer :: m, m_idx
        real(dp) :: rho_pow

        rout = rcos(1, iz)
        zout = zcos(1, iz)
        rho_pow = rho

        do m = 1, mmax
            m_idx = m + 1
            rout = rout + rho_pow * (rcos(m_idx, iz) * cos(real(m, dp) * theta) &
                + rsin(m_idx, iz) * sin(real(m, dp) * theta))
            zout = zout + rho_pow * (zcos(m_idx, iz) * cos(real(m, dp) * theta) &
                + zsin(m_idx, iz) * sin(real(m, dp) * theta))
            rho_pow = rho_pow * rho
        end do
    end subroutine reconstruct_point

    subroutine fourier_project(boundary, mmax, rcos, rsin, zcos, zsin)
        type(gframe_boundary_t), intent(in) :: boundary
        integer, intent(in) :: mmax
        real(dp), intent(out) :: rcos(:, :), rsin(:, :)
        real(dp), intent(out) :: zcos(:, :), zsin(:, :)

        integer :: jz, jt, m, m_idx, ntheta_eff
        real(dp) :: theta_val, scale

        rcos = 0.0_dp
        rsin = 0.0_dp
        zcos = 0.0_dp
        zsin = 0.0_dp
        ntheta_eff = boundary%ntheta - 1
        scale = 2.0_dp / real(ntheta_eff, dp)

        do jz = 1, boundary%nzeta
            do m = 0, mmax
                m_idx = m + 1
                do jt = 1, ntheta_eff
                    theta_val = boundary%theta(jt)
                    if (m == 0) then
                        rcos(m_idx, jz) = rcos(m_idx, jz) + boundary%Rb(jt, jz)
                        zcos(m_idx, jz) = zcos(m_idx, jz) + boundary%Zb(jt, jz)
                    else
                        rcos(m_idx, jz) = rcos(m_idx, jz) + boundary%Rb(jt, jz) &
                            * cos(real(m, dp) * theta_val)
                        rsin(m_idx, jz) = rsin(m_idx, jz) + boundary%Rb(jt, jz) &
                            * sin(real(m, dp) * theta_val)
                        zcos(m_idx, jz) = zcos(m_idx, jz) + boundary%Zb(jt, jz) &
                            * cos(real(m, dp) * theta_val)
                        zsin(m_idx, jz) = zsin(m_idx, jz) + boundary%Zb(jt, jz) &
                            * sin(real(m, dp) * theta_val)
                    end if
                end do
                if (m == 0) then
                    rcos(m_idx, jz) = rcos(m_idx, jz) / real(ntheta_eff, dp)
                    zcos(m_idx, jz) = zcos(m_idx, jz) / real(ntheta_eff, dp)
                else
                    rcos(m_idx, jz) = scale * rcos(m_idx, jz)
                    rsin(m_idx, jz) = scale * rsin(m_idx, jz)
                    zcos(m_idx, jz) = scale * zcos(m_idx, jz)
                    zsin(m_idx, jz) = scale * zsin(m_idx, jz)
                end if
            end do
        end do
    end subroutine fourier_project

    subroutine evaluate_fourier_map(self, u, R, Z, dR, dZ)
        class(gframe_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: R, Z
        real(dp), intent(out), optional :: dR(3), dZ(3)

        integer :: m
        real(dp) :: rho, theta
        real(dp) :: rho_pow, rho_der
        real(dp) :: cos_m, sin_m

        rho = u(1)
        theta = u(2)

        R = self%rcos(1, 1)
        Z = self%zcos(1, 1)
        if (present(dR)) then
            dR = 0.0_dp
            dZ = 0.0_dp
        end if

        rho_pow = 1.0_dp
        do m = 1, self%mmax
            rho_pow = rho_pow * rho
            cos_m = cos(real(m, dp) * theta)
            sin_m = sin(real(m, dp) * theta)

            R = R + rho_pow * (self%rcos(m + 1, 1) * cos_m + &
                self%rsin(m + 1, 1) * sin_m)
            Z = Z + rho_pow * (self%zcos(m + 1, 1) * cos_m + &
                self%zsin(m + 1, 1) * sin_m)

            if (present(dR)) then
                if (rho == 0.0_dp) then
                    if (m == 1) then
                        rho_der = 1.0_dp
                    else
                        rho_der = 0.0_dp
                    end if
                else
                    rho_der = real(m, dp) * rho_pow / rho
                end if

                dR(1) = dR(1) + rho_der * (self%rcos(m + 1, 1) * cos_m &
                    + self%rsin(m + 1, 1) * sin_m)
                dZ(1) = dZ(1) + rho_der * (self%zcos(m + 1, 1) * cos_m &
                    + self%zsin(m + 1, 1) * sin_m)

                dR(2) = dR(2) + rho_pow * real(m, dp) * (-self%rcos(m + 1, 1) &
                    * sin_m + self%rsin(m + 1, 1) * cos_m)
                dZ(2) = dZ(2) + rho_pow * real(m, dp) * (-self%zcos(m + 1, 1) &
                    * sin_m + self%zsin(m + 1, 1) * cos_m)
            end if
        end do

        if (present(dR)) then
            dR(3) = 0.0_dp
            dZ(3) = 0.0_dp
        end if
    end subroutine evaluate_fourier_map

    subroutine gframe_evaluate_point(self, u, x)
        class(gframe_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: R, Z

        if (self%axisymmetric) then
            call evaluate_fourier_map(self, u, R, Z)
        else
            call evaluate_splines_3d(self%spl_r, u, R)
            call evaluate_splines_3d(self%spl_z, u, Z)
        end if

        x(1) = R
        x(2) = u(3)
        x(3) = Z
    end subroutine gframe_evaluate_point

    subroutine gframe_covariant_basis(self, u, e_cov)
        class(gframe_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        real(dp) :: R, Z
        real(dp) :: dR(3), dZ(3)
        real(dp) :: cos_phi, sin_phi

        if (self%axisymmetric) then
            call evaluate_fourier_map(self, u, R, Z, dR, dZ)
        else
            call evaluate_splines_3d_der(self%spl_r, u, R, dR)
            call evaluate_splines_3d_der(self%spl_z, u, Z, dZ)
        end if

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
    end subroutine gframe_covariant_basis

    subroutine gframe_metric_tensor(self, u, g, ginv, sqrtg)
        class(gframe_coordinate_system_t), intent(in) :: self
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
    end subroutine gframe_metric_tensor

    subroutine gframe_from_cyl(self, xcyl, u, ierr)
        class(gframe_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: rho_theta(2)

        call newton_slice(self, xcyl(1), xcyl(3), xcyl(2), rho_theta, ierr)
        u(1) = min(max(rho_theta(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(rho_theta(2), two_pi())
        u(3) = modulo(xcyl(2), two_pi())
    end subroutine gframe_from_cyl

    subroutine newton_slice(self, R_target, Z_target, zeta, rho_theta, ierr)
        class(gframe_coordinate_system_t), intent(in) :: self
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
        if (self%axisymmetric) then
            call evaluate_fourier_map(self, uvec, axis_R, axis_Z)
        else
            call evaluate_splines_3d(self%spl_r, uvec, axis_R)
            call evaluate_splines_3d(self%spl_z, uvec, axis_Z)
        end if
        theta = atan2(Z_target - axis_Z, R_target - axis_R)
        uvec = [1.0_dp, theta, zeta]
        if (self%axisymmetric) then
            call evaluate_fourier_map(self, uvec, bound_R, bound_Z)
        else
            call evaluate_splines_3d(self%spl_r, uvec, bound_R)
            call evaluate_splines_3d(self%spl_z, uvec, bound_Z)
        end if
        rho = sqrt((R_target - axis_R) ** 2 + (Z_target - axis_Z) ** 2) &
            / max(1.0e-12_dp, sqrt((bound_R - axis_R) ** 2 + (bound_Z - axis_Z) ** 2))
        rho = min(max(rho, 0.0_dp), 1.0_dp)

        do iter = 1, 30
            uvec = [rho, theta, zeta]
            if (self%axisymmetric) then
                call evaluate_fourier_map(self, uvec, R_val, Z_val, dR, dZ)
            else
                call evaluate_splines_3d_der(self%spl_r, uvec, R_val, dR)
                call evaluate_splines_3d_der(self%spl_z, uvec, Z_val, dZ)
            end if

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
            theta = modulo(theta + delta_theta, two_pi())

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

    subroutine destroy_volume(rgrid, zgrid)
        real(dp), allocatable, intent(inout) :: rgrid(:, :, :)
        real(dp), allocatable, intent(inout) :: zgrid(:, :, :)

        if (allocated(rgrid)) deallocate(rgrid)
        if (allocated(zgrid)) deallocate(zgrid)
    end subroutine destroy_volume

    pure real(dp) function two_pi()
        two_pi = 2.0_dp * acos(-1.0_dp)
    end function two_pi

end submodule libneo_coordinates_gframe
