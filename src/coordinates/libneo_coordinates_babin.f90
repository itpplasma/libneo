submodule (libneo_coordinates) libneo_coordinates_babin
    use babin_boundary, only: babin_boundary_t
    use interpolate, only: SplineData3D, construct_splines_3d, &
        evaluate_splines_3d, evaluate_splines_3d_der, destroy_splines_3d
    implicit none

contains

    module subroutine make_babin_coordinate_system(cs, boundary, nrho)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        type(babin_boundary_t), intent(in) :: boundary
        integer, intent(in), optional :: nrho

        type(babin_coordinate_system_t), pointer :: bcs
        integer :: nrho_eff

        if (present(nrho)) then
            nrho_eff = nrho
        else
            nrho_eff = 33
        end if
        allocate(babin_coordinate_system_t :: cs)

        select type (bcs => cs)
        type is (babin_coordinate_system_t)
            call initialize_babin(bcs, boundary, nrho_eff)
        class default
            error stop "make_babin_coordinate_system: allocation failed"
        end select
    end subroutine make_babin_coordinate_system

    subroutine initialize_babin(bcs, boundary, nrho)
        type(babin_coordinate_system_t), intent(inout) :: bcs
        type(babin_boundary_t), intent(in) :: boundary
        integer, intent(in) :: nrho

        real(dp), allocatable :: rgrid(:, :, :)
        real(dp), allocatable :: zgrid(:, :, :)
        real(dp) :: x_min(3), x_max(3)
        logical :: periodic(3)
        integer :: order(3)

        order = [3, 3, 3]
        periodic(1) = .false.
        periodic(2) = .true.
        periodic(3) = (boundary%nzeta > 1)
        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max = [1.0_dp, two_pi(), two_pi()]

        call build_harmonic_volume(boundary, nrho, rgrid, zgrid)

        call construct_splines_3d(x_min, x_max, rgrid, order, periodic, &
            bcs%spl_r)
        call construct_splines_3d(x_min, x_max, zgrid, order, periodic, &
            bcs%spl_z)

        bcs%nrho = size(rgrid, 1)
        bcs%ntheta = size(rgrid, 2)
        bcs%nzeta = size(rgrid, 3)

        call destroy_volume(rgrid, zgrid)
    end subroutine initialize_babin

    subroutine build_harmonic_volume(boundary, nrho, rgrid, zgrid)
        type(babin_boundary_t), intent(in) :: boundary
        integer, intent(in) :: nrho
        real(dp), allocatable, intent(out) :: rgrid(:, :, :)
        real(dp), allocatable, intent(out) :: zgrid(:, :, :)

        integer :: mmax
        integer :: jt, jz, ir, m
        real(dp) :: theta_val, rho_val
        real(dp) :: r_acc, z_acc
        real(dp), allocatable :: rcos(:, :)
        real(dp), allocatable :: rsin(:, :)
        real(dp), allocatable :: zcos(:, :)
        real(dp), allocatable :: zsin(:, :)

        mmax = (boundary%ntheta - 1)/2
        allocate(rgrid(nrho, boundary%ntheta, boundary%nzeta))
        allocate(zgrid(nrho, boundary%ntheta, boundary%nzeta))

        allocate(rcos(mmax + 1, boundary%nzeta))
        allocate(rsin(mmax + 1, boundary%nzeta))
        allocate(zcos(mmax + 1, boundary%nzeta))
        allocate(zsin(mmax + 1, boundary%nzeta))
        call fourier_project(boundary, mmax, rcos, rsin, zcos, zsin)

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
        deallocate(rcos, rsin, zcos, zsin)
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
        type(babin_boundary_t), intent(in) :: boundary
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
        u(2) = modulo(rho_theta(2), two_pi())
        u(3) = modulo(xcyl(2), two_pi())
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

end submodule libneo_coordinates_babin
