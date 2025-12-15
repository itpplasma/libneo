submodule (libneo_coordinates) libneo_coordinates_vmec
    use spline_vmec_sub, only: splint_vmec_data
    use cylindrical_cartesian, only: cyl_to_cart
    implicit none

contains

    module subroutine make_vmec_coordinate_system(cs)
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

        error stop "vmec_from_cyl: not implemented"
    end subroutine vmec_from_cyl

end submodule libneo_coordinates_vmec
