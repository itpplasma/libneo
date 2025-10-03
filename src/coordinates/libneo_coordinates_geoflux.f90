submodule (libneo_coordinates) libneo_coordinates_geoflux
    use geoflux_coordinates, only: geoflux_to_cyl, assign_geoflux_to_cyl_jacobian
    implicit none

    type, extends(coordinate_system_t) :: geoflux_coordinate_system_t
    contains
        procedure :: evaluate_point => geoflux_evaluate_point
        procedure :: covariant_basis => geoflux_covariant_basis
        procedure :: metric_tensor => geoflux_metric_tensor
    end type geoflux_coordinate_system_t

contains

    module subroutine make_geoflux_coordinate_system(cs)
        class(coordinate_system_t), allocatable, intent(out) :: cs
        allocate(geoflux_coordinate_system_t :: cs)
    end subroutine make_geoflux_coordinate_system

    subroutine geoflux_evaluate_point(self, u, x)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: xgeo(3), xcyl(3)

        xgeo(1) = u(1)
        xgeo(2) = u(2)
        xgeo(3) = u(3)

        call geoflux_to_cyl(xgeo, xcyl)

        x(1) = xcyl(1)
        x(2) = xcyl(2)
        x(3) = xcyl(3)
    end subroutine geoflux_evaluate_point

    subroutine geoflux_covariant_basis(self, u, e_cov)
        class(geoflux_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3,3)

        real(dp) :: s, theta, phi
        real(dp) :: R, Z
        real(dp) :: jac_cyl(3,3)
        real(dp) :: cos_phi, sin_phi
        real(dp) :: xgeo(3), xcyl(3)

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

end submodule libneo_coordinates_geoflux
