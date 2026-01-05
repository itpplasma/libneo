module libneo_coordinates_vmec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates_base, only: coordinate_system_t
    use spline_vmec_sub, only: splint_vmec_data
    use cylindrical_cartesian, only: cyl_to_cart
    use math_constants, only: TWOPI
    implicit none
    private

    public :: vmec_coordinate_system_t, make_vmec_coordinate_system

    type, extends(coordinate_system_t) :: vmec_coordinate_system_t
    contains
        procedure :: evaluate_cart => vmec_evaluate_cart
        procedure :: evaluate_cyl => vmec_evaluate_cyl
        procedure :: covariant_basis => vmec_covariant_basis
        procedure :: metric_tensor => vmec_metric_tensor
        procedure :: from_cyl => vmec_from_cyl
    end type vmec_coordinate_system_t

contains

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
        real(dp), intent(out) :: e_cov(3, 3)

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

end module libneo_coordinates_vmec
