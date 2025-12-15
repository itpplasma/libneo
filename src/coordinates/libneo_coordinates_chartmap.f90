submodule(libneo_coordinates) libneo_coordinates_chartmap
    use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
    use math_constants, only: TWOPI
    use interpolate, only: construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    implicit none

contains

    module subroutine make_chartmap_coordinate_system(cs, filename)
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
        real(dp) :: x_min(3), x_max(3)
        logical :: periodic(3)
        integer :: order(3)

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
        call nc_close(ncid)

        order = [3, 3, 3]
        x_min = [rho(1), theta(1), zeta(1)]
        x_max = [rho(nrho), theta(ntheta), zeta(nzeta)]

        periodic(1) = .false.
        periodic(2) = .true.
        periodic(3) = (nzeta > 1)

        allocate (pos_batch(nrho, ntheta, nzeta, 3))
        pos_batch(:, :, :, 1) = x
        pos_batch(:, :, :, 2) = y
        pos_batch(:, :, :, 3) = z

        call construct_batch_splines_3d(x_min, x_max, pos_batch, order, periodic, &
                                        ccs%spl_xyz)

        ccs%nrho = nrho
        ccs%ntheta = ntheta
        ccs%nzeta = nzeta
    end subroutine initialize_chartmap

    subroutine chartmap_evaluate_point(self, u, x)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: x(3)

        real(dp) :: vals(3)

        call evaluate_batch_splines_3d(self%spl_xyz, u, vals)
        x = vals
    end subroutine chartmap_evaluate_point

    subroutine chartmap_covariant_basis(self, u, e_cov)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: e_cov(3, 3)

        real(dp) :: vals(3)
        real(dp) :: dvals(3, 3)

        call evaluate_batch_splines_3d_der(self%spl_xyz, u, vals, dvals)

        e_cov(1, 1) = dvals(1, 1)
        e_cov(1, 2) = dvals(2, 1)
        e_cov(1, 3) = dvals(3, 1)

        e_cov(2, 1) = dvals(1, 2)
        e_cov(2, 2) = dvals(2, 2)
        e_cov(2, 3) = dvals(3, 2)

        e_cov(3, 1) = dvals(1, 3)
        e_cov(3, 2) = dvals(2, 3)
        e_cov(3, 3) = dvals(3, 3)
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

    subroutine chartmap_from_cyl(self, xcyl, u, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: xcyl(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: x_target(3), rho_theta(2)

        x_target(1) = xcyl(1)*cos(xcyl(2))
        x_target(2) = xcyl(1)*sin(xcyl(2))
        x_target(3) = xcyl(3)

        call newton_slice(self, x_target, xcyl(2), rho_theta, ierr)
        u(1) = min(max(rho_theta(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(rho_theta(2), TWOPI)
        u(3) = modulo(xcyl(2), TWOPI)
    end subroutine chartmap_from_cyl

    subroutine newton_slice(self, x_target, zeta, rho_theta, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3), zeta
        real(dp), intent(out) :: rho_theta(2)
        integer, intent(out) :: ierr

        integer :: iter
        real(dp) :: uvec(3)
        real(dp) :: x_val(3), dx_drho(3), dx_dtheta(3)
        real(dp) :: residual(3), jac(3, 2), jtj(2, 2), jtr(2)
        real(dp) :: det, delta(2)
        real(dp) :: axis_x(3), bound_x(3)
        real(dp) :: rho, theta
        real(dp) :: vals(3)
        real(dp) :: dvals(3, 3)

        ierr = 0

        uvec = [0.0_dp, 0.0_dp, zeta]
        call evaluate_batch_splines_3d(self%spl_xyz, uvec, vals)
        axis_x = vals

        theta = atan2(x_target(2) - axis_x(2), x_target(1) - axis_x(1))

        uvec = [1.0_dp, theta, zeta]
        call evaluate_batch_splines_3d(self%spl_xyz, uvec, vals)
        bound_x = vals

        rho = sqrt(sum((x_target - axis_x)**2)) &
              /max(1.0e-12_dp, sqrt(sum((bound_x - axis_x)**2)))
        rho = min(max(rho, 0.0_dp), 1.0_dp)

        do iter = 1, 30
            uvec = [rho, theta, zeta]

            call evaluate_batch_splines_3d_der(self%spl_xyz, uvec, vals, dvals)

            x_val = vals
            dx_drho(1) = dvals(1, 1)
            dx_drho(2) = dvals(1, 2)
            dx_drho(3) = dvals(1, 3)

            dx_dtheta(1) = dvals(2, 1)
            dx_dtheta(2) = dvals(2, 2)
            dx_dtheta(3) = dvals(2, 3)

            residual = x_target - x_val
            jac(:, 1) = dx_drho
            jac(:, 2) = dx_dtheta

            jtj(1, 1) = dot_product(jac(:, 1), jac(:, 1))
            jtj(1, 2) = dot_product(jac(:, 1), jac(:, 2))
            jtj(2, 1) = jtj(1, 2)
            jtj(2, 2) = dot_product(jac(:, 2), jac(:, 2))

            jtr(1) = dot_product(jac(:, 1), residual)
            jtr(2) = dot_product(jac(:, 2), residual)

            det = jtj(1, 1)*jtj(2, 2) - jtj(1, 2)*jtj(2, 1)
            if (abs(det) < 1.0e-18_dp) then
                ierr = 2
                exit
            end if

            delta(1) = (jtj(2, 2)*jtr(1) - jtj(1, 2)*jtr(2))/det
            delta(2) = (jtj(1, 1)*jtr(2) - jtj(2, 1)*jtr(1))/det

            rho = rho + delta(1)
            theta = modulo(theta + delta(2), TWOPI)

            if (maxval(abs(delta)) < self%tol_newton) exit

            if (rho < -0.05_dp .or. rho > 1.05_dp) then
                ierr = 3
                exit
            end if
        end do

        if (ierr == 0 .and. maxval(abs(delta)) >= self%tol_newton) ierr = 4

        rho_theta = [rho, theta]
    end subroutine newton_slice

end submodule libneo_coordinates_chartmap
