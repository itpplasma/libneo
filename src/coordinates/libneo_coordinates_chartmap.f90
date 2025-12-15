submodule(libneo_coordinates) libneo_coordinates_chartmap
    use nctools_module, only: nc_open, nc_close, nc_inq_dim, nc_get
    use math_constants, only: TWOPI
    use interpolate, only: construct_batch_splines_3d, &
                           evaluate_batch_splines_3d, evaluate_batch_splines_3d_der
    use netcdf, only: NF90_CHAR, NF90_GLOBAL, NF90_NOERR, nf90_get_att, nf90_get_var, &
                      nf90_inq_varid, nf90_inquire_attribute
    implicit none

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
        real(dp), allocatable :: theta_spl(:), zeta_spl(:)
        real(dp), allocatable :: x_th(:, :, :), y_th(:, :, :), z_th(:, :, :)
        real(dp), allocatable :: x_spl(:, :, :), y_spl(:, :, :), z_spl(:, :, :)
        real(dp) :: x_min(3), x_max(3)
        real(dp) :: zeta_period
        logical :: periodic(3)
        integer :: order(3)
        integer :: num_field_periods
        integer :: zeta_convention

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
        call nc_close(ncid)

        order = [3, 3, 3]

        periodic(1) = .false.
        periodic(2) = .true.
        periodic(3) = (nzeta > 1)

        ccs%num_field_periods = num_field_periods
        ccs%zeta_convention = zeta_convention
        zeta_period = TWOPI/real(num_field_periods, dp)

        call chartmap_extend_theta(nrho, ntheta, nzeta, theta, x, y, z, theta_spl, &
                                   x_th, y_th, &
                                   z_th)

        if (periodic(3)) then
            call chartmap_extend_zeta(nrho, ntheta + 1, nzeta, zeta, zeta_period, &
                                      x_th, y_th, &
                                      z_th, zeta_spl, x_spl, y_spl, z_spl)
        else
            allocate (zeta_spl(nzeta))
            zeta_spl = zeta
            call move_alloc(x_th, x_spl)
            call move_alloc(y_th, y_spl)
            call move_alloc(z_th, z_spl)
        end if

        x_min = [rho(1), theta_spl(1), zeta_spl(1)]
        x_max = [rho(nrho), theta_spl(size(theta_spl)), zeta_spl(size(zeta_spl))]

        allocate (pos_batch(nrho, size(theta_spl), size(zeta_spl), 3))
        pos_batch(:, :, :, 1) = x_spl
        pos_batch(:, :, :, 2) = y_spl
        pos_batch(:, :, :, 3) = z_spl

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
        real(dp) :: zeta
        real(dp) :: zeta_period

        x_target(1) = xcyl(1)*cos(xcyl(2))
        x_target(2) = xcyl(1)*sin(xcyl(2))
        x_target(3) = xcyl(3)

        zeta_period = TWOPI/real(self%num_field_periods, dp)
        if (self%zeta_convention == CYL) then
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

    subroutine chartmap_invert_cart(self, x_target, u, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: ierr

        real(dp) :: u_try(3)
        real(dp) :: res_norm, best_res
        real(dp) :: zeta_period
        real(dp) :: phi
        real(dp) :: zeta_seed(5)
        logical :: trace
        integer :: k
        integer :: ierr_try

        ierr = chartmap_from_cyl_err_max_iter
        best_res = huge(1.0_dp)
        u = 0.0_dp

        zeta_period = TWOPI/real(self%num_field_periods, dp)
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

        call evaluate_batch_splines_3d(self%spl_xyz, u, x_round)
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
        real(dp) :: rho_new, theta_new, zeta_new
        real(dp) :: res_norm, res_norm_new
        real(dp) :: tol_res, tol_step
        real(dp) :: step_norm
        real(dp) :: zeta_period
        integer :: iter
        logical :: success

        ierr = chartmap_from_cyl_ok
        zeta_period = TWOPI/real(self%num_field_periods, dp)
        zeta = modulo(zeta, zeta_period)
        tol_res = max(1.0e-10_dp, self%tol_newton)
        tol_step = max(1.0e-10_dp, self%tol_newton)
        do iter = 1, 60
            call chartmap_newton_delta_cart(self, x_target, rho, theta, zeta, delta, &
                                            res_norm, ierr)
            if (ierr /= chartmap_from_cyl_ok) return
            if (res_norm < tol_res) return
            step_norm = sqrt(sum(delta**2))
            if (step_norm < tol_step) then
                ierr = chartmap_from_cyl_err_max_iter
                return
            end if
            call chartmap_line_search_cart(self, x_target, rho, theta, zeta, &
                                           zeta_period, &
                                           delta, res_norm, rho_new, theta_new, &
                                           zeta_new, &
                                           res_norm_new, success)
            if (.not. success) then
                ierr = chartmap_from_cyl_err_max_iter
                return
            end if
            if (trace) print *, "chartmap_from_cart iter=", iter, " res=", res_norm, &
                " step=", step_norm
            rho = rho_new
            theta = theta_new
            zeta = zeta_new
        end do
        ierr = chartmap_from_cyl_err_max_iter
    end subroutine chartmap_solve_cart

    subroutine chartmap_newton_delta_cart(self, x_target, rho, theta, zeta, delta, &
                                          res_norm, ierr)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(out) :: delta(3)
        real(dp), intent(out) :: res_norm
        integer, intent(out) :: ierr

        real(dp) :: residual(3)
        real(dp) :: dx_drho(3), dx_dtheta(3), dx_dzeta(3)
        real(dp) :: jtj(3, 3), jtr(3)

        ierr = chartmap_from_cyl_ok
        delta = 0.0_dp

        call chartmap_eval_residual_and_partials_cart(self, x_target, rho, &
                                                      theta, zeta, &
                                                      residual, res_norm, dx_drho, &
                                                      dx_dtheta, dx_dzeta)

        jtj(1, 1) = dot_product(dx_drho, dx_drho)
        jtj(1, 2) = dot_product(dx_drho, dx_dtheta)
        jtj(1, 3) = dot_product(dx_drho, dx_dzeta)
        jtj(2, 1) = jtj(1, 2)
        jtj(2, 2) = dot_product(dx_dtheta, dx_dtheta)
        jtj(2, 3) = dot_product(dx_dtheta, dx_dzeta)
        jtj(3, 1) = jtj(1, 3)
        jtj(3, 2) = jtj(2, 3)
        jtj(3, 3) = dot_product(dx_dzeta, dx_dzeta)

        jtr(1) = dot_product(dx_drho, residual)
        jtr(2) = dot_product(dx_dtheta, residual)
        jtr(3) = dot_product(dx_dzeta, residual)

        call chartmap_solve_normal_eq3(jtj, jtr, delta, ierr)
    end subroutine chartmap_newton_delta_cart

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
        call evaluate_batch_splines_3d_der(self%spl_xyz, uvec, vals, dvals)

        residual = x_target - vals
        res_norm = sqrt(sum(residual**2))
        dx_drho = dvals(1, :)
        dx_dtheta = dvals(2, :)
        dx_dzeta = dvals(3, :)
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

    subroutine chartmap_line_search_cart(self, x_target, rho, theta, zeta, &
                                         zeta_period, &
                                         delta, res_norm, rho_new, theta_new, &
                                         zeta_new, &
                                         res_norm_new, success)
        class(chartmap_coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: x_target(3)
        real(dp), intent(in) :: rho
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: zeta
        real(dp), intent(in) :: zeta_period
        real(dp), intent(in) :: delta(3)
        real(dp), intent(in) :: res_norm
        real(dp), intent(out) :: rho_new
        real(dp), intent(out) :: theta_new
        real(dp), intent(out) :: zeta_new
        real(dp), intent(out) :: res_norm_new
        logical, intent(out) :: success

        real(dp) :: alpha
        real(dp) :: uvec(3)
        real(dp) :: x_trial(3)
        real(dp) :: residual(3)
        integer :: k

        success = .false.; rho_new = rho; theta_new = theta; zeta_new = zeta
        res_norm_new = res_norm
        alpha = 1.0_dp
        do k = 1, 12
            rho_new = rho + alpha*delta(1)
            if (rho_new < -0.02_dp .or. rho_new > 1.02_dp) then
                alpha = 0.5_dp*alpha
                cycle
            end if

            theta_new = modulo(theta + alpha*delta(2), TWOPI)
            zeta_new = modulo(zeta + alpha*delta(3), zeta_period)
            uvec = [rho_new, theta_new, zeta_new]
            call evaluate_batch_splines_3d(self%spl_xyz, uvec, x_trial)
            residual = x_target - x_trial
            res_norm_new = sqrt(sum(residual**2))

            if (res_norm_new < res_norm) then
                success = .true.
                return
            end if

            alpha = 0.5_dp*alpha
        end do
    end subroutine chartmap_line_search_cart

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

        uvec = [0.0_dp, 0.0_dp, zeta]
        call evaluate_batch_splines_3d(self%spl_xyz, uvec, vals)
        axis_x = vals

        e_r(1) = cos(zeta)
        e_r(2) = sin(zeta)
        r_rel = e_r(1)*(x_target(1) - axis_x(1)) + e_r(2)*(x_target(2) - axis_x(2))
        z_rel = x_target(3) - axis_x(3)
        theta = atan2(z_rel, r_rel)

        uvec = [1.0_dp, theta, zeta]
        call evaluate_batch_splines_3d(self%spl_xyz, uvec, vals)
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
        call evaluate_batch_splines_3d_der(self%spl_xyz, uvec, vals, dvals)

        residual = x_target - vals
        res_norm = sqrt(sum(residual**2))
        dx_drho = dvals(1, :)
        dx_dtheta = dvals(2, :)
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
            call evaluate_batch_splines_3d(self%spl_xyz, uvec, x_trial)
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

end submodule libneo_coordinates_chartmap
