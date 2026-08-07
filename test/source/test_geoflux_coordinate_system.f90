program test_geoflux_coordinate_system
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_geoflux_coordinate_system, &
                                  geoflux_coordinate_system_t
    use geoflux_coordinates, only: geoflux_to_cyl, cyl_to_geoflux, &
        assign_geoflux_to_cyl_jacobian
    use geoflux_field, only: spline_geoflux_data
    use cylindrical_cartesian, only: cyl_to_cart
    use math_constants, only: pi
    implicit none

    character(len=*), parameter :: fallback_geqdsk = '../../python/tests/test.geqdsk'
    character(len=*), parameter :: env_geqdsk = 'LIBNEO_TEST_GEQDSK'

    character(len=512) :: geqdsk_file
    character(len=512) :: arg_buffer
    integer :: arg_status, arg_len
    integer :: nerrors
    class(coordinate_system_t), allocatable :: cs

    nerrors = 0

    geqdsk_file = fallback_geqdsk
    call get_command_argument(1, value=arg_buffer, length=arg_len, status=arg_status)
    if (arg_status == 0 .and. arg_len > 0 .and. len_trim(arg_buffer) > 0) then
        geqdsk_file = trim(arg_buffer)
    else
        call get_environment_variable(env_geqdsk, value=arg_buffer, &
                                      length=arg_len, status=arg_status)
        if (arg_status == 0 .and. arg_len > 0 .and. len_trim(arg_buffer) > 0) then
            geqdsk_file = trim(arg_buffer)
        end if
    end if

    call spline_geoflux_data(trim(geqdsk_file), 64, 128)

    call make_geoflux_coordinate_system(cs)
    if (.not. allocated(cs)) error stop "Failed to create geoflux coordinate system"

    select type (gcs => cs)
    type is (geoflux_coordinate_system_t)
        call test_evaluate_cyl_matches_module(gcs, nerrors)
        call test_evaluate_cart_consistency(gcs, nerrors)
        call test_coordinate_roundtrip(gcs, nerrors)
        call test_covariant_basis_finite_difference(gcs, nerrors)
        call test_metric_tensor_properties(gcs, nerrors)
        call test_near_axis_jacobian(gcs, nerrors)
        call test_near_axis_volume_derivative(nerrors)
    class default
        error stop "Unexpected coordinate system type"
    end select

    call test_analytic_circular_near_axis_oracle(nerrors)
    call test_analytic_shaped_near_axis_oracle(nerrors)

    if (nerrors > 0) then
        print *, "FAILED:", nerrors, "error(s) in geoflux coordinate system tests"
        error stop 1
    end if
    print *, "All geoflux coordinate system tests passed!"

contains

    subroutine test_evaluate_cyl_matches_module(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), xcyl_oop(3), xcyl_mod(3)
        real(dp) :: tol, err
        integer :: i, n_passed
        real(dp) :: test_s(3), test_theta(3), test_phi(2)
        integer :: is, it, ip

        tol = 1.0e-14_dp
        n_passed = 0

        test_s = [0.1_dp, 0.5_dp, 0.9_dp]
        test_theta = [0.3_dp, 1.5_dp, -0.8_dp]
        test_phi = [0.0_dp, 1.2_dp]

        do ip = 1, 2
            do it = 1, 3
                do is = 1, 3
                    u = [test_s(is), test_theta(it), test_phi(ip)]

                    call cs%evaluate_cyl(u, xcyl_oop)
                    call geoflux_to_cyl(u, xcyl_mod)

                    err = maxval(abs(xcyl_oop - xcyl_mod))
                    if (err < tol) then
                        n_passed = n_passed + 1
                    else
                        print *, "  FAIL: evaluate_cyl mismatch at u=", u
                        print *, "    OOP:", xcyl_oop
                        print *, "    mod:", xcyl_mod
                    end if
                end do
            end do
        end do

        if (n_passed == 18) then
            print *, "  PASS: evaluate_cyl matches geoflux_to_cyl (18 points)"
        else
            print *, "  FAIL: evaluate_cyl mismatch", 18 - n_passed, "/18"
            nerrors = nerrors + 1
        end if
    end subroutine test_evaluate_cyl_matches_module

    subroutine test_evaluate_cart_consistency(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), xcyl(3), xcart_oop(3), xcart_manual(3)
        real(dp) :: tol, err
        integer :: n_passed
        real(dp) :: test_points(3, 4)
        integer :: i

        tol = 1.0e-14_dp
        n_passed = 0

        test_points(:, 1) = [0.2_dp, 0.5_dp, 0.0_dp]
        test_points(:, 2) = [0.6_dp, -1.0_dp, 0.8_dp]
        test_points(:, 3) = [0.4_dp, 2.0_dp, 1.5_dp]
        test_points(:, 4) = [0.8_dp, 0.1_dp, 2.3_dp]

        do i = 1, 4
            u = test_points(:, i)

            call cs%evaluate_cart(u, xcart_oop)
            call cs%evaluate_cyl(u, xcyl)
            call cyl_to_cart(xcyl, xcart_manual)

            err = maxval(abs(xcart_oop - xcart_manual))
            if (err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  FAIL: evaluate_cart inconsistent at u=", u
                print *, "    OOP:", xcart_oop
                print *, "    manual:", xcart_manual
            end if
        end do

        if (n_passed == 4) then
            print *, "  PASS: evaluate_cart = cyl_to_cart(evaluate_cyl)"
        else
            print *, "  FAIL: evaluate_cart consistency", 4 - n_passed, "/4"
            nerrors = nerrors + 1
        end if
    end subroutine test_evaluate_cart_consistency

    subroutine test_coordinate_roundtrip(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), xcyl(3), u_back(3), xcyl_back(3)
        real(dp), parameter :: tol = 1.0e-8_dp
        real(dp) :: err
        integer :: n_passed
        real(dp) :: test_points(3, 5)
        integer :: i

        n_passed = 0

        test_points(:, 1) = [0.2_dp, 0.3_dp, 0.0_dp]
        test_points(:, 2) = [0.5_dp, 1.0_dp, 0.5_dp]
        test_points(:, 3) = [0.7_dp, -0.5_dp, 1.2_dp]
        test_points(:, 4) = [0.3_dp, 2.0_dp, 0.8_dp]
        test_points(:, 5) = [0.9_dp, 0.1_dp, 0.0_dp]

        do i = 1, 5
            u = test_points(:, i)

            call cs%evaluate_cyl(u, xcyl)
            call cyl_to_geoflux(xcyl, u_back)
            call cs%evaluate_cyl(u_back, xcyl_back)

            err = maxval(abs(xcyl_back - xcyl))
            if (err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  FAIL: coordinate roundtrip at u=", u
                print *, "    u_back=", u_back, " err=", err
            end if
        end do

        if (n_passed == 5) then
            print *, "  PASS: u->cyl->u->cyl roundtrip (5 points)"
        else
            print *, "  FAIL: coordinate roundtrip", 5 - n_passed, "/5"
            nerrors = nerrors + 1
        end if
    end subroutine test_coordinate_roundtrip

    subroutine test_covariant_basis_finite_difference(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3), e_cov_fd(3, 3)
        real(dp) :: u_plus(3), u_minus(3), x_plus(3), x_minus(3)
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-4_dp
        real(dp) :: rel_err, max_err
        integer :: k, n_tests, n_passed
        real(dp) :: test_points(3, 4)
        integer :: i

        n_tests = 0
        n_passed = 0
        max_err = 0.0_dp

        test_points(:, 1) = [0.3_dp, 0.5_dp, 0.2_dp]
        test_points(:, 2) = [0.5_dp, 1.2_dp, 0.0_dp]
        test_points(:, 3) = [0.7_dp, -0.3_dp, 1.0_dp]
        test_points(:, 4) = [0.4_dp, 2.5_dp, 0.7_dp]

        do i = 1, 4
            u = test_points(:, i)
            n_tests = n_tests + 1

            call cs%covariant_basis(u, e_cov)

            do k = 1, 3
                u_plus = u
                u_minus = u
                u_plus(k) = u(k) + h
                u_minus(k) = u(k) - h
                call cs%evaluate_cart(u_plus, x_plus)
                call cs%evaluate_cart(u_minus, x_minus)
                e_cov_fd(:, k) = (x_plus - x_minus) / (2.0_dp * h)
            end do

            rel_err = maxval(abs(e_cov - e_cov_fd)) / max(maxval(abs(e_cov_fd)), 1.0e-10_dp)
            if (rel_err > max_err) max_err = rel_err

            if (rel_err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  covariant_basis FD mismatch at u=", u, " rel_err=", rel_err
            end if
        end do

        if (n_passed == n_tests) then
            print *, "  PASS: covariant_basis matches FD (", n_tests, " points, max_err=", max_err, ")"
        else
            print *, "  FAIL: covariant_basis FD check", n_tests - n_passed, "/", n_tests
            nerrors = nerrors + 1
        end if
    end subroutine test_covariant_basis_finite_difference

    subroutine test_metric_tensor_properties(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: identity(3, 3), prod(3, 3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: i, j, k
        logical :: sym_ok, inv_ok, pos_ok
        real(dp) :: test_points(3, 3)
        integer :: ip, n_passed

        n_passed = 0

        test_points(:, 1) = [0.3_dp, 0.8_dp, 0.1_dp]
        test_points(:, 2) = [0.6_dp, 1.5_dp, 0.5_dp]
        test_points(:, 3) = [0.5_dp, -0.5_dp, 1.2_dp]

        do ip = 1, 3
            u = test_points(:, ip)
            call cs%metric_tensor(u, g, ginv, sqrtg)

            sym_ok = .true.
            do i = 1, 3
                do j = i + 1, 3
                    if (abs(g(i, j) - g(j, i)) > tol) sym_ok = .false.
                    if (abs(ginv(i, j) - ginv(j, i)) > tol) sym_ok = .false.
                end do
            end do

            identity = 0.0_dp
            do i = 1, 3
                identity(i, i) = 1.0_dp
            end do

            prod = 0.0_dp
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        prod(i, j) = prod(i, j) + g(i, k) * ginv(k, j)
                    end do
                end do
            end do

            inv_ok = maxval(abs(prod - identity)) < tol
            pos_ok = sqrtg > 0.0_dp

            if (sym_ok .and. inv_ok .and. pos_ok) then
                n_passed = n_passed + 1
            else
                print *, "  metric_tensor fail at u=", u
                if (.not. sym_ok) print *, "    not symmetric"
                if (.not. inv_ok) print *, "    g*ginv /= I"
                if (.not. pos_ok) print *, "    sqrtg <= 0"
            end if
        end do

        if (n_passed == 3) then
            print *, "  PASS: metric_tensor properties verified (3 points)"
        else
            print *, "  FAIL: metric_tensor properties", 3 - n_passed, "/3"
            nerrors = nerrors + 1
        end if
    end subroutine test_metric_tensor_properties

    subroutine test_near_axis_jacobian(cs, nerrors)
        class(geoflux_coordinate_system_t), intent(in) :: cs
        integer, intent(inout) :: nerrors

        real(dp), parameter :: s_values(5) = [ &
            2.44140625e-10_dp, 9.76562500e-10_dp, 3.90625000e-9_dp, &
            1.56250000e-8_dp, 6.25000000e-8_dp]
        real(dp), parameter :: theta = 0.73_dp
        real(dp), parameter :: phi = 0.41_dp
        real(dp), parameter :: max_spread = 0.25_dp
        real(dp), parameter :: max_first_pair_ratio = 1.25_dp
        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: sqrtg_values(size(s_values))
        real(dp) :: spread, first_pair_ratio
        integer :: i

        do i = 1, size(s_values)
            u = [s_values(i), theta, phi]
            call cs%metric_tensor(u, g, ginv, sqrtg)
            sqrtg_values(i) = sqrtg
        end do

        if (.not. all(ieee_is_finite(sqrtg_values)) .or. &
            minval(sqrtg_values) <= 0.0_dp) then
            print *, "  FAIL: near-axis Jacobian is not finite and positive"
            nerrors = nerrors + 1
            return
        end if

        spread = (maxval(sqrtg_values) - minval(sqrtg_values)) / &
            maxval(sqrtg_values)
        first_pair_ratio = max(sqrtg_values(1), sqrtg_values(2)) / &
            min(sqrtg_values(1), sqrtg_values(2))
        if (spread > max_spread .or. first_pair_ratio > max_first_pair_ratio) then
            print *, "  FAIL: near-axis Jacobian does not approach a finite limit"
            print *, "    sqrtg=", sqrtg_values, " spread=", spread, &
                " first_pair_ratio=", first_pair_ratio
            nerrors = nerrors + 1
        else
            print *, "  PASS: near-axis Jacobian has a finite limit"
        end if
    end subroutine test_near_axis_jacobian

    subroutine test_near_axis_volume_derivative(nerrors)
        integer, intent(inout) :: nerrors

        integer, parameter :: ntheta = 1000
        real(dp), parameter :: s_inner = (0.001_dp/64.0_dp)**2
        real(dp), parameter :: s_outer = (0.002_dp/64.0_dp)**2
        real(dp), parameter :: max_relative_change = 1.0e-2_dp
        real(dp) :: dVds_inner, dVds_outer, relative_change

        call integrate_volume_derivative(s_inner, ntheta, dVds_inner)
        call integrate_volume_derivative(s_outer, ntheta, dVds_outer)

        if (.not. all(ieee_is_finite([dVds_inner, dVds_outer]))) then
            print *, "  FAIL: near-axis dV/ds is not finite"
            nerrors = nerrors + 1
            return
        end if
        if (min(dVds_inner, dVds_outer) <= 0.0_dp) then
            print *, "  FAIL: near-axis dV/ds is not positive"
            nerrors = nerrors + 1
            return
        end if

        relative_change = abs(dVds_outer/dVds_inner - 1.0_dp)
        if (relative_change > max_relative_change) then
            print *, "  FAIL: flux-surface-averaged dV/ds has no axis limit"
            print *, "    dV/ds=", dVds_inner, dVds_outer, &
                " relative_change=", relative_change
            nerrors = nerrors + 1
        else
            print *, "  PASS: flux-surface-averaged dV/ds has a finite axis limit"
        end if
    end subroutine test_near_axis_volume_derivative

    subroutine integrate_volume_derivative(s_value, ntheta, dVds_value)
        real(dp), intent(in) :: s_value
        integer, intent(in) :: ntheta
        real(dp), intent(out) :: dVds_value

        real(dp) :: theta, dtheta, xcyl(3), jac(3,3), sqrtg
        integer :: k

        dtheta = 2.0_dp*pi/real(ntheta, dp)
        dVds_value = 0.0_dp
        do k = 1, ntheta
            theta = -pi + real(k, dp)*dtheta
            call geoflux_to_cyl([s_value, theta, 0.0_dp], xcyl)
            call assign_geoflux_to_cyl_jacobian(s_value, theta, jac)
            sqrtg = xcyl(1)*(jac(1,1)*jac(3,2) - jac(1,2)*jac(3,1))
            dVds_value = dVds_value + abs(sqrtg)*dtheta
        end do
        dVds_value = 2.0_dp*pi*dVds_value
    end subroutine integrate_volume_derivative

    subroutine test_analytic_circular_near_axis_oracle(nerrors)
        integer, intent(inout) :: nerrors

        real(dp), parameter :: R0 = 3.0_dp
        real(dp), parameter :: Z0 = -0.4_dp
        real(dp), parameter :: a = 0.7_dp
        real(dp), parameter :: theta = 0.61_dp
        real(dp), parameter :: s_values(7) = [ &
            2.0_dp**(-8), 2.0_dp**(-12), 2.0_dp**(-16), 2.0_dp**(-20), &
            2.0_dp**(-24), 2.0_dp**(-28), 2.0_dp**(-32)]
        real(dp), parameter :: tol = 1.0e-9_dp
        real(dp) :: rho, R, Z, R_s, Z_s, R_theta, Z_theta
        real(dp) :: det, det_expected, s_back, theta_back
        real(dp) :: derivative_error, determinant_error, roundtrip_error
        integer :: i

        derivative_error = 0.0_dp
        determinant_error = 0.0_dp
        roundtrip_error = 0.0_dp

        do i = 1, size(s_values)
            rho = sqrt(s_values(i))
            call circular_surface_oracle(R0, Z0, a, s_values(i), theta, &
                R, Z, R_s, Z_s, R_theta, Z_theta)

            det = R * (R_s * Z_theta - R_theta * Z_s)
            det_expected = R * a**2 / 2.0_dp
            derivative_error = max(derivative_error, &
                abs(2.0_dp * rho * R_s - a * cos(theta)))
            derivative_error = max(derivative_error, &
                abs(2.0_dp * rho * Z_s - a * sin(theta)))
            derivative_error = max(derivative_error, &
                abs(R_theta / rho + a * sin(theta)))
            derivative_error = max(derivative_error, &
                abs(Z_theta / rho - a * cos(theta)))
            determinant_error = max(determinant_error, &
                abs(det - det_expected) / det_expected)

            s_back = ((R - R0)**2 + (Z - Z0)**2) / a**2
            theta_back = atan2(Z - Z0, R - R0)
            roundtrip_error = max(roundtrip_error, &
                abs(s_back - s_values(i)) / s_values(i))
            roundtrip_error = max(roundtrip_error, &
                abs(theta_back - theta))
        end do

        if (derivative_error > tol .or. determinant_error > tol .or. &
            roundtrip_error > tol) then
            print *, "  FAIL: analytic circular near-axis oracle"
            print *, "    derivative=", derivative_error, &
                " determinant=", determinant_error, " roundtrip=", roundtrip_error
            nerrors = nerrors + 1
        else
            print *, "  PASS: analytic circular near-axis oracle"
        end if
    end subroutine test_analytic_circular_near_axis_oracle

    subroutine test_analytic_shaped_near_axis_oracle(nerrors)
        integer, intent(inout) :: nerrors

        real(dp), parameter :: R0 = 3.0_dp
        real(dp), parameter :: Z0 = -0.4_dp
        real(dp), parameter :: a = 0.75_dp
        real(dp), parameter :: b = 0.18_dp
        real(dp), parameter :: theta = 0.61_dp
        real(dp), parameter :: s_values(7) = [ &
            2.0_dp**(-8), 2.0_dp**(-12), 2.0_dp**(-16), 2.0_dp**(-20), &
            2.0_dp**(-24), 2.0_dp**(-28), 2.0_dp**(-32)]
        real(dp), parameter :: tol = 2.0e-9_dp
        real(dp) :: rho, R, Z, R_s, Z_s, R_theta, Z_theta
        real(dp) :: radius, radial_s, r1, r2, r1_theta, r2_theta
        real(dp) :: det, det_expected, rho_back, s_back, theta_back
        real(dp) :: derivative_error, determinant_error, roundtrip_error
        integer :: i

        derivative_error = 0.0_dp
        determinant_error = 0.0_dp
        roundtrip_error = 0.0_dp

        do i = 1, size(s_values)
            rho = sqrt(s_values(i))
            call shaped_surface_oracle(R0, Z0, a, b, s_values(i), theta, &
                R, Z, R_s, Z_s, R_theta, Z_theta, r1, r2, r1_theta, r2_theta)

            radius = hypot(R - R0, Z - Z0)
            radial_s = (r1 + 2.0_dp*r2*rho) / (2.0_dp*rho)
            ! R*(R_s*Z_theta-R_theta*Z_s)=R*r*dr/ds tends to
            ! R0*r1(theta)**2/2, independently of the angular derivatives.
            det = R * (R_s * Z_theta - R_theta * Z_s)
            det_expected = R * radial_s * radius
            derivative_error = max(derivative_error, &
                abs(2.0_dp*rho*R_s - (r1 + 2.0_dp*r2*rho)*cos(theta)))
            derivative_error = max(derivative_error, &
                abs(2.0_dp*rho*Z_s - (r1 + 2.0_dp*r2*rho)*sin(theta)))
            derivative_error = max(derivative_error, &
                abs(R_theta/rho - ((r1_theta + r2_theta*rho)*cos(theta) &
                - (r1 + r2*rho)*sin(theta))))
            derivative_error = max(derivative_error, &
                abs(Z_theta/rho - ((r1_theta + r2_theta*rho)*sin(theta) &
                + (r1 + r2*rho)*cos(theta))))
            determinant_error = max(determinant_error, &
                abs(det - det_expected) / max(abs(det_expected), 1.0e-30_dp))

            ! Invert r=r1*rho+r2*rho**2 with the stable positive-root form.
            rho_back = 2.0_dp*radius / &
                (r1 + sqrt(r1*r1 + 4.0_dp*r2*radius))
            s_back = rho_back*rho_back
            theta_back = atan2(Z - Z0, R - R0)
            roundtrip_error = max(roundtrip_error, &
                abs(s_back - s_values(i)) / s_values(i))
            roundtrip_error = max(roundtrip_error, &
                abs(theta_back - theta))
        end do

        if (derivative_error > tol .or. determinant_error > tol .or. &
            roundtrip_error > tol) then
            print *, "  FAIL: analytic shaped near-axis oracle"
            print *, "    derivative=", derivative_error, &
                " determinant=", determinant_error, " roundtrip=", roundtrip_error
            nerrors = nerrors + 1
        else
            print *, "  PASS: analytic shaped near-axis oracle"
        end if
    end subroutine test_analytic_shaped_near_axis_oracle

    pure subroutine circular_surface_oracle(R0, Z0, a, s, theta, R, Z, R_s, Z_s, &
            R_theta, Z_theta)
        real(dp), intent(in) :: R0, Z0, a, s, theta
        real(dp), intent(out) :: R, Z, R_s, Z_s, R_theta, Z_theta
        real(dp) :: rho

        rho = sqrt(s)
        R = R0 + a * rho * cos(theta)
        Z = Z0 + a * rho * sin(theta)
        R_s = a * cos(theta) / (2.0_dp * rho)
        Z_s = a * sin(theta) / (2.0_dp * rho)
        R_theta = -a * rho * sin(theta)
        Z_theta = a * rho * cos(theta)
    end subroutine circular_surface_oracle

    pure subroutine shaped_surface_oracle(R0, Z0, a, b, s, theta, R, Z, R_s, Z_s, &
            R_theta, Z_theta, r1, r2, r1_theta, r2_theta)
        real(dp), intent(in) :: R0, Z0, a, b, s, theta
        real(dp), intent(out) :: R, Z, R_s, Z_s, R_theta, Z_theta
        real(dp), intent(out) :: r1, r2, r1_theta, r2_theta
        real(dp) :: rho, radius, radius_rho, radius_theta

        rho = sqrt(s)
        r1 = a * (1.0_dp + 0.2_dp*cos(2.0_dp*theta))
        r1_theta = -0.4_dp*a*sin(2.0_dp*theta)
        r2 = b * (0.4_dp + 0.1_dp*sin(3.0_dp*theta))
        r2_theta = 0.3_dp*b*cos(3.0_dp*theta)
        radius = r1*rho + r2*rho*rho
        radius_rho = r1 + 2.0_dp*r2*rho
        radius_theta = r1_theta*rho + r2_theta*rho*rho

        R = R0 + radius*cos(theta)
        Z = Z0 + radius*sin(theta)
        R_s = radius_rho*cos(theta) / (2.0_dp*rho)
        Z_s = radius_rho*sin(theta) / (2.0_dp*rho)
        R_theta = radius_theta*cos(theta) - radius*sin(theta)
        Z_theta = radius_theta*sin(theta) + radius*cos(theta)
    end subroutine shaped_surface_oracle

end program test_geoflux_coordinate_system
