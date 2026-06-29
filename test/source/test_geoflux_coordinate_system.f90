program test_geoflux_coordinate_system
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_geoflux_coordinate_system, &
                                  geoflux_coordinate_system_t
    use geoflux_coordinates, only: geoflux_to_cyl, cyl_to_geoflux
    use geoflux_field, only: spline_geoflux_data
    use cylindrical_cartesian, only: cyl_to_cart
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
    class default
        error stop "Unexpected coordinate system type"
    end select

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

end program test_geoflux_coordinate_system
