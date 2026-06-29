program test_vmec_coordinate_system
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: vmec_coordinate_system_t
    use cylindrical_cartesian, only: cyl_to_cart
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    type(vmec_coordinate_system_t) :: vmec
    real(dp) :: zeta_period
    integer :: nerrors

    nerrors = 0

    netcdffile = "wout.nc"
    call spline_vmec_data
    if (nper < 1) then
        print *, "  FAIL: VMEC nper invalid"
        nerrors = nerrors + 1
    end if

    zeta_period = TWOPI/real(max(1, nper), dp)

    call test_basic_evaluate_cyl(vmec, zeta_period, nerrors)
    call test_evaluate_cart_consistency(vmec, zeta_period, nerrors)
    call test_covariant_basis_fd(vmec, zeta_period, nerrors)
    call test_metric_tensor_properties(vmec, zeta_period, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) in VMEC coordinate system test"
        error stop 1
    end if

    print *, "All VMEC coordinate system tests passed!"

contains

    subroutine test_basic_evaluate_cyl(vmec, zeta_period, nerrors)
        type(vmec_coordinate_system_t), intent(in) :: vmec
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), xcyl(3)

        u = [0.25_dp, 0.13_dp*TWOPI, 0.2_dp*zeta_period]
        call vmec%evaluate_cyl(u, xcyl)

        if (.not. (xcyl(1) > 100.0_dp .and. xcyl(1) < 2.0e5_dp)) then
            print *, "  FAIL: VMEC R not in cm range, R=", xcyl(1)
            nerrors = nerrors + 1
        end if
        if (abs(xcyl(2) - u(3)) > 1.0e-14_dp) then
            print *, "  FAIL: VMEC phi mismatch"
            nerrors = nerrors + 1
        end if
        if (.not. (abs(xcyl(3)) < 2.0e5_dp)) then
            print *, "  FAIL: VMEC Z not finite/plausible, Z=", xcyl(3)
            nerrors = nerrors + 1
        else
            print *, "  PASS: evaluate_cyl basic sanity"
        end if
    end subroutine test_basic_evaluate_cyl

    subroutine test_evaluate_cart_consistency(vmec, zeta_period, nerrors)
        type(vmec_coordinate_system_t), intent(in) :: vmec
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), xcyl(3), xcart_oop(3), xcart_manual(3)
        real(dp), parameter :: tol = 1.0e-14_dp
        real(dp) :: err
        integer :: i, n_passed
        real(dp) :: test_points(3, 4)

        n_passed = 0

        test_points(:, 1) = [0.2_dp, 0.5_dp, 0.1_dp*zeta_period]
        test_points(:, 2) = [0.5_dp, 2.0_dp, 0.3_dp*zeta_period]
        test_points(:, 3) = [0.8_dp, 4.5_dp, 0.7_dp*zeta_period]
        test_points(:, 4) = [0.3_dp, 1.0_dp, 0.0_dp]

        do i = 1, 4
            u = test_points(:, i)

            call vmec%evaluate_cart(u, xcart_oop)
            call vmec%evaluate_cyl(u, xcyl)
            call cyl_to_cart(xcyl, xcart_manual)

            err = maxval(abs(xcart_oop - xcart_manual))
            if (err < tol) then
                n_passed = n_passed + 1
            else
                print *, "  FAIL: evaluate_cart inconsistent at u=", u
                print *, "    err=", err
            end if
        end do

        if (n_passed == 4) then
            print *, "  PASS: evaluate_cart = cyl_to_cart(evaluate_cyl)"
        else
            print *, "  FAIL: evaluate_cart consistency", 4 - n_passed, "/4"
            nerrors = nerrors + 1
        end if
    end subroutine test_evaluate_cart_consistency

    subroutine test_covariant_basis_fd(vmec, zeta_period, nerrors)
        type(vmec_coordinate_system_t), intent(in) :: vmec
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), e_cov(3, 3), e_cov_fd(3, 3)
        real(dp) :: u_plus(3), u_minus(3), x_plus(3), x_minus(3)
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol = 1.0e-4_dp
        real(dp) :: rel_err, max_err
        integer :: k, i, n_passed
        real(dp) :: test_points(3, 4)

        n_passed = 0
        max_err = 0.0_dp

        test_points(:, 1) = [0.3_dp, 1.0_dp, 0.2_dp*zeta_period]
        test_points(:, 2) = [0.5_dp, 3.0_dp, 0.5_dp*zeta_period]
        test_points(:, 3) = [0.7_dp, 5.0_dp, 0.8_dp*zeta_period]
        test_points(:, 4) = [0.4_dp, 0.5_dp, 0.0_dp]

        do i = 1, 4
            u = test_points(:, i)

            call vmec%covariant_basis(u, e_cov)

            do k = 1, 3
                u_plus = u
                u_minus = u
                u_plus(k) = u(k) + h
                u_minus(k) = u(k) - h
                call vmec%evaluate_cart(u_plus, x_plus)
                call vmec%evaluate_cart(u_minus, x_minus)
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

        if (n_passed == 4) then
            print *, "  PASS: covariant_basis matches FD (4 points, max_err=", max_err, ")"
        else
            print *, "  FAIL: covariant_basis FD", 4 - n_passed, "/4"
            nerrors = nerrors + 1
        end if
    end subroutine test_covariant_basis_fd

    subroutine test_metric_tensor_properties(vmec, zeta_period, nerrors)
        type(vmec_coordinate_system_t), intent(in) :: vmec
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: identity(3, 3), prod(3, 3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: i, j, k, ip, n_passed
        logical :: sym_ok, inv_ok, pos_ok
        real(dp) :: test_points(3, 3)

        n_passed = 0

        test_points(:, 1) = [0.25_dp, 0.13_dp*TWOPI, 0.2_dp*zeta_period]
        test_points(:, 2) = [0.6_dp, 2.5_dp, 0.4_dp*zeta_period]
        test_points(:, 3) = [0.8_dp, 4.0_dp, 0.9_dp*zeta_period]

        do ip = 1, 3
            u = test_points(:, ip)
            call vmec%metric_tensor(u, g, ginv, sqrtg)

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
            pos_ok = sqrtg > 0.0_dp .and. sqrtg < huge(1.0_dp)

            if (sym_ok .and. inv_ok .and. pos_ok) then
                n_passed = n_passed + 1
            else
                print *, "  metric_tensor fail at u=", u
                if (.not. sym_ok) print *, "    not symmetric"
                if (.not. inv_ok) print *, "    g*ginv /= I, err=", maxval(abs(prod - identity))
                if (.not. pos_ok) print *, "    sqrtg invalid:", sqrtg
            end if
        end do

        if (n_passed == 3) then
            print *, "  PASS: metric_tensor properties (3 points)"
        else
            print *, "  FAIL: metric_tensor properties", 3 - n_passed, "/3"
            nerrors = nerrors + 1
        end if
    end subroutine test_metric_tensor_properties

end program test_vmec_coordinate_system

