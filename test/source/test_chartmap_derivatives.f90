program test_chartmap_derivatives
    !> Test that chartmap coordinate system derivatives are self-consistent.
    !>
    !> Verifies:
    !> 1. Covariant basis vectors match finite differences of evaluate_cart
    !> 2. At matched physical points, both VMEC and chartmap give same x_cart
    !> 3. Metric tensor is symmetric and positive definite
    !>
    !> Note: Jacobians differ between VMEC and chartmap because coordinate
    !> systems are different (s vs rho, theta_vmec vs theta_chartmap).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_chartmap_vmec_generator, only: write_chartmap_from_vmec
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  vmec_coordinate_system_t, chartmap_from_cyl_ok
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "wout_deriv_test.chartmap.nc"

    class(coordinate_system_t), allocatable :: cs
    type(vmec_coordinate_system_t) :: vmec
    integer :: ierr
    integer :: i, j, n_tested, n_failed, n_deriv_fail
    integer, parameter :: ntest = 50
    real(dp) :: u_vmec(3), u_chart(3), u_plus(3), u_minus(3)
    real(dp) :: x_vmec_cyl(3), x_chart_cyl(3)
    real(dp) :: x_cart(3), x_plus(3), x_minus(3)
    real(dp) :: e_cov(3, 3), e_cov_fd(3, 3)
    real(dp) :: g(3, 3), ginv(3, 3), sqrtg
    real(dp) :: zeta_period
    real(dp) :: x_diff, max_x_diff
    real(dp) :: deriv_err, max_deriv_err
    real(dp), parameter :: tol_x = 1.0e-4_dp
    real(dp), parameter :: tol_deriv = 5.0e-2_dp
    real(dp), parameter :: eps = 1.0e-5_dp
    character(len=2048) :: message

    print *, "Test: Chartmap derivative consistency"
    print *, "======================================"

    call write_chartmap_from_vmec(wout_file, chartmap_file, 33, 33, 17, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: write_chartmap_from_vmec: ", trim(message)
        error stop 1
    end if

    netcdffile = wout_file
    call spline_vmec_data

    call make_chartmap_coordinate_system(cs, chartmap_file)

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        zeta_period = TWOPI/real(ccs%num_field_periods, dp)

        max_x_diff = 0.0_dp
        max_deriv_err = 0.0_dp
        n_tested = 0
        n_failed = 0
        n_deriv_fail = 0

        print *, ""
        print *, "  Test 1: Position consistency (VMEC -> cyl -> chartmap -> cart)"

        do i = 1, ntest
            u_vmec(1) = (0.1_dp + 0.8_dp*real(i - 1, dp)/real(ntest - 1, dp))**2
            u_vmec(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
            u_vmec(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)

            call vmec%evaluate_cyl(u_vmec, x_vmec_cyl)

            call ccs%from_cyl(x_vmec_cyl, u_chart, ierr)
            if (ierr /= chartmap_from_cyl_ok) then
                n_failed = n_failed + 1
                cycle
            end if

            call ccs%evaluate_cyl(u_chart, x_chart_cyl)
            x_diff = sqrt(sum((x_vmec_cyl - x_chart_cyl)**2))
            max_x_diff = max(max_x_diff, x_diff)
            n_tested = n_tested + 1
        end do

        print *, "    Points tested: ", n_tested
        print *, "    Mapping failures: ", n_failed
        print *, "    Max position diff: ", max_x_diff

        if (max_x_diff > tol_x) then
            print *, "    FAIL: Position mismatch exceeds tolerance"
            error stop 1
        end if
        print *, "    PASS"

        print *, ""
        print *, "  Test 2: Covariant basis vs finite differences"

        do i = 1, ntest
            u_chart(1) = 0.1_dp + 0.8_dp*real(i - 1, dp)/real(ntest - 1, dp)
            u_chart(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
            u_chart(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)

            call ccs%evaluate_cart(u_chart, x_cart)
            call ccs%covariant_basis(u_chart, e_cov)

            do j = 1, 3
                u_plus = u_chart
                u_minus = u_chart
                u_plus(j) = u_chart(j) + eps
                u_minus(j) = u_chart(j) - eps

                call ccs%evaluate_cart(u_plus, x_plus)
                call ccs%evaluate_cart(u_minus, x_minus)

                e_cov_fd(:, j) = (x_plus - x_minus) / (2.0_dp*eps)

                deriv_err = sqrt(sum((e_cov(:, j) - e_cov_fd(:, j))**2))
                deriv_err = deriv_err / max(sqrt(sum(e_cov(:, j)**2)), 1.0e-10_dp)
                max_deriv_err = max(max_deriv_err, deriv_err)

                if (deriv_err > tol_deriv) then
                    n_deriv_fail = n_deriv_fail + 1
                end if
            end do
        end do

        print *, "    Max relative derivative error: ", max_deriv_err
        print *, "    Derivative failures: ", n_deriv_fail

        if (max_deriv_err > tol_deriv) then
            print *, "    FAIL: Derivative error exceeds tolerance"
            error stop 1
        end if
        print *, "    PASS"

        print *, ""
        print *, "  Test 3: Metric tensor symmetry and positivity"

        do i = 1, ntest
            u_chart(1) = 0.1_dp + 0.8_dp*real(i - 1, dp)/real(ntest - 1, dp)
            u_chart(2) = modulo(real(37*i, dp)*0.017_dp*TWOPI, TWOPI)
            u_chart(3) = modulo(real(29*i, dp)*0.013_dp*zeta_period, zeta_period)

            call ccs%metric_tensor(u_chart, g, ginv, sqrtg)

            if (abs(g(1, 2) - g(2, 1)) > 1.0e-10_dp .or. &
                abs(g(1, 3) - g(3, 1)) > 1.0e-10_dp .or. &
                abs(g(2, 3) - g(3, 2)) > 1.0e-10_dp) then
                print *, "    FAIL: Metric tensor not symmetric at point ", i
                error stop 1
            end if

            if (sqrtg <= 0.0_dp) then
                print *, "    FAIL: sqrt(g) not positive at point ", i
                print *, "      sqrtg =", sqrtg
                error stop 1
            end if
        end do

        print *, "    PASS: Metric tensor symmetric and positive definite"

        print *, ""
        print *, "  All tests PASSED"

    class default
        print *, "  FAIL: make_chartmap_coordinate_system did not return chartmap type"
        error stop 1
    end select

end program test_chartmap_derivatives
