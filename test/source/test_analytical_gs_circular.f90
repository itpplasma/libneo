program test_analytical_gs_circular
    use iso_fortran_env, only: dp => real64
    use analytical_gs_circular
    implicit none

    integer :: test_count, pass_count

    test_count = 0
    pass_count = 0

    print *, "Testing analytical_gs_circular module..."
    print *, ""

    call test_basis_functions()
    call test_basis_derivatives()

    print *, ""
    print *, "===================="
    print *, "Test Summary:"
    print *, "  Total tests: ", test_count
    print *, "  Passed:      ", pass_count
    print *, "  Failed:      ", test_count - pass_count
    print *, "===================="

    if (pass_count == test_count) then
        print *, "All tests passed!"
    else
        error stop "Some tests failed!"
    end if

contains

    subroutine test_basis_functions()
        real(dp) :: x, y, expected, result, tol
        integer :: i

        print *, "Testing basis functions..."
        tol = 1.0e-12_dp

        ! Test psi_0 at origin (x=1, y=0)
        x = 1.0_dp
        y = 0.0_dp
        result = psi_0(x, y)
        expected = 1.0_dp
        call assert_close(result, expected, tol, "psi_0 at (1,0)")

        ! Test psi_1 at origin
        result = psi_1(x, y)
        expected = 1.0_dp
        call assert_close(result, expected, tol, "psi_1 at (1,0)")

        ! Test psi_2 at origin
        result = psi_2(x, y)
        expected = 0.0_dp
        call assert_close(result, expected, tol, "psi_2 at (1,0)")

        ! Test psi_3 at (x=1.2, y=0.1)
        x = 1.2_dp
        y = 0.1_dp
        result = psi_3(x, y)
        expected = x**2 * y**2
        call assert_close(result, expected, tol, "psi_3 at (1.2,0.1)")

        ! Test psi_4 at (x=1.2, y=0.1)
        result = psi_4(x, y)
        expected = y**4
        call assert_close(result, expected, tol, "psi_4 at (1.2,0.1)")

        ! Test psi_5 at (x=1.2, y=0.1)
        result = psi_5(x, y)
        expected = x**4 * y**2
        call assert_close(result, expected, tol, "psi_5 at (1.2,0.1)")

        ! Test psi_6 at (x=1.2, y=0.1)
        result = psi_6(x, y)
        expected = x**2 * y**4
        call assert_close(result, expected, tol, "psi_6 at (1.2,0.1)")

        print *, ""
    end subroutine test_basis_functions

    subroutine test_basis_derivatives()
        real(dp) :: x, y, dx, dy
        real(dp) :: dpsi_dx_fd, dpsi_dy_fd
        real(dp) :: dpsi_dx, dpsi_dy
        real(dp) :: tol_fd
        integer :: i

        print *, "Testing basis function derivatives (finite differences)..."

        ! Finite difference step
        dx = 1.0e-8_dp
        dy = 1.0e-8_dp
        tol_fd = 1.0e-6_dp  ! Looser tolerance for finite differences

        ! Test point
        x = 1.2_dp
        y = 0.15_dp

        ! Test dpsi_0/dx
        dpsi_dx_fd = (psi_0(x+dx, y) - psi_0(x-dx, y)) / (2.0_dp * dx)
        call dpsi_0_derivatives(x, y, dpsi_dx, dpsi_dy)
        call assert_close(dpsi_dx, dpsi_dx_fd, tol_fd, "dpsi_0/dx finite diff")

        ! Test dpsi_0/dy
        dpsi_dy_fd = (psi_0(x, y+dy) - psi_0(x, y-dy)) / (2.0_dp * dy)
        call assert_close(dpsi_dy, dpsi_dy_fd, tol_fd, "dpsi_0/dy finite diff")

        ! Test dpsi_1/dx
        dpsi_dx_fd = (psi_1(x+dx, y) - psi_1(x-dx, y)) / (2.0_dp * dx)
        call dpsi_1_derivatives(x, y, dpsi_dx, dpsi_dy)
        call assert_close(dpsi_dx, dpsi_dx_fd, tol_fd, "dpsi_1/dx finite diff")

        ! Test dpsi_1/dy
        dpsi_dy_fd = (psi_1(x, y+dy) - psi_1(x, y-dy)) / (2.0_dp * dy)
        call assert_close(dpsi_dy, dpsi_dy_fd, tol_fd, "dpsi_1/dy finite diff")

        ! Test dpsi_2/dx
        dpsi_dx_fd = (psi_2(x+dx, y) - psi_2(x-dx, y)) / (2.0_dp * dx)
        call dpsi_2_derivatives(x, y, dpsi_dx, dpsi_dy)
        call assert_close(dpsi_dx, dpsi_dx_fd, tol_fd, "dpsi_2/dx finite diff")

        ! Test dpsi_2/dy
        dpsi_dy_fd = (psi_2(x, y+dy) - psi_2(x, y-dy)) / (2.0_dp * dy)
        call assert_close(dpsi_dy, dpsi_dy_fd, tol_fd, "dpsi_2/dy finite diff")

        print *, ""
    end subroutine test_basis_derivatives

    subroutine assert_close(actual, expected, tolerance, test_name)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: test_name
        real(dp) :: diff

        test_count = test_count + 1
        diff = abs(actual - expected)

        if (diff < tolerance) then
            print *, "  PASS: ", trim(test_name)
            pass_count = pass_count + 1
        else
            print *, "  FAIL: ", trim(test_name)
            print *, "    Expected: ", expected
            print *, "    Got:      ", actual
            print *, "    Diff:     ", diff
        end if
    end subroutine assert_close

end program test_analytical_gs_circular
