program test_plag_coeff
    use libneo_kinds, only: dp
    use plag_coeff_sub, only: plag_coeff

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_interpolation_property(all_passed)
    call test_partition_of_unity(all_passed)
    call test_polynomial_interpolation(all_passed)
    call test_derivative_computation(all_passed)
    call test_different_orders(all_passed)
    call test_nonuniform_grid(all_passed)

    if (.not. all_passed) then
        error stop "One or more plag_coeff tests failed"
    end if

contains

    subroutine test_interpolation_property(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        integer, parameter :: npoi = 5
        real(dp) :: xp(npoi), coef(0:0,npoi)
        integer :: i, j

        xp = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]

        do i = 1, npoi
            call plag_coeff(npoi, 0, xp(i), xp, coef)
            do j = 1, npoi
                if (i == j) then
                    if (abs(coef(0,j) - 1.0d0) > tol) then
                        write(*,*) "FAIL: interpolation_property L_", j, "(x_", i, &
                            ") = ", coef(0,j), " expected 1.0"
                        passed = .false.
                        return
                    end if
                else
                    if (abs(coef(0,j)) > tol) then
                        write(*,*) "FAIL: interpolation_property L_", j, "(x_", i, &
                            ") = ", coef(0,j), " expected 0.0"
                        passed = .false.
                        return
                    end if
                end if
            end do
        end do

        write(*,*) "PASS: test_interpolation_property"
    end subroutine test_interpolation_property

    subroutine test_partition_of_unity(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        integer, parameter :: npoi = 5
        real(dp) :: xp(npoi), coef(0:0,npoi)
        real(dp) :: x, coef_sum
        integer :: i

        xp = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]

        do i = 1, 10
            x = 0.37d0 + 0.41d0 * i
            call plag_coeff(npoi, 0, x, xp, coef)
            coef_sum = sum(coef(0,:))
            if (abs(coef_sum - 1.0d0) > tol) then
                write(*,*) "FAIL: partition_of_unity at x=", x, " sum=", coef_sum
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_partition_of_unity"
    end subroutine test_partition_of_unity

    subroutine test_polynomial_interpolation(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: npoi = 5
        real(dp) :: xp(npoi), coef(0:1,npoi)
        real(dp) :: fvals(npoi)
        real(dp) :: x, f_interp, df_interp, f_exact, df_exact
        integer :: i

        xp = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]

        do i = 1, npoi
            fvals(i) = xp(i)**3 - 2.0d0*xp(i)**2 + xp(i) - 1.0d0
        end do

        x = 1.7d0
        f_exact = x**3 - 2.0d0*x**2 + x - 1.0d0
        df_exact = 3.0d0*x**2 - 4.0d0*x + 1.0d0

        call plag_coeff(npoi, 1, x, xp, coef)

        f_interp = sum(fvals * coef(0,:))
        df_interp = sum(fvals * coef(1,:))

        if (abs(f_interp - f_exact) > tol) then
            write(*,*) "FAIL: polynomial_interpolation f(", x, ") = ", f_interp, &
                " expected ", f_exact
            passed = .false.
            return
        end if

        if (abs(df_interp - df_exact) > tol) then
            write(*,*) "FAIL: polynomial_interpolation df(", x, ") = ", df_interp, &
                " expected ", df_exact
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_polynomial_interpolation"
    end subroutine test_polynomial_interpolation

    subroutine test_derivative_computation(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: npoi = 6
        real(dp) :: xp(npoi), coef(0:1,npoi)
        real(dp) :: fvals(npoi)
        real(dp) :: x, df_interp, df_exact
        integer :: i

        xp = [-1.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]

        do i = 1, npoi
            fvals(i) = sin(xp(i))
        end do

        x = 1.5d0
        df_exact = cos(x)

        call plag_coeff(npoi, 1, x, xp, coef)
        df_interp = sum(fvals * coef(1,:))

        if (abs(df_interp - df_exact) > 1.0d-4) then
            write(*,*) "FAIL: derivative_computation df(", x, ") = ", df_interp, &
                " expected ", df_exact
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_derivative_computation"
    end subroutine test_derivative_computation

    subroutine test_different_orders(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), allocatable :: xp(:), coef(:,:)
        real(dp) :: x, coef_sum
        integer :: npoi

        do npoi = 2, 8
            allocate(xp(npoi), coef(0:0,npoi))

            call linspace_local(0.0d0, 1.0d0, npoi, xp)

            x = 0.37d0
            call plag_coeff(npoi, 0, x, xp, coef)
            coef_sum = sum(coef(0,:))

            if (abs(coef_sum - 1.0d0) > tol) then
                write(*,*) "FAIL: different_orders npoi=", npoi, " sum=", coef_sum
                passed = .false.
                deallocate(xp, coef)
                return
            end if

            deallocate(xp, coef)
        end do

        write(*,*) "PASS: test_different_orders"
    end subroutine test_different_orders

    subroutine test_nonuniform_grid(passed)
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: npoi = 5
        real(dp) :: xp(npoi), coef(0:1,npoi)
        real(dp) :: fvals(npoi)
        real(dp) :: x, f_interp, f_exact
        integer :: i

        xp = [0.0d0, 0.1d0, 0.5d0, 0.9d0, 1.0d0]

        do i = 1, npoi
            fvals(i) = xp(i)**2
        end do

        x = 0.73d0
        f_exact = x**2

        call plag_coeff(npoi, 1, x, xp, coef)
        f_interp = sum(fvals * coef(0,:))

        if (abs(f_interp - f_exact) > tol) then
            write(*,*) "FAIL: nonuniform_grid f(", x, ") = ", f_interp, &
                " expected ", f_exact
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_nonuniform_grid"
    end subroutine test_nonuniform_grid

    subroutine linspace_local(a, b, n, x)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n)
        real(dp) :: dx
        integer :: i

        if (n == 1) then
            x(1) = a
            return
        end if
        dx = (b - a) / (n - 1)
        do i = 1, n
            x(i) = a + (i - 1) * dx
        end do
    end subroutine linspace_local

end program test_plag_coeff
