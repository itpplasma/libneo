program test_simpson
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simpson_integration, only: simpson_nonequi, simpson_prefix

    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    logical :: all_passed

    all_passed = .true.

    call test_sin_integral(all_passed)
    call test_polynomial_x2(all_passed)
    call test_polynomial_x3(all_passed)
    call test_polynomial_x4(all_passed)
    call test_exponential(all_passed)
    call test_gaussian(all_passed)
    call test_constant_function(all_passed)
    call test_linear_function(all_passed)
    call test_prefix_monotonicity(all_passed)
    call test_nonuniform_grid(all_passed)
    call test_many_points(all_passed)
    call test_few_points(all_passed)

    if (.not. all_passed) then
        error stop "One or more simpson tests failed"
    end if

contains

    subroutine test_sin_integral(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 55
        real(dp) :: approx, x(n), f(n), prefix(n)
        integer :: i

        do i = 1, n
            x(i) = 0.5_dp * pi * (1.0_dp - cos((i - 1.0_dp) * pi / (n - 1.0_dp)))
            f(i) = sin(x(i))
        end do

        call simpson_prefix(prefix, x, f)
        if (abs(prefix(n) - 2.0_dp) > 1.0e-6_dp) then
            write(*,*) "FAIL: sin_integral prefix. Value:", prefix(n), " expected 2.0"
            passed = .false.
            return
        end if

        call simpson_nonequi(approx, x, f)
        if (abs(approx - 2.0_dp) > 1.0e-6_dp) then
            write(*,*) "FAIL: sin_integral total. Value:", approx, " expected 2.0"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_sin_integral"
    end subroutine test_sin_integral

    subroutine test_polynomial_x2(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 101
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = x(i)**2
        end do

        expected = 1.0_dp / 3.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-8_dp) then
            write(*,*) "FAIL: polynomial_x2. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_polynomial_x2"
    end subroutine test_polynomial_x2

    subroutine test_polynomial_x3(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 101
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = x(i)**3
        end do

        expected = 1.0_dp / 4.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-8_dp) then
            write(*,*) "FAIL: polynomial_x3. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_polynomial_x3"
    end subroutine test_polynomial_x3

    subroutine test_polynomial_x4(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 201
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = x(i)**4
        end do

        expected = 1.0_dp / 5.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-6_dp) then
            write(*,*) "FAIL: polynomial_x4. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_polynomial_x4"
    end subroutine test_polynomial_x4

    subroutine test_exponential(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 101
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = exp(x(i))
        end do

        expected = exp(1.0_dp) - 1.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-8_dp) then
            write(*,*) "FAIL: exponential. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_exponential"
    end subroutine test_exponential

    subroutine test_gaussian(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 201
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = -3.0_dp + 6.0_dp * (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = exp(-x(i)**2)
        end do

        expected = sqrt(pi) * erf(3.0_dp)

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-6_dp) then
            write(*,*) "FAIL: gaussian. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_gaussian"
    end subroutine test_gaussian

    subroutine test_constant_function(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 17
        real(dp) :: approx, x(n), f(n)
        real(dp), parameter :: c = 3.7_dp
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = 2.0_dp + 5.0_dp * (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = c
        end do

        expected = c * 5.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-12_dp) then
            write(*,*) "FAIL: constant_function. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_constant_function"
    end subroutine test_constant_function

    subroutine test_linear_function(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 11
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = 1.0_dp + 3.0_dp * (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = 2.0_dp * x(i) + 1.0_dp
        end do

        expected = (4.0_dp**2 + 4.0_dp) - (1.0_dp**2 + 1.0_dp)

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-12_dp) then
            write(*,*) "FAIL: linear_function. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linear_function"
    end subroutine test_linear_function

    subroutine test_prefix_monotonicity(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 51
        real(dp) :: x(n), f(n), prefix(n)
        integer :: i

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = x(i)**2
        end do

        call simpson_prefix(prefix, x, f)

        do i = 2, n
            if (prefix(i) < prefix(i-1) - 1.0e-15_dp) then
                write(*,*) "FAIL: prefix_monotonicity at i=", i, &
                    " prefix(i)=", prefix(i), " prefix(i-1)=", prefix(i-1)
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_prefix_monotonicity"
    end subroutine test_prefix_monotonicity

    subroutine test_nonuniform_grid(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 51
        real(dp) :: approx, x(n), f(n)
        real(dp) :: expected
        integer :: i

        do i = 1, n
            x(i) = ((i - 1.0_dp) / (n - 1.0_dp))**2
            f(i) = exp(x(i))
        end do

        expected = exp(1.0_dp) - 1.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-4_dp) then
            write(*,*) "FAIL: nonuniform_grid. Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_nonuniform_grid"
    end subroutine test_nonuniform_grid

    subroutine test_many_points(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 10001
        real(dp) :: approx
        real(dp), allocatable :: x(:), f(:)
        real(dp) :: expected
        integer :: i

        allocate(x(n), f(n))

        do i = 1, n
            x(i) = (i - 1.0_dp) / (n - 1.0_dp)
            f(i) = sin(pi * x(i))
        end do

        expected = 2.0_dp / pi

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-6_dp) then
            write(*,*) "FAIL: many_points. Value:", approx, " expected ", expected
            passed = .false.
            deallocate(x, f)
            return
        end if

        deallocate(x, f)
        write(*,*) "PASS: test_many_points"
    end subroutine test_many_points

    subroutine test_few_points(passed)
        logical, intent(inout) :: passed

        real(dp) :: approx, x(3), f(3)
        real(dp) :: expected
        integer :: i

        do i = 1, 3
            x(i) = (i - 1.0_dp) / 2.0_dp
            f(i) = x(i)**2
        end do

        expected = 1.0_dp / 3.0_dp

        call simpson_nonequi(approx, x, f)
        if (abs(approx - expected) > 1.0e-14_dp) then
            write(*,*) "FAIL: few_points (3). Value:", approx, " expected ", expected
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_few_points"
    end subroutine test_few_points

end program test_simpson
