program test_rng
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_continuous_mean(all_passed)
    call test_continuous_variance(all_passed)
    call test_continuous_range(all_passed)
    call test_discrete_values(all_passed)
    call test_discrete_balance(all_passed)
    call test_continuous_distribution(all_passed)

    if (.not. all_passed) then
        error stop "One or more rng tests failed"
    end if

contains

    subroutine test_continuous_mean(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 100000
        real :: ur, total
        integer :: i

        total = 0.0

        do i = 1, n
            call getran(0, ur)
            total = total + ur
        end do

        if (abs(total / n) > 0.02) then
            write(*,*) "FAIL: continuous mean should be ~0. Got:", total / n
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_continuous_mean"
    end subroutine test_continuous_mean

    subroutine test_continuous_variance(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 100000
        real :: ur, total_sq
        real :: expected_variance
        integer :: i

        total_sq = 0.0
        expected_variance = 1.0

        do i = 1, n
            call getran(0, ur)
            total_sq = total_sq + ur**2
        end do

        if (abs(total_sq / n - expected_variance) > 0.05) then
            write(*,*) "FAIL: continuous variance should be ~1. Got:", total_sq / n
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_continuous_variance"
    end subroutine test_continuous_variance

    subroutine test_continuous_range(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 10000
        real :: ur
        real :: sqrt3
        integer :: i

        sqrt3 = sqrt(3.0)

        do i = 1, n
            call getran(0, ur)
            if (ur < -sqrt3 .or. ur > sqrt3) then
                write(*,*) "FAIL: continuous value out of range [-sqrt(3), sqrt(3)]. Got:", ur
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_continuous_range"
    end subroutine test_continuous_range

    subroutine test_discrete_values(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 1000
        real :: ur
        integer :: i

        do i = 1, n
            call getran(1, ur)
            if (abs(ur - 1.0) > 1.0e-6 .and. abs(ur + 1.0) > 1.0e-6) then
                write(*,*) "FAIL: discrete value should be +1 or -1. Got:", ur
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_discrete_values"
    end subroutine test_discrete_values

    subroutine test_discrete_balance(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 100000
        real :: ur
        integer :: count_plus, count_minus
        real :: ratio
        integer :: i

        count_plus = 0
        count_minus = 0

        do i = 1, n
            call getran(1, ur)
            if (ur > 0.0) then
                count_plus = count_plus + 1
            else
                count_minus = count_minus + 1
            end if
        end do

        ratio = real(count_plus) / real(count_minus)

        if (abs(ratio - 1.0) > 0.05) then
            write(*,*) "FAIL: discrete +1/-1 ratio should be ~1. Got:", ratio
            write(*,*) "  +1 count:", count_plus, " -1 count:", count_minus
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_discrete_balance"
    end subroutine test_discrete_balance

    subroutine test_continuous_distribution(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 100000
        integer, parameter :: nbins = 10
        real :: ur
        integer :: bins(nbins)
        real :: bin_width, sqrt3
        integer :: i, bin_idx
        real :: expected_count, chi_sq
        integer :: min_count, max_count

        sqrt3 = sqrt(3.0)
        bin_width = 2.0 * sqrt3 / nbins
        bins = 0

        do i = 1, n
            call getran(0, ur)
            bin_idx = int((ur + sqrt3) / bin_width) + 1
            if (bin_idx < 1) bin_idx = 1
            if (bin_idx > nbins) bin_idx = nbins
            bins(bin_idx) = bins(bin_idx) + 1
        end do

        expected_count = real(n) / real(nbins)
        chi_sq = 0.0
        min_count = bins(1)
        max_count = bins(1)

        do i = 1, nbins
            chi_sq = chi_sq + (bins(i) - expected_count)**2 / expected_count
            if (bins(i) < min_count) min_count = bins(i)
            if (bins(i) > max_count) max_count = bins(i)
        end do

        if (chi_sq > 30.0) then
            write(*,*) "FAIL: chi-squared too high for uniform distribution:", chi_sq
            write(*,*) "  Expected ~", nbins - 1, " (degrees of freedom)"
            passed = .false.
            return
        end if

        if (min_count < expected_count * 0.8 .or. max_count > expected_count * 1.2) then
            write(*,*) "FAIL: bin counts too uneven"
            write(*,*) "  Min:", min_count, " Max:", max_count, " Expected:", int(expected_count)
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_continuous_distribution"
    end subroutine test_continuous_distribution

end program test_rng
