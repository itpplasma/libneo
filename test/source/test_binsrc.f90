program test_binsrc
    use binsrc_sub, only: binsrc
    use libneo_kinds, only: dp

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_value_between_bounds(all_passed)
    call test_value_below_minimum(all_passed)
    call test_value_above_maximum(all_passed)
    call test_exact_match(all_passed)
    call test_exact_match_first(all_passed)
    call test_exact_match_last(all_passed)
    call test_negative_values(all_passed)
    call test_large_array(all_passed)
    call test_small_spacing(all_passed)
    call test_nonuniform_spacing(all_passed)
    call test_two_element_array(all_passed)
    call test_fractional_values(all_passed)
    call test_wide_range(all_passed)

    if (.not. all_passed) then
        error stop "One or more binsrc tests failed"
    end if

contains

    subroutine test_value_between_bounds(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 5.5_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= 6) then
            write(*,*) "FAIL: value_between_bounds. Found:", found_index, " expected 6"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_value_between_bounds"
    end subroutine test_value_between_bounds

    subroutine test_value_below_minimum(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 0.0_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= min_index + 1) then
            write(*,*) "FAIL: value_below_minimum. Found:", found_index, &
                " expected", min_index + 1
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_value_below_minimum"
    end subroutine test_value_below_minimum

    subroutine test_value_above_maximum(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 21.0_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= max_index) then
            write(*,*) "FAIL: value_above_maximum. Found:", found_index, &
                " expected", max_index
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_value_above_maximum"
    end subroutine test_value_above_maximum

    subroutine test_exact_match(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 10.0_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= 11) then
            write(*,*) "FAIL: exact_match. Found:", found_index, " expected 11"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_exact_match"
    end subroutine test_exact_match

    subroutine test_exact_match_first(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 1.0_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= 2) then
            write(*,*) "FAIL: exact_match_first. Found:", found_index, " expected 2"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_exact_match_first"
    end subroutine test_exact_match_first

    subroutine test_exact_match_last(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 20
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k, dp), k=min_index, max_index)]

        value_to_find = 20.0_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= 20) then
            write(*,*) "FAIL: exact_match_last. Found:", found_index, " expected 20"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_exact_match_last"
    end subroutine test_exact_match_last

    subroutine test_negative_values(passed)
        logical, intent(inout) :: passed

        integer, parameter :: min_index = 1, max_index = 11
        real(dp) :: values(min_index:max_index)
        real(dp) :: value_to_find
        integer :: found_index, k

        values = [(real(k - 6, dp), k=min_index, max_index)]

        value_to_find = -2.5_dp
        call binsrc(values, min_index, max_index, value_to_find, found_index)
        if (found_index /= 4) then
            write(*,*) "FAIL: negative_values. Found:", found_index, " expected 4"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_negative_values"
    end subroutine test_negative_values

    subroutine test_large_array(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 10000
        real(dp), allocatable :: values(:)
        real(dp) :: value_to_find
        integer :: found_index, k

        allocate(values(n))
        values = [(real(k, dp), k=1, n)]

        value_to_find = 7777.5_dp
        call binsrc(values, 1, n, value_to_find, found_index)
        if (found_index /= 7778) then
            write(*,*) "FAIL: large_array. Found:", found_index, " expected 7778"
            passed = .false.
            deallocate(values)
            return
        end if

        deallocate(values)
        write(*,*) "PASS: test_large_array"
    end subroutine test_large_array

    subroutine test_small_spacing(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 100
        real(dp) :: values(n)
        real(dp) :: value_to_find
        integer :: found_index, k

        do k = 1, n
            values(k) = 1.0_dp + (k - 1) * 1.0d-10
        end do

        value_to_find = 1.0_dp + 50.5d-10
        call binsrc(values, 1, n, value_to_find, found_index)
        if (found_index /= 52) then
            write(*,*) "FAIL: small_spacing. Found:", found_index, " expected 52"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_small_spacing"
    end subroutine test_small_spacing

    subroutine test_nonuniform_spacing(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 10
        real(dp) :: values(n)
        real(dp) :: value_to_find
        integer :: found_index, k

        do k = 1, n
            values(k) = real(k, dp)**2
        end do

        value_to_find = 50.0_dp
        call binsrc(values, 1, n, value_to_find, found_index)
        if (found_index /= 8) then
            write(*,*) "FAIL: nonuniform_spacing. Found:", found_index, " expected 8"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_nonuniform_spacing"
    end subroutine test_nonuniform_spacing

    subroutine test_two_element_array(passed)
        logical, intent(inout) :: passed

        real(dp) :: values(2)
        real(dp) :: value_to_find
        integer :: found_index

        values = [1.0_dp, 2.0_dp]

        value_to_find = 1.5_dp
        call binsrc(values, 1, 2, value_to_find, found_index)
        if (found_index /= 2) then
            write(*,*) "FAIL: two_element_array middle. Found:", found_index, &
                " expected 2"
            passed = .false.
            return
        end if

        value_to_find = 2.5_dp
        call binsrc(values, 1, 2, value_to_find, found_index)
        if (found_index /= 2) then
            write(*,*) "FAIL: two_element_array above. Found:", found_index, &
                " expected 2"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_two_element_array"
    end subroutine test_two_element_array

    subroutine test_fractional_values(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 10
        real(dp) :: values(n)
        real(dp) :: value_to_find
        integer :: found_index, k

        do k = 1, n
            values(k) = 0.1_dp * k
        end do

        value_to_find = 0.55_dp
        call binsrc(values, 1, n, value_to_find, found_index)
        if (found_index /= 6) then
            write(*,*) "FAIL: fractional_values. Found:", found_index, " expected 6"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_fractional_values"
    end subroutine test_fractional_values

    subroutine test_wide_range(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 5
        real(dp) :: values(n)
        real(dp) :: value_to_find
        integer :: found_index

        values = [1.0d-10, 1.0d-5, 1.0_dp, 1.0d5, 1.0d10]

        value_to_find = 1.0d3
        call binsrc(values, 1, n, value_to_find, found_index)
        if (found_index /= 4) then
            write(*,*) "FAIL: wide_range. Found:", found_index, " expected 4"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_wide_range"
    end subroutine test_wide_range

end program test_binsrc
