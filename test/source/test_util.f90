program test_util
    use libneo_kinds, only: dp

    implicit none

    logical :: all_passed

    all_passed = .true.

    call test_linspace_basic(all_passed)
    call test_linspace_single_point(all_passed)
    call test_linspace_two_points(all_passed)
    call test_linspace_large(all_passed)
    call test_linspace_negative_range(all_passed)
    call test_linspace_reversed(all_passed)
    call test_linspace_uniform_spacing(all_passed)

    if (.not. all_passed) then
        error stop "One or more util tests failed"
    end if

contains

    subroutine test_linspace_basic(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), dimension(5) :: expected, actual

        expected = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
        call linspace(0.0d0, 4.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > tol) then
            write(*,*) "FAIL: linspace_basic"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linspace_basic"
    end subroutine test_linspace_basic

    subroutine test_linspace_single_point(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), dimension(1) :: actual

        call linspace(5.0d0, 10.0d0, 1, actual)

        if (abs(actual(1) - 5.0d0) > tol) then
            write(*,*) "FAIL: linspace_single_point. Got:", actual(1), " expected 5.0"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linspace_single_point"
    end subroutine test_linspace_single_point

    subroutine test_linspace_two_points(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), dimension(2) :: actual

        call linspace(1.0d0, 3.0d0, 2, actual)

        if (abs(actual(1) - 1.0d0) > tol) then
            write(*,*) "FAIL: linspace_two_points first. Got:", actual(1)
            passed = .false.
            return
        end if
        if (abs(actual(2) - 3.0d0) > tol) then
            write(*,*) "FAIL: linspace_two_points second. Got:", actual(2)
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linspace_two_points"
    end subroutine test_linspace_two_points

    subroutine test_linspace_large(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-10
        integer, parameter :: n = 10001
        real(dp), dimension(n) :: actual
        real(dp) :: expected_first, expected_last
        integer :: i

        call linspace(0.0d0, 100.0d0, n, actual)

        expected_first = 0.0d0
        expected_last = 100.0d0

        if (abs(actual(1) - expected_first) > tol) then
            write(*,*) "FAIL: linspace_large first. Got:", actual(1)
            passed = .false.
            return
        end if
        if (abs(actual(n) - expected_last) > tol) then
            write(*,*) "FAIL: linspace_large last. Got:", actual(n)
            passed = .false.
            return
        end if

        do i = 1, n-1
            if (actual(i+1) <= actual(i)) then
                write(*,*) "FAIL: linspace_large monotonicity at i=", i
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_linspace_large"
    end subroutine test_linspace_large

    subroutine test_linspace_negative_range(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), dimension(5) :: expected, actual

        expected = [-2.0d0, -1.0d0, 0.0d0, 1.0d0, 2.0d0]
        call linspace(-2.0d0, 2.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > tol) then
            write(*,*) "FAIL: linspace_negative_range"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linspace_negative_range"
    end subroutine test_linspace_negative_range

    subroutine test_linspace_reversed(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        real(dp), dimension(5) :: expected, actual

        expected = [4.0d0, 3.0d0, 2.0d0, 1.0d0, 0.0d0]
        call linspace(4.0d0, 0.0d0, 5, actual)

        if (maxval(abs(expected - actual)) > tol) then
            write(*,*) "FAIL: linspace_reversed"
            passed = .false.
            return
        end if

        write(*,*) "PASS: test_linspace_reversed"
    end subroutine test_linspace_reversed

    subroutine test_linspace_uniform_spacing(passed)
        use libneo_util, only: linspace
        logical, intent(inout) :: passed

        real(dp), parameter :: tol = 1.0d-12
        integer, parameter :: n = 101
        real(dp), dimension(n) :: actual
        real(dp) :: expected_spacing, actual_spacing
        integer :: i

        call linspace(0.0d0, 10.0d0, n, actual)

        expected_spacing = 10.0d0 / (n - 1)

        do i = 1, n-1
            actual_spacing = actual(i+1) - actual(i)
            if (abs(actual_spacing - expected_spacing) > tol) then
                write(*,*) "FAIL: linspace_uniform_spacing at i=", i, &
                    " spacing=", actual_spacing, " expected=", expected_spacing
                passed = .false.
                return
            end if
        end do

        write(*,*) "PASS: test_linspace_uniform_spacing"
    end subroutine test_linspace_uniform_spacing

end program test_util
