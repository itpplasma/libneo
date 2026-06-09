program test_dawson

    use libneo_kinds, only: dp
    use neo_dawson, only: dawson
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_reference_values
    call test_zero
    call test_odd_symmetry
    call test_asymptotic_tail_exact
    call test_elemental_interface

contains

    subroutine test_reference_values

        ! reference values from mpmath at >= 60 digits, evaluated at the exact
        ! binary double value of each x; for x > 20 via the asymptotic series
        ! with remainder < e^(-x^2) < 1e-170

        integer, parameter :: n = 54
        real(dp), dimension(n), parameter :: x_ref = [ &
            1d-308, 1d-300, &
            1d-200, 1d-100, &
            1d-30, 1d-15, &
            1d-08, 0.0001d0, &
            0.01d0, 0.1d0, &
            0.25d0, 0.5d0, &
            0.75d0, 0.9d0, &
            0.9241388730045917d0, 0.999999d0, &
            1.0d0, 1.0000001d0, &
            1.3d0, 1.5d0, &
            1.7d0, 2.0d0, &
            2.5d0, 3.0d0, &
            3.5d0, 4.0d0, &
            4.5d0, 5.0d0, &
            5.5d0, 6.0d0, &
            6.5d0, 7.0d0, &
            7.5d0, 8.0d0, &
            8.5d0, 9.0d0, &
            9.5d0, 9.999999d0, &
            10.0d0, 10.5d0, &
            11.0d0, 12.0d0, &
            15.0d0, 20.0d0, &
            50.0d0, 100.0d0, &
            1000.0d0, 1000000.0d0, &
            100000000.0d0, 10000000000.0d0, &
            1d+20, 1d+100, &
            1d+200, 1d+300]
        real(dp), dimension(n), parameter :: f_ref = [ &
            9.9999999999999991d-309, 1.0d-300, &
            9.9999999999999998d-201, 1.0d-100, &
            1.0000000000000001d-30, 1.0000000000000001d-15, &
            9.9999999999999995d-9, 9.9999999333333341d-5, &
            0.0099993333599992383d0, 0.099335992397852867d0, &
            0.23983916356289821d0, 0.4244363835020223d0, &
            0.52301276774451825d0, 0.54072431872629868d0, &
            0.5410442246351817d0, 0.53807958307132033d0, &
            0.53807950691276842d0, 0.53807949929686241d0, &
            0.48339751738482412d0, 0.42824907108539863d0, &
            0.37255934897407884d0, 0.30134038892379197d0, &
            0.22308372216743548d0, 0.17827103061055829d0, &
            0.14962159308075648d0, 0.12934800123600512d0, &
            0.11408861022682498d0, 0.10213407442427684d0, &
            0.09249323231075476d0, 0.084542688974543852d0, &
            0.077867818986069871d0, 0.072180974658236292d0, &
            0.067275811644630616d0, 0.063000198707553388d0, &
            0.05923937177997214d0, 0.055905046724350461d0, &
            0.052928152705625646d0, 0.050253852264542792d0, &
            0.050253847187598528d0, 0.047838014074213438d0, &
            0.04564475216411602d0, 0.04181287645398826d0, &
            0.033407906808639226d0, 0.025031367926403672d0, &
            0.010002001201201683d0, 0.0050002500375093783d0, &
            0.000500000250000375d0, 5.0000000000025d-7, &
            5.0000000000000002d-9, 5.0d-11, &
            5.0d-21, 4.9999999999999999d-101, &
            5.0000000000000002d-201, 4.9999999999999997d-301]
        real(dp), parameter :: rel_tol = 1.0d-14

        real(dp) :: f, rel_err
        integer :: i

        call print_test("test_reference_values")

        do i = 1, n
            f = dawson(x_ref(i))
            rel_err = abs(f - f_ref(i))/abs(f_ref(i))
            if (rel_err > rel_tol) then
                call print_fail
                print *, "x = ", x_ref(i)
                print *, "got ", f, ", expected ", f_ref(i), ", rel err ", rel_err
                stop "dawson deviates from mpmath reference"
            end if
        end do
        call print_ok

    end subroutine test_reference_values

    subroutine test_zero

        call print_test("test_zero")
        if (dawson(0.0d0) /= 0.0d0) then
            call print_fail
            stop "dawson(0) /= 0"
        end if
        call print_ok

    end subroutine test_zero

    subroutine test_odd_symmetry

        integer, parameter :: n = 8
        real(dp), dimension(n), parameter :: x_ref = [ &
            1.0d-12, 0.3d0, 0.99d0, 1.0d0, 2.7d0, 9.5d0, 42.0d0, 1.0d250]
        integer :: i

        call print_test("test_odd_symmetry")

        do i = 1, n
            if (dawson(-x_ref(i)) /= -dawson(x_ref(i))) then
                call print_fail
                print *, "x = ", x_ref(i)
                stop "dawson is not odd to the bit"
            end if
        end do
        call print_ok

    end subroutine test_odd_symmetry

    subroutine test_asymptotic_tail_exact

        ! for x >= 1e9 the correction 1/(2x^2) is below 2^-53, so the
        ! double-precision result must equal 0.5/x exactly

        integer, parameter :: n = 5
        real(dp), dimension(n), parameter :: x_ref = [ &
            1.0d9, 1.0d16, 1.0d100, 1.0d200, 1.0d300]
        integer :: i

        call print_test("test_asymptotic_tail_exact")

        do i = 1, n
            if (dawson(x_ref(i)) /= 0.5d0/x_ref(i)) then
                call print_fail
                print *, "x = ", x_ref(i), ", got ", dawson(x_ref(i))
                stop "asymptotic tail not exact to fp"
            end if
        end do
        call print_ok

    end subroutine test_asymptotic_tail_exact

    subroutine test_elemental_interface

        real(dp), dimension(4), parameter :: x_arr = [0.1d0, 1.5d0, 5.0d0, 20.0d0]
        real(dp), dimension(4) :: f_arr
        integer :: i

        call print_test("test_elemental_interface")

        f_arr = dawson(x_arr)
        do i = 1, 4
            if (f_arr(i) /= dawson(x_arr(i))) then
                call print_fail
                stop "elemental array result differs from scalar call"
            end if
        end do
        call print_ok

    end subroutine test_elemental_interface

end program test_dawson
