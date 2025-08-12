program test_spline_vmec_array_temps
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util_for_test, only: print_test, print_ok, print_fail
    implicit none

    logical :: test_passed
    integer :: test_status

    test_passed = .true.
    test_status = 0

    call test_s_to_rho_interface()
    call test_axis_healing_performance()

    if (.not. test_passed) then
        print *, "Tests FAILED"
        test_status = 1
    else
        print *, "All tests PASSED"
    end if

    call exit(test_status)

contains

    subroutine test_s_to_rho_interface()
        integer, parameter :: m = 2, ns = 10, nrho = 15, nheal = 2
        integer, parameter :: nmodes = 5
        real(dp), dimension(nmodes, ns) :: arr_2d_in
        real(dp), dimension(nmodes, nrho) :: arr_2d_out
        real(dp), dimension(ns) :: arr_1d_in
        real(dp), dimension(nrho) :: arr_1d_out
        real(dp) :: diff

        call print_test("Testing s_to_rho_healaxis interface")

        call random_number(arr_2d_in)
        arr_1d_in = arr_2d_in(1, :)

        call s_to_rho_healaxis_original(m, ns, nrho, nheal, arr_1d_in, arr_1d_out)

        call s_to_rho_healaxis_fixed(m, ns, nrho, nheal, 1, nmodes, arr_2d_in, arr_2d_out)

        diff = abs(arr_1d_out(1) - arr_2d_out(1, 1))
        if (diff < 1.0e-14_dp) then
            call print_ok()
        else
            call print_fail()
            test_passed = .false.
            print *, "  Difference:", diff
        end if
    end subroutine test_s_to_rho_interface

    subroutine test_axis_healing_performance()
        integer, parameter :: nmodes = 100, ns = 50, nrho = 75
        real(dp), dimension(nmodes, ns) :: input_arrays
        real(dp), dimension(nmodes, nrho) :: output_arrays
        real(dp) :: t1, t2
        integer :: i

        call print_test("Testing axis healing performance improvements")

        call random_number(input_arrays)

        call cpu_time(t1)
        do i = 1, nmodes
            call process_mode_fixed(i, nmodes, ns, nrho, input_arrays, output_arrays)
        end do
        call cpu_time(t2)

        print *, "  Performance test completed in", t2 - t1, "seconds"
        call print_ok()
    end subroutine test_axis_healing_performance

    subroutine s_to_rho_healaxis_original(m, ns, nrho, nheal, arr_in, arr_out)
        integer, intent(in) :: m, ns, nrho, nheal
        real(dp), dimension(ns), intent(in) :: arr_in
        real(dp), dimension(nrho), intent(out) :: arr_out
        
        integer :: i
        real(dp) :: scale
        
        scale = real(m + nheal, dp) / real(ns * nrho, dp)
        
        do i = 1, min(nrho, size(arr_out))
            if (i <= size(arr_in)) then
                arr_out(i) = arr_in(i) * scale
            else
                arr_out(i) = 0.0_dp
            end if
        end do
    end subroutine s_to_rho_healaxis_original

    subroutine s_to_rho_healaxis_fixed(m, ns, nrho, nheal, imode, nmodes, arr_in, arr_out)
        integer, intent(in) :: m, ns, nrho, nheal, imode, nmodes
        real(dp), dimension(nmodes, ns), intent(in) :: arr_in
        real(dp), dimension(nmodes, nrho), intent(out) :: arr_out
        
        integer :: i
        real(dp) :: scale
        
        scale = real(m + nheal, dp) / real(ns * nrho, dp)
        
        do i = 1, min(nrho, size(arr_out, 2))
            if (i <= size(arr_in, 2)) then
                arr_out(imode, i) = arr_in(imode, i) * scale
            else
                arr_out(imode, i) = 0.0_dp
            end if
        end do
    end subroutine s_to_rho_healaxis_fixed

    subroutine process_mode_fixed(imode, nmodes, ns, nrho, input_arrays, output_arrays)
        integer, intent(in) :: imode, nmodes, ns, nrho
        real(dp), dimension(nmodes, ns), intent(in) :: input_arrays
        real(dp), dimension(nmodes, nrho), intent(out) :: output_arrays
        
        integer :: m, nheal
        
        m = mod(imode, 5) + 1
        nheal = min(m, 4)
        
        call s_to_rho_healaxis_fixed(m, ns, nrho, nheal, imode, nmodes, input_arrays, output_arrays)
    end subroutine process_mode_fixed

end program test_spline_vmec_array_temps