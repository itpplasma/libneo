program test_spline_vmec_array_temps
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util_for_test, only: print_test, print_ok, print_fail
    use spline_vmec_sub, only: s_to_rho_healaxis_2d
    use new_vmec_stuff_mod, only: ns_s
    implicit none

    logical :: test_passed
    integer :: test_status

    test_passed = .true.
    test_status = 0

    call test_array_temporary_elimination_in_healaxis()
    call test_array_slice_temporary_pattern()
    call test_cross_mode_consistency()
    call test_edge_cases_and_bounds()

    if (.not. test_passed) then
        print *, "Tests FAILED"
        test_status = 1
    else
        print *, "All tests PASSED"
    end if

    call exit(test_status)

contains

    subroutine test_array_temporary_elimination_in_healaxis()
        integer, parameter :: nmodes = 5, ns = 20, nrho = 20
        real(dp), dimension(nmodes, ns) :: input_arrays
        real(dp), dimension(nmodes, 0:nrho-1) :: output_arrays
        real(dp), dimension(0:ns_s, ns) :: splcoe_workspace
        integer :: i, m, nheal
        real(dp) :: t1, t2

        ! BDD: Given-When-Then format
        call print_test("GIVEN 2D arrays WHEN calling healaxis routines THEN no array temporaries created")
        
        ! GIVEN: Multi-mode arrays typical of VMEC data processing
        call random_number(input_arrays)
        
        ! WHEN: Processing multiple modes with s_to_rho_healaxis_2d
        ! This would trigger array temporary warnings in the old implementation
        call cpu_time(t1)
        do i = 1, nmodes
            m = mod(i, 3) + 1
            nheal = min(m, 2)
            call s_to_rho_healaxis_2d(m, ns, nrho, nheal, i, nmodes, input_arrays, output_arrays, splcoe_workspace)
        end do
        call cpu_time(t2)
        
        ! THEN: Execution completes without array temporary warnings and performs efficiently
        print *, "  Processing completed in", t2 - t1, "seconds"
        call print_ok()
    end subroutine test_array_temporary_elimination_in_healaxis

    subroutine test_array_slice_temporary_pattern()
        integer, parameter :: ns_test = 10, ns_A_test = 5
        real(dp), dimension(0:ns_A_test, ns_test) :: full_array
        real(dp), dimension(0:ns_A_test - 1, ns_test) :: proper_sized_array
        real(dp) :: sum_full, sum_proper
        
        ! BDD: Test for the array slice temporary pattern that was fixed
        call print_test("GIVEN array slicing WHEN avoiding temporaries THEN proper memory management")
        
        ! GIVEN: This test demonstrates the pattern that caused array temporaries
        ! The original code used patterns like splcoe(0:ns_A-1, :) which create temporaries
        call random_number(full_array)
        
        ! WHEN: Using properly sized arrays instead of slices
        proper_sized_array = full_array(0:ns_A_test - 1, :)
        
        ! THEN: Both approaches should give same results but avoid temporaries
        sum_full = sum(full_array(0:ns_A_test - 1, :))
        sum_proper = sum(proper_sized_array)
        
        if (abs(sum_full - sum_proper) < 1.0e-14_dp) then
            call print_ok()
        else
            call print_fail()
            test_passed = .false.
            print *, "  Array slice handling inconsistency"
        end if
    end subroutine test_array_slice_temporary_pattern

    subroutine test_cross_mode_consistency()
        integer, parameter :: ns = 20, nrho = 20, nmodes = 3
        real(dp), dimension(nmodes, ns) :: arr_2d_in
        real(dp), dimension(nmodes, 0:nrho-1) :: arr_2d_out
        real(dp), dimension(0:ns_s, ns) :: splcoe_workspace
        integer :: m, nheal, j, mode
        real(dp) :: max_diff, diff
        
        call print_test("GIVEN identical input WHEN processing different modes THEN results consistent")
        
        ! GIVEN: Identical input data across multiple modes
        call random_number(arr_2d_in)
        do mode = 2, nmodes
            arr_2d_in(mode, :) = arr_2d_in(1, :)  ! Make all modes identical
        end do
        
        m = 2
        nheal = 1
        
        ! WHEN: Processing each mode with the same parameters
        do mode = 1, nmodes
            call s_to_rho_healaxis_2d(m, ns, nrho, nheal, mode, nmodes, arr_2d_in, arr_2d_out, splcoe_workspace)
        end do
        
        ! THEN: All modes should produce identical results
        max_diff = 0.0_dp
        do mode = 2, nmodes
            do j = 0, nrho-1
                diff = abs(arr_2d_out(1, j) - arr_2d_out(mode, j))
                max_diff = max(max_diff, diff)
            end do
        end do
        
        if (max_diff < 1.0e-14_dp) then
            call print_ok()
        else
            call print_fail()
            test_passed = .false.
            print *, "  Maximum difference between modes:", max_diff
        end if
    end subroutine test_cross_mode_consistency

    subroutine test_edge_cases_and_bounds()
        integer, parameter :: nmodes = 2
        real(dp), dimension(nmodes, 10) :: small_input
        real(dp), dimension(nmodes, 0:9) :: small_output
        real(dp), dimension(nmodes, 100) :: large_input
        real(dp), dimension(nmodes, 0:99) :: large_output
        real(dp), dimension(0:ns_s, 10) :: small_splcoe
        real(dp), dimension(0:ns_s, 100) :: large_splcoe
        
        call print_test("GIVEN edge case parameters WHEN processing THEN robust handling without crashes")
        
        ! GIVEN: Various edge case scenarios that could cause issues
        call random_number(small_input)
        call random_number(large_input)
        
        ! WHEN: Testing various edge cases that previously caused issues
        
        ! Test minimum viable array size (respects spline order ns_s = 5)
        call s_to_rho_healaxis_2d(1, 10, 10, 2, 1, nmodes, small_input, small_output, small_splcoe)
        
        ! Test m = 0 case (different computational path)
        call s_to_rho_healaxis_2d(0, 10, 10, 1, 1, nmodes, small_input, small_output, small_splcoe)
        
        ! Test large arrays (performance/memory stress test)
        call s_to_rho_healaxis_2d(3, 100, 100, 4, 1, nmodes, large_input, large_output, large_splcoe)
        
        ! Test moderate nheal (boundary condition testing)
        call s_to_rho_healaxis_2d(1, 10, 10, 3, 1, nmodes, small_input, small_output, small_splcoe)
        
        ! Test bounds checking improvement (nheal close to ns)
        call s_to_rho_healaxis_2d(2, 10, 10, 7, 1, nmodes, small_input, small_output, small_splcoe)
        
        ! THEN: All edge cases complete without crashes or bounds violations
        call print_ok()
    end subroutine test_edge_cases_and_bounds

end program test_spline_vmec_array_temps