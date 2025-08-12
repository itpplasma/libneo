program test_spline_vmec_array_temps
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util_for_test, only: print_test, print_ok, print_fail
    use spline_vmec_sub, only: s_to_rho_healaxis_2d
    implicit none

    logical :: test_passed
    integer :: test_status

    test_passed = .true.
    test_status = 0

    call test_array_temporary_elimination()
    call test_consistency()
    call test_edge_cases()

    if (.not. test_passed) then
        print *, "Tests FAILED"
        test_status = 1
    else
        print *, "All tests PASSED"
    end if

    call exit(test_status)

contains

    subroutine test_array_temporary_elimination()
        integer, parameter :: nmodes = 5, ns = 20, nrho = 20
        real(dp), dimension(nmodes, ns) :: input_arrays
        real(dp), dimension(nmodes, 0:nrho-1) :: output_arrays
        integer :: i, m, nheal
        real(dp) :: t1, t2

        call print_test("Testing array temporary elimination")

        ! Initialize with test data
        call random_number(input_arrays)
        
        ! This test primarily validates that the new interface works without
        ! creating array temporaries when compiled with -fcheck=array-temps
        call cpu_time(t1)
        do i = 1, nmodes
            m = mod(i, 3) + 1
            nheal = min(m, 2)
            call s_to_rho_healaxis_2d(m, ns, nrho, nheal, i, nmodes, input_arrays, output_arrays)
        end do
        call cpu_time(t2)
        
        print *, "  2D version completed in", t2 - t1, "seconds"
        call print_ok()
    end subroutine test_array_temporary_elimination

    subroutine test_consistency()
        integer, parameter :: ns = 20, nrho = 20, nmodes = 3
        real(dp), dimension(nmodes, ns) :: arr_2d_in
        real(dp), dimension(nmodes, 0:nrho-1) :: arr_2d_out1
        integer :: m, nheal, j, mode
        real(dp) :: max_diff, diff
        
        call print_test("Testing consistency across modes")
        
        ! Initialize test data with identical values across modes
        call random_number(arr_2d_in)
        do mode = 2, nmodes
            arr_2d_in(mode, :) = arr_2d_in(1, :)
        end do
        
        m = 2
        nheal = 1
        
        ! Process each mode
        do mode = 1, nmodes
            call s_to_rho_healaxis_2d(m, ns, nrho, nheal, mode, nmodes, arr_2d_in, arr_2d_out1)
        end do
        
        ! Verify all modes produce identical results (since input was identical)
        max_diff = 0.0_dp
        do mode = 2, nmodes
            do j = 0, nrho-1
                diff = abs(arr_2d_out1(1, j) - arr_2d_out1(mode, j))
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
    end subroutine test_consistency

    subroutine test_edge_cases()
        integer, parameter :: nmodes = 2
        real(dp), dimension(nmodes, 10) :: small_input
        real(dp), dimension(nmodes, 0:9) :: small_output
        real(dp), dimension(nmodes, 100) :: large_input
        real(dp), dimension(nmodes, 0:99) :: large_output
        
        call print_test("Testing edge cases")
        
        ! Test with small arrays (minimum size for spline order 5)
        call random_number(small_input)
        call s_to_rho_healaxis_2d(1, 10, 10, 2, 1, nmodes, small_input, small_output)
        
        ! Test with m = 0 (different code path)
        call s_to_rho_healaxis_2d(0, 10, 10, 1, 1, nmodes, small_input, small_output)
        
        ! Test with large arrays
        call random_number(large_input)
        call s_to_rho_healaxis_2d(3, 100, 100, 4, 1, nmodes, large_input, large_output)
        
        ! Test with moderate nheal value
        call s_to_rho_healaxis_2d(1, 10, 10, 3, 1, nmodes, small_input, small_output)
        
        ! If we get here without crashes, edge cases pass
        call print_ok()
    end subroutine test_edge_cases

end program test_spline_vmec_array_temps