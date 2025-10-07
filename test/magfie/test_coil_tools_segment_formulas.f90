program test_coil_tools_segment_formulas
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    use coil_tools, only: coil_t, coil_init, coil_deinit, &
                          segment_vector_potential_contribution, &
                          segment_eccentricity, &
                          segment_gradient_kernel, &
                          segment_gradient_contribution
    use libneo_kinds, only: dp
    use math_constants, only: pi
    use util_for_test, only: print_fail, print_ok, print_test

    implicit none

    real(dp), parameter :: abs_tol = 1.0e-10_dp
    real(dp), parameter :: rel_tol = 1.0e-4_dp
    real(dp), parameter :: fd_delta = 1.0e-7_dp
    real(dp), parameter :: max_eccentricity = 0.999_dp
    character(len=*), parameter :: plot_directory = COIL_TOOLS_PLOT_DIR

    call print_test('coil_tools Biot-Savart segment formula validation')

    call test_segment_eccentricity_formula()
    call test_vector_potential_contribution_formula()
    call test_gradient_kernel_formula()
    call test_gradient_contribution_formula()

    call print_ok

contains

    subroutine test_segment_eccentricity_formula()
        real(dp) :: dist_seg, dist_i, dist_f, ecc
        write (output_unit, '(a)') '  Testing segment_eccentricity formula...'

        dist_seg = 0.1_dp
        dist_i = 0.5_dp
        dist_f = 0.6_dp
        ecc = segment_eccentricity(dist_seg, dist_i, dist_f, max_eccentricity)

        ! Expected: min(0.999, 0.1/(0.5+0.6)) = min(0.999, 0.0909...) = 0.0909...
        if (abs(ecc - 0.1_dp/(dist_i + dist_f)) > abs_tol) then
            write (error_unit, '(a,es12.4,a,es12.4)') '  ERROR: ecc=', ecc, &
                ' expected=', 0.1_dp/(dist_i + dist_f)
            error stop
        end if
        write (output_unit, '(a)') '    segment_eccentricity: OK'
    end subroutine test_segment_eccentricity_formula

    subroutine test_vector_potential_contribution_formula()
        write (output_unit, '(a)') '  Testing vector_potential_contribution formula...'
        ! This will be validated with finite differences in integration test
        write (output_unit, '(a)') '    (deferred to integration test)'
    end subroutine test_vector_potential_contribution_formula

    subroutine test_gradient_kernel_formula()
        write (output_unit, '(a)') '  Testing gradient_kernel formula...'
        ! This will be validated with finite differences
        write (output_unit, '(a)') '    (deferred to finite difference validation)'
    end subroutine test_gradient_kernel_formula

    subroutine test_gradient_contribution_formula()
        real(dp) :: XYZ_segment(3), grad_kernel(3)
        real(dp) :: grad_AX(3), grad_AY(3), grad_AZ(3)

        write (output_unit, '(a)') '  Testing gradient_contribution formula...'

        XYZ_segment = [1.0_dp, 2.0_dp, 3.0_dp]
        grad_kernel = [0.1_dp, 0.2_dp, 0.3_dp]
        grad_AX = 0.0_dp
        grad_AY = 0.0_dp
        grad_AZ = 0.0_dp

        call segment_gradient_contribution(XYZ_segment, grad_kernel, grad_AX, grad_AY, grad_AZ)

        ! Expected:
        ! grad_AX = XYZ_segment(1) * grad_kernel = 1.0 * [0.1, 0.2, 0.3]
        ! grad_AY = XYZ_segment(2) * grad_kernel = 2.0 * [0.1, 0.2, 0.3]
        ! grad_AZ = XYZ_segment(3) * grad_kernel = 3.0 * [0.1, 0.2, 0.3]

        if (any(abs(grad_AX - XYZ_segment(1) * grad_kernel) > abs_tol)) then
            write (error_unit, '(a)') '  ERROR: grad_AX mismatch'
            error stop
        end if
        if (any(abs(grad_AY - XYZ_segment(2) * grad_kernel) > abs_tol)) then
            write (error_unit, '(a)') '  ERROR: grad_AY mismatch'
            error stop
        end if
        if (any(abs(grad_AZ - XYZ_segment(3) * grad_kernel) > abs_tol)) then
            write (error_unit, '(a)') '  ERROR: grad_AZ mismatch'
            error stop
        end if

        write (output_unit, '(a)') '    gradient_contribution: OK'
    end subroutine test_gradient_contribution_formula

end program test_coil_tools_segment_formulas
