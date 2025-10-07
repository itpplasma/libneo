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
    real(dp), parameter :: segment_start(3) = [-0.35_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: segment_end(3) = [0.35_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: observation_point(3) = [0.12_dp, 0.27_dp, -0.18_dp]
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
        real(dp) :: obs(3)
        real(dp) :: XYZ_segment(3)
        real(dp) :: dist_segment, dist_i, dist_f, ecc
        real(dp) :: contribution(3), expected(3)
        real(dp) :: expected_prefactor, tol

        write (output_unit, '(a)') '  Testing vector_potential_contribution formula...'

        obs = observation_point
        XYZ_segment = segment_end - segment_start
        dist_segment = sqrt(sum(XYZ_segment ** 2))
        dist_i = sqrt(sum((segment_start - obs) ** 2))
        dist_f = sqrt(sum((segment_end - obs) ** 2))
        ecc = segment_eccentricity(dist_segment, dist_i, dist_f, max_eccentricity)
        contribution = segment_vector_potential_contribution(XYZ_segment, dist_segment, ecc)

        expected_prefactor = log((dist_i + dist_f + dist_segment) / &
          (dist_i + dist_f - dist_segment))
        expected = XYZ_segment / dist_segment * expected_prefactor

        tol = abs_tol + rel_tol * maxval(abs(expected))
        if (maxval(abs(contribution - expected)) > tol) then
            write (error_unit, '(a)') '  ERROR: segment vector potential formula mismatch'
            write (error_unit, '(a,3es13.5)') '    contribution = ', contribution
            write (error_unit, '(a,3es13.5)') '    expected      = ', expected
            error stop
        end if

        write (output_unit, '(a)') '    segment_vector_potential_contribution: OK'
    end subroutine test_vector_potential_contribution_formula

    subroutine test_gradient_kernel_formula()
        real(dp) :: obs(3), obs_plus(3), obs_minus(3)
        real(dp) :: grad_AX(3), grad_AY(3), grad_AZ(3)
        real(dp) :: fd_grad_AX(3), fd_grad_AY(3), fd_grad_AZ(3)
        real(dp) :: A_plus(3), A_minus(3)
        real(dp) :: tol_AX, tol_AY, tol_AZ
        integer :: idx

        write (output_unit, '(a)') '  Testing gradient_kernel formula...'

        obs = observation_point
        call compute_gradient(obs, grad_AX, grad_AY, grad_AZ)

        do idx = 1, 3
            obs_plus = obs
            obs_minus = obs
            obs_plus(idx) = obs_plus(idx) + fd_delta
            obs_minus(idx) = obs_minus(idx) - fd_delta
            A_plus = compute_vector_potential(obs_plus)
            A_minus = compute_vector_potential(obs_minus)
            fd_grad_AX(idx) = (A_plus(1) - A_minus(1)) / (2.0_dp * fd_delta)
            fd_grad_AY(idx) = (A_plus(2) - A_minus(2)) / (2.0_dp * fd_delta)
            fd_grad_AZ(idx) = (A_plus(3) - A_minus(3)) / (2.0_dp * fd_delta)
        end do

        tol_AX = abs_tol + rel_tol * max(1.0_dp, maxval(abs(fd_grad_AX)))
        tol_AY = abs_tol + rel_tol * max(1.0_dp, maxval(abs(fd_grad_AY)))
        tol_AZ = abs_tol + rel_tol * max(1.0_dp, maxval(abs(fd_grad_AZ)))

        if (maxval(abs(grad_AX - fd_grad_AX)) > tol_AX) then
            write (error_unit, '(a)') '  ERROR: grad_AX mismatch against finite difference benchmark'
            error stop
        end if
        if (maxval(abs(grad_AY - fd_grad_AY)) > tol_AY) then
            write (error_unit, '(a)') '  ERROR: grad_AY mismatch against finite difference benchmark'
            error stop
        end if
        if (maxval(abs(grad_AZ - fd_grad_AZ)) > tol_AZ) then
            write (error_unit, '(a)') '  ERROR: grad_AZ mismatch against finite difference benchmark'
            error stop
        end if

        write (output_unit, '(a)') '    segment_gradient_kernel: OK'
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

    pure function compute_vector_potential(obs) result(A_segment)
        real(dp), intent(in) :: obs(3)
        real(dp) :: A_segment(3)
        real(dp) :: XYZ_segment(3)
        real(dp) :: dist_segment, dist_i, dist_f, ecc
        real(dp) :: XYZ_i(3), XYZ_f(3)

        XYZ_segment = segment_end - segment_start
        dist_segment = sqrt(sum(XYZ_segment ** 2))
        XYZ_i = segment_start - obs
        dist_i = sqrt(sum(XYZ_i ** 2))
        XYZ_f = segment_end - obs
        dist_f = sqrt(sum(XYZ_f ** 2))
        ecc = segment_eccentricity(dist_segment, dist_i, dist_f, max_eccentricity)
        A_segment = segment_vector_potential_contribution(XYZ_segment, dist_segment, ecc)
    end function compute_vector_potential

    pure subroutine compute_gradient(obs, grad_AX, grad_AY, grad_AZ)
        real(dp), intent(in) :: obs(3)
        real(dp), intent(out) :: grad_AX(3), grad_AY(3), grad_AZ(3)
        real(dp) :: XYZ_segment(3)
        real(dp) :: XYZ_i(3), XYZ_f(3)
        real(dp) :: dist_i, dist_f
        real(dp) :: grad_kernel(3)

        XYZ_segment = segment_end - segment_start
        XYZ_i = segment_start - obs
        dist_i = sqrt(sum(XYZ_i ** 2))
        XYZ_f = segment_end - obs
        dist_f = sqrt(sum(XYZ_f ** 2))
        grad_kernel = segment_gradient_kernel(XYZ_i, XYZ_f, dist_i, dist_f)
        grad_AX = 0.0_dp
        grad_AY = 0.0_dp
        grad_AZ = 0.0_dp
        call segment_gradient_contribution(XYZ_segment, grad_kernel, grad_AX, grad_AY, grad_AZ)
    end subroutine compute_gradient

end program test_coil_tools_segment_formulas
