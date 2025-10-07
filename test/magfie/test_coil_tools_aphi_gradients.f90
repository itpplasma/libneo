program test_coil_tools_aphi_gradients
    use, intrinsic :: iso_fortran_env, only: output_unit
    use coil_tools, only: coil_t, coil_init, coil_deinit
    use libneo_kinds, only: dp
    use math_constants, only: pi
    use util_for_test, only: print_fail, print_ok, print_test

    implicit none

    real(dp), parameter :: fd_step = 1.0e-6_dp
    real(dp), parameter :: abs_tol = 1.0e-10_dp
    real(dp), parameter :: rel_tol = 1.0e-4_dp
    real(dp), parameter :: min_distance = 1.0e-6_dp
    real(dp), parameter :: max_eccentricity = 0.999_dp
    real(dp), parameter :: coil_radius = 1.0_dp
    real(dp), parameter :: exclusion_R = 0.08_dp
    integer, parameter :: n_test_points = 20
    character(len=*), parameter :: plot_directory = COIL_TOOLS_PLOT_DIR

    type(coil_t) :: coil
    real(dp) :: test_R(n_test_points), test_Z(n_test_points)
    real(dp) :: R, Z, phi
    real(dp) :: Aphi_center, Aphi_R_plus, Aphi_R_minus, Aphi_Z_plus, Aphi_Z_minus
    real(dp) :: dAphi_dR_analytic, dAphi_dZ_analytic
    real(dp) :: dAphi_dR_fd, dAphi_dZ_fd
    real(dp) :: error_R(n_test_points), error_Z(n_test_points)
    real(dp) :: max_error_R, max_error_Z
    integer :: k
    logical :: test_passed

    call print_test('coil_tools Aphi gradients (point-by-point finite difference)')

    ! Create circular coil at R=1.0
    call coil_init(coil, 128, 1)
    call create_circular_coil(coil, coil_radius, 128)

    ! Generate test points in R, avoiding near-coil region
    do k = 1, n_test_points
        test_R(k) = 0.5_dp + 1.1_dp * real(k - 1, dp) / real(n_test_points - 1, dp)
        test_Z(k) = -0.5_dp + 1.0_dp * real(k - 1, dp) / real(n_test_points - 1, dp)
    end do

    ! Test at phi = 0 for simplicity
    phi = 0.0_dp
    max_error_R = 0.0_dp
    max_error_Z = 0.0_dp

    ! Test dAphi/dR at each point
    do k = 1, n_test_points
        R = test_R(k)
        Z = 0.0_dp

        if (abs(R - coil_radius) < exclusion_R) cycle

        ! Compute at center and perturbed R
        call compute_Aphi_and_gradient(R, Z, phi, Aphi_center, dAphi_dR_analytic, dAphi_dZ_analytic)
        call compute_Aphi_and_gradient(R + fd_step, Z, phi, Aphi_R_plus)
        call compute_Aphi_and_gradient(R - fd_step, Z, phi, Aphi_R_minus)

        dAphi_dR_fd = (Aphi_R_plus - Aphi_R_minus) / (2.0_dp * fd_step)
        error_R(k) = relative_error(dAphi_dR_fd, dAphi_dR_analytic, abs_tol)
        max_error_R = max(max_error_R, error_R(k))
    end do

    ! Test dAphi/dZ at each point
    do k = 1, n_test_points
        R = 1.0_dp
        Z = test_Z(k)

        call compute_Aphi_and_gradient(R, Z, phi, Aphi_center, dAphi_dR_analytic, dAphi_dZ_analytic)
        call compute_Aphi_and_gradient(R, Z + fd_step, phi, Aphi_Z_plus)
        call compute_Aphi_and_gradient(R, Z - fd_step, phi, Aphi_Z_minus)

        dAphi_dZ_fd = (Aphi_Z_plus - Aphi_Z_minus) / (2.0_dp * fd_step)
        error_Z(k) = relative_error(dAphi_dZ_fd, dAphi_dZ_analytic, abs_tol)
        max_error_Z = max(max_error_Z, error_Z(k))
    end do

    ! Save and plot results
    call ensure_plot_directory(plot_directory)
    call save_error_data(test_R, error_R, trim(plot_directory)//'/aphi_dR_error_pointwise.csv')
    call save_error_data(test_Z, error_Z, trim(plot_directory)//'/aphi_dZ_error_pointwise.csv')
    call generate_plots(plot_directory)

    ! Check tolerance
    test_passed = (max_error_R <= rel_tol) .and. (max_error_Z <= rel_tol)

    if (test_passed) then
        call print_ok
    else
        call print_fail
        write (output_unit, '(a, es12.5)') 'max rel err dAphi/dR:', max_error_R
        write (output_unit, '(a, es12.5)') 'max rel err dAphi/dZ:', max_error_Z
        call coil_deinit(coil)
        error stop 'Aphi gradient validation failed'
    end if

    call coil_deinit(coil)

contains

    subroutine create_circular_coil(coil, radius, nseg)
        type(coil_t), intent(inout) :: coil
        real(dp), intent(in) :: radius
        integer, intent(in) :: nseg
        integer :: k
        real(dp) :: theta

        do k = 1, nseg
            theta = 2.0_dp * pi * real(k - 1, dp) / real(nseg, dp)
            coil%XYZ(1, k) = radius * cos(theta)
            coil%XYZ(2, k) = radius * sin(theta)
            coil%XYZ(3, k) = 0.0_dp
        end do
    end subroutine create_circular_coil

    subroutine compute_Aphi_and_gradient(R, Z, phi, Aphi, dAphi_dR, dAphi_dZ)
        use coil_tools, only: segment_eccentricity, segment_vector_potential_contribution, &
                              segment_gradient_kernel, segment_gradient_contribution
        real(dp), intent(in) :: R, Z, phi
        real(dp), intent(out) :: Aphi
        real(dp), intent(out), optional :: dAphi_dR, dAphi_dZ
        real(dp) :: XYZ_r(3), XYZ_i(3), XYZ_f(3), XYZ_if(3)
        real(dp) :: AXYZ(3), grad_AX(3), grad_AY(3), grad_AZ(3)
        real(dp) :: dist_i, dist_f, dist_if, eccentricity
        real(dp) :: common_gradient_term(3)
        real(dp) :: cosphi, sinphi
        integer :: ks, ks_prev

        cosphi = cos(phi)
        sinphi = sin(phi)
        XYZ_r = [R * cosphi, R * sinphi, Z]

        AXYZ = 0.0_dp
        grad_AX = 0.0_dp
        grad_AY = 0.0_dp
        grad_AZ = 0.0_dp

        ks_prev = coil%nseg
        XYZ_f = coil%XYZ(:, ks_prev) - XYZ_r
        dist_f = max(min_distance, sqrt(sum(XYZ_f * XYZ_f)))

        do ks = 1, coil%nseg
            XYZ_i = XYZ_f
            dist_i = dist_f
            XYZ_f = coil%XYZ(:, ks) - XYZ_r
            dist_f = max(min_distance, sqrt(sum(XYZ_f * XYZ_f)))
            XYZ_if = coil%XYZ(:, ks) - coil%XYZ(:, ks_prev)
            dist_if = sqrt(sum(XYZ_if * XYZ_if))

            eccentricity = segment_eccentricity(dist_if, dist_i, dist_f, max_eccentricity)
            AXYZ = AXYZ + segment_vector_potential_contribution(XYZ_if, dist_if, eccentricity)

            if (present(dAphi_dR) .or. present(dAphi_dZ)) then
                common_gradient_term = segment_gradient_kernel(XYZ_i, XYZ_f, dist_i, dist_f)
                call segment_gradient_contribution(XYZ_if, common_gradient_term, grad_AX, grad_AY, grad_AZ)
            end if

            ks_prev = ks
        end do

        Aphi = AXYZ(2) * cosphi - AXYZ(1) * sinphi

        if (present(dAphi_dR)) then
            dAphi_dR = grad_AY(1) * cosphi**2 - grad_AX(2) * sinphi**2 + &
                       (grad_AY(2) - grad_AX(1)) * cosphi * sinphi
        end if

        if (present(dAphi_dZ)) then
            dAphi_dZ = grad_AY(3) * cosphi - grad_AX(3) * sinphi
        end if
    end subroutine compute_Aphi_and_gradient

    subroutine compute_Aphi_and_gradient_interface(R, Z, phi, Aphi)
        real(dp), intent(in) :: R, Z, phi
        real(dp), intent(out) :: Aphi
        call compute_Aphi_and_gradient(R, Z, phi, Aphi)
    end subroutine compute_Aphi_and_gradient_interface

    pure function relative_error(fd_val, analytic_val, abs_floor) result(err)
        real(dp), intent(in) :: fd_val, analytic_val, abs_floor
        real(dp) :: err
        real(dp) :: denom
        denom = max(abs(analytic_val), abs_floor)
        err = abs(fd_val - analytic_val) / denom
    end function relative_error

    subroutine ensure_plot_directory(path)
        character(len=*), intent(in) :: path
        integer :: istat
        call execute_command_line('mkdir -p "'//trim(path)//'"', exitstat=istat)
        if (istat /= 0) then
            call print_fail
            write (output_unit, '(a)') 'failed to create plot directory: '//trim(path)
            error stop 'failed to create plot directory'
        end if
    end subroutine ensure_plot_directory

    subroutine save_error_data(axis_vals, error_vals, filename)
        real(dp), intent(in) :: axis_vals(:), error_vals(:)
        character(len=*), intent(in) :: filename
        integer :: fid, k
        open (newunit=fid, file=filename, status='replace', action='write')
        write (fid, '(a)') 'coord,error'
        do k = 1, size(axis_vals)
            write (fid, '(es24.16e3,",",es24.16e3)') axis_vals(k), error_vals(k)
        end do
        close (fid)
    end subroutine save_error_data

    subroutine generate_plots(plot_dir)
        character(len=*), intent(in) :: plot_dir
        character(len=1024) :: command
        integer :: istat
        command = 'python3 "'//trim(COIL_TOOLS_SCRIPT_DIR)// &
                  '/plot_aphi_gradients.py" "'//trim(plot_dir)//'"'
        call execute_command_line(trim(command), exitstat=istat)
        if (istat /= 0) then
            write (output_unit, '(a)') 'Warning: failed to generate plots'
        end if
    end subroutine generate_plots

end program test_coil_tools_aphi_gradients
