program test_coil_tools_vector_potential_derivs
    use, intrinsic :: iso_fortran_env, only: output_unit
    use coil_tools, only: coil_t, coil_init, coil_deinit, grid_from_bounding_box, &
                          vector_potential_biot_savart_fourier
    use libneo_kinds, only: dp
    use math_constants, only: pi
    use util_for_test, only: print_fail, print_ok, print_test

    implicit none

    integer, parameter :: nphi = 128
    integer, parameter :: nmax = 16
    integer, parameter :: nR = 101
    integer, parameter :: nZ = 41
    real(dp), parameter :: Rmin = 0.505_dp
    real(dp), parameter :: Rmax = 1.505_dp
    real(dp), parameter :: Zmin = -0.6_dp
    real(dp), parameter :: Zmax = 0.6_dp
    real(dp), parameter :: min_distance = 1.0e-6_dp
    real(dp), parameter :: max_eccentricity = 0.999_dp
    real(dp), parameter :: abs_tol = 1.0e-8_dp
    real(dp), parameter :: rel_tol = 1.0e-1_dp
    logical, parameter :: use_convex_wall = .false.
    character(len=*), parameter :: plot_directory = COIL_TOOLS_PLOT_DIR

    type(coil_t), allocatable :: coils(:)
    complex(dp), allocatable :: AnR(:, :, :, :), Anphi(:, :, :, :), AnZ(:, :, :, :)
    complex(dp), allocatable :: dAnphi_dR(:, :, :, :), dAnphi_dZ(:, :, :, :)
    complex(dp), allocatable :: dAnphi_dx(:, :, :, :), dAnphi_dy(:, :, :, :)
    complex(dp), allocatable :: fd_dAnphi_dR(:, :, :, :), fd_dAnphi_dZ(:, :, :, :)
    complex(dp) :: fd_value
    real(dp), allocatable :: R(:), Z(:), phi(:)
    real(dp), allocatable :: cosphi(:), sinphi(:)
    real(dp) :: rel_err, inv_R
    real(dp) :: max_rel_err_R, max_rel_err_Z, max_rel_err_dx, max_rel_err_dy
    real(dp), allocatable :: dphi_error_vs_R(:), dphi_error_vs_Z(:)
    real(dp), allocatable :: dphi_error_dx_vs_R(:), dphi_error_dy_vs_R(:)
    real(dp), allocatable :: Aphi_samples(:, :), dAphi_dx_samples(:, :), dAphi_dy_samples(:, :)
    real(dp), allocatable :: dAphi_dx_fd(:, :), dAphi_dy_fd(:, :), dAphi_dR_fd(:, :), dAphi_dphi_fd(:, :)
    real(dp), allocatable :: tmp_phi_values(:)
    integer, allocatable :: best_mode_R(:)
    integer :: mid_R, mid_Z, mid_mode
    integer :: kc, kR, kZ, nmode, fid_debug, kphi
    complex(dp), allocatable :: tmp_complex(:)
    real(dp), parameter :: coil_radius = 1.0_dp
    real(dp), parameter :: exclusion_R = 0.08_dp
    logical :: skip_point

    call print_test('coil_tools vector potential derivative consistency '// &
                    '(analytic vs finite diff)')

    allocate (coils(1))
    call create_circular_coil(coils(1), 1.0_dp, 128)

    allocate (R(nR), Z(nZ), phi(nphi))
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)
    allocate (cosphi(nphi), sinphi(nphi))
    cosphi = cos(phi)
    sinphi = sin(phi)

    call vector_potential_biot_savart_fourier( &
        coils, nmax, min_distance, max_eccentricity, use_convex_wall, &
        Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ, &
        dAnphi_dx, dAnphi_dy)

    allocate (fd_dAnphi_dR(0:nmax, nR, nZ, size(coils)))
    allocate (fd_dAnphi_dZ(0:nmax, nR, nZ, size(coils)))
    fd_dAnphi_dR = (0.0_dp, 0.0_dp)
    fd_dAnphi_dZ = (0.0_dp, 0.0_dp)

    max_rel_err_R = 0.0_dp
    max_rel_err_Z = 0.0_dp
    max_rel_err_dx = 0.0_dp
    max_rel_err_dy = 0.0_dp
    mid_R = (nR + 1)/2
    mid_Z = (nZ + 1)/2
    mid_mode = 2
    allocate (dphi_error_vs_R(nR))
    allocate (best_mode_R(nR))
    allocate (dphi_error_vs_Z(nZ))
    allocate (dphi_error_dx_vs_R(nR))
    allocate (dphi_error_dy_vs_R(nR))
    allocate (Aphi_samples(nR, nphi))
    allocate (dAphi_dx_samples(nR, nphi))
    allocate (dAphi_dy_samples(nR, nphi))
    allocate (dAphi_dx_fd(nR, nphi))
    allocate (dAphi_dy_fd(nR, nphi))
    allocate (dAphi_dR_fd(nR, nphi))
    allocate (dAphi_dphi_fd(nR, nphi))
    allocate (tmp_phi_values(nphi))
    allocate (tmp_complex(nR))
    dphi_error_vs_R = 0.0_dp
    best_mode_R = -1
    dphi_error_vs_Z = 0.0_dp
    dphi_error_dx_vs_R = 0.0_dp
    dphi_error_dy_vs_R = 0.0_dp

    do kc = 1, size(coils)
        do nmode = 0, nmax
            do kZ = 1, nZ
                do kR = 1, nR
                    fd_value = finite_difference_R(Anphi(nmode, :, kZ, kc), R, kR)
                    fd_dAnphi_dR(nmode, kR, kZ, kc) = fd_value
                    skip_point = abs(R(kR) - coil_radius) < exclusion_R
                    if (.not. skip_point) then
                        rel_err = derivative_error(fd_value, dAnphi_dR(nmode, kR, kZ, kc), abs_tol)
                        max_rel_err_R = max(max_rel_err_R, rel_err)
                        if (kZ == mid_Z) then
                            if (rel_err > dphi_error_vs_R(kR)) then
                                dphi_error_vs_R(kR) = rel_err
                                best_mode_R(kR) = nmode
                            end if
                        end if
                    end if
                end do
            end do
            do kR = 1, nR
                do kZ = 1, nZ
                    fd_value = finite_difference_Z(Anphi(nmode, kR, :, kc), Z, kZ)
                    fd_dAnphi_dZ(nmode, kR, kZ, kc) = fd_value
                    skip_point = abs(R(kR) - coil_radius) < exclusion_R
                    if (.not. skip_point) then
                        rel_err = derivative_error(fd_value, dAnphi_dZ(nmode, kR, kZ, kc), abs_tol)
                        max_rel_err_Z = max(max_rel_err_Z, rel_err)
                        if (kR == mid_R) then
                            dphi_error_vs_Z(kZ) = max(dphi_error_vs_Z(kZ), rel_err)
                        end if
                    end if
                end do
            end do
        end do
    end do

    call ensure_plot_directory(plot_directory)
    ! reconstruct physical-space samples for additional diagnostics (mid-plane of coil 1)
    do kR = 1, nR
        call reconstruct_real_series(Anphi(:, kR, mid_Z, 1), phi, tmp_phi_values)
        Aphi_samples(kR, :) = tmp_phi_values
        call reconstruct_real_series(dAnphi_dx(:, kR, mid_Z, 1), phi, tmp_phi_values)
        dAphi_dx_samples(kR, :) = tmp_phi_values
        call reconstruct_real_series(dAnphi_dy(:, kR, mid_Z, 1), phi, tmp_phi_values)
        dAphi_dy_samples(kR, :) = tmp_phi_values
    end do

    do kphi = 1, nphi
        tmp_complex = cmplx(Aphi_samples(:, kphi), 0.0_dp)
        do kR = 1, nR
            dAphi_dR_fd(kR, kphi) = real(finite_difference_R(tmp_complex, R, kR))
        end do
    end do

    do kR = 1, nR
        do kphi = 1, nphi
            dAphi_dphi_fd(kR, kphi) = finite_difference_phi(Aphi_samples(kR, :), kphi, phi(2) - phi(1))
        end do
    end do

    do kR = 1, nR
        if (R(kR) <= 0.0_dp) then
            inv_R = 0.0_dp
        else
            inv_R = 1.0_dp/R(kR)
        end if
        do kphi = 1, nphi
            dAphi_dx_fd(kR, kphi) = cosphi(kphi) * dAphi_dR_fd(kR, kphi) - &
                sinphi(kphi) * inv_R * dAphi_dphi_fd(kR, kphi)
            dAphi_dy_fd(kR, kphi) = sinphi(kphi) * dAphi_dR_fd(kR, kphi) + &
                cosphi(kphi) * inv_R * dAphi_dphi_fd(kR, kphi)
        end do
    end do

    do kR = 1, nR
        do kphi = 1, nphi
            rel_err = abs(dAphi_dx_fd(kR, kphi) - dAphi_dx_samples(kR, kphi))/ &
                max(abs(dAphi_dx_samples(kR, kphi)), abs_tol)
            max_rel_err_dx = max(max_rel_err_dx, rel_err)
            dphi_error_dx_vs_R(kR) = max(dphi_error_dx_vs_R(kR), rel_err)

            rel_err = abs(dAphi_dy_fd(kR, kphi) - dAphi_dy_samples(kR, kphi))/ &
                max(abs(dAphi_dy_samples(kR, kphi)), abs_tol)
            max_rel_err_dy = max(max_rel_err_dy, rel_err)
            dphi_error_dy_vs_R(kR) = max(dphi_error_dy_vs_R(kR), rel_err)
        end do
    end do

    call write_error_series(trim(plot_directory)//'/coil_tools_dAphi_dR_error.dat', &
        R, dphi_error_vs_R)
    call write_error_series(trim(plot_directory)//'/coil_tools_dAphi_dZ_error.dat', &
        Z, dphi_error_vs_Z)
    call write_error_series(trim(plot_directory)//'/coil_tools_dAphi_dx_error.dat', &
        R, dphi_error_dx_vs_R)
    call write_error_series(trim(plot_directory)//'/coil_tools_dAphi_dy_error.dat', &
        R, dphi_error_dy_vs_R)
    ! store maximum errors across slices for diagnostics if needed
    open(newunit=fid_debug, file=trim(plot_directory)//'/coil_tools_dAphi_debug.dat', status='replace', action='write')
    write(fid_debug, '(a)') '# R  mode  analytic_real  analytic_imag  fd_real  fd_imag  rel_err'
    do kR = 1, nR
        mid_mode = best_mode_R(kR)
        fd_value = fd_dAnphi_dR(mid_mode, kR, mid_Z, 1)
        write(fid_debug, '(es24.16,1x,i4,1x,5es24.16)') R(kR), mid_mode, &
            real(dAnphi_dR(mid_mode, kR, mid_Z, 1)), aimag(dAnphi_dR(mid_mode, kR, mid_Z, 1)), &
            real(fd_value), aimag(fd_value), &
            dphi_error_vs_R(kR)
        if (kR <= 5) then
            write (output_unit, '(a,1x,i0,1x,3es12.5)') 'debug mode/analytic/fd', mid_mode, &
                real(dAnphi_dR(mid_mode, kR, mid_Z, 1)), real(fd_value), &
                real(Anphi(mid_mode, kR, mid_Z, 1))
        end if
    end do
    close(fid_debug)

    if (max_rel_err_R > rel_tol .or. max_rel_err_Z > rel_tol) then
        call print_fail
        write (output_unit, '(a, 1x, es12.5)') 'max rel err dAphi/dR:', max_rel_err_R
        write (output_unit, '(a, 1x, es12.5)') 'max rel err dAphi/dZ:', max_rel_err_Z
        call coil_deinit(coils(1))
        error stop 'Vector potential derivative comparison failed beyond tolerance'
    end if

    call print_ok
    write (output_unit, '(a, 1x, es12.5)') 'diag max rel err dAphi/dx:', max_rel_err_dx
    write (output_unit, '(a, 1x, es12.5)') 'diag max rel err dAphi/dy:', max_rel_err_dy
    call coil_deinit(coils(1))

contains

    subroutine create_circular_coil(coil, radius, nseg)
        type(coil_t), intent(inout) :: coil
        real(dp), intent(in) :: radius
        integer, intent(in) :: nseg
        integer :: k
        real(dp) :: theta

        call coil_init(coil, nseg, 1)
        do k = 1, nseg
            theta = 2.0_dp*pi*real(k - 1, dp)/real(nseg, dp)
            coil%XYZ(1, k) = radius*cos(theta)
            coil%XYZ(2, k) = radius*sin(theta)
            coil%XYZ(3, k) = 0.0_dp
    end do
    end subroutine create_circular_coil

    subroutine reconstruct_real_series(coeffs, phi_vals, output)
        complex(dp), intent(in) :: coeffs(0:)
        real(dp), intent(in) :: phi_vals(:)
        real(dp), intent(out) :: output(:)
        integer :: nmode, kphi

        if (size(output) /= size(phi_vals)) then
            call print_fail
            write (output_unit, '(a)') 'reconstruct_real_series: array size mismatch'
            error stop 'reconstruct_real_series: array size mismatch'
        end if

        output(:) = real(coeffs(0))
        do nmode = 1, ubound(coeffs, 1)
            do kphi = 1, size(phi_vals)
                output(kphi) = output(kphi) + 2.0_dp * ( &
                    real(coeffs(nmode)) * cos(real(nmode, dp) * phi_vals(kphi)) - &
                    aimag(coeffs(nmode)) * sin(real(nmode, dp) * phi_vals(kphi)) )
            end do
        end do
    end subroutine reconstruct_real_series

    pure function finite_difference_R(values, coords, idx) result(df)
        complex(dp), intent(in) :: values(:)
        real(dp), intent(in) :: coords(:)
        integer, intent(in) :: idx
        complex(dp) :: df
        real(dp) :: delta

        if (idx == 1) then
            delta = coords(2) - coords(1)
            df = (values(2) - values(1))/delta
        else if (idx == size(coords)) then
            delta = coords(idx) - coords(idx - 1)
            df = (values(idx) - values(idx - 1))/delta
        else
            delta = coords(idx + 1) - coords(idx - 1)
            df = (values(idx + 1) - values(idx - 1))/delta
        end if
    end function finite_difference_R

    pure function finite_difference_Z(values, coords, idx) result(df)
        complex(dp), intent(in) :: values(:)
        real(dp), intent(in) :: coords(:)
        integer, intent(in) :: idx
        complex(dp) :: df
        real(dp) :: delta

        if (idx == 1) then
            delta = coords(2) - coords(1)
            df = (values(2) - values(1))/delta
        else if (idx == size(coords)) then
            delta = coords(idx) - coords(idx - 1)
            df = (values(idx) - values(idx - 1))/delta
        else
            delta = coords(idx + 1) - coords(idx - 1)
            df = (values(idx + 1) - values(idx - 1))/delta
        end if
    end function finite_difference_Z

    pure function finite_difference_phi(values, idx, dphi) result(df)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: idx
        real(dp), intent(in) :: dphi
        real(dp) :: df
        integer :: npts, prev_idx, next_idx

        npts = size(values)
        if (npts < 2) then
            df = 0.0_dp
            return
        end if

        if (idx == 1) then
            prev_idx = npts
        else
            prev_idx = idx - 1
        end if

        if (idx == npts) then
            next_idx = 1
        else
            next_idx = idx + 1
        end if

        df = (values(next_idx) - values(prev_idx))/(2.0_dp * dphi)
    end function finite_difference_phi

    pure function derivative_error(fd_val, analytic_val, abs_floor) result(err)
        complex(dp), intent(in) :: fd_val
        complex(dp), intent(in) :: analytic_val
        real(dp), intent(in) :: abs_floor
        real(dp) :: err
        real(dp) :: denom

        denom = max(abs(analytic_val), abs_floor)
        err = abs(fd_val - analytic_val)/denom
    end function derivative_error

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

    subroutine write_error_series(filename, axis_vals, error_vals)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: axis_vals(:)
        real(dp), intent(in) :: error_vals(:)
        integer :: fid, idx

        if (size(axis_vals) /= size(error_vals)) then
            call print_fail
            write (output_unit, '(a)') 'write_error_series: array size mismatch'
            error stop 'write_error_series: array size mismatch'
        end if

        open(newunit=fid, file=filename, status='replace', action='write')
        write(fid, '(a)') '# coord  relative_error'
        do idx = 1, size(axis_vals)
            write(fid, '(es24.16,1x,es24.16)') axis_vals(idx), error_vals(idx)
        end do
        close(fid)
    end subroutine write_error_series

end program test_coil_tools_vector_potential_derivs
