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
    real(dp), parameter :: abs_tol = 1.0e-10_dp
    real(dp), parameter :: rel_tol = 1.0e-2_dp
    logical, parameter :: use_convex_wall = .false.
    character(len=*), parameter :: plot_directory = COIL_TOOLS_PLOT_DIR
    character(len=*), parameter :: script_directory = COIL_TOOLS_SCRIPT_DIR

    type(coil_t), allocatable :: coils(:)
    complex(dp), allocatable :: AnR(:, :, :, :), Anphi(:, :, :, :), AnZ(:, :, :, :)
    complex(dp), allocatable :: dAnphi_dR(:, :, :, :), dAnphi_dZ(:, :, :, :)
    complex(dp), allocatable :: fd_dAnphi_dR(:, :, :, :), fd_dAnphi_dZ(:, :, :, :)
    complex(dp) :: fd_value
    real(dp), allocatable :: R(:), Z(:), phi(:)
    real(dp) :: rel_err
    real(dp) :: max_rel_err_R, max_rel_err_Z
    real(dp), allocatable :: dphi_error_vs_R(:), dphi_error_vs_Z(:)
    integer :: mid_R, mid_Z
    integer :: kc, kR, kZ, nmode
    real(dp), parameter :: coil_radius = 1.0_dp
    real(dp), parameter :: exclusion_R = 0.08_dp
    logical :: skip_point

    call print_test('coil_tools vector potential derivative consistency '// &
                    '(analytic vs finite diff)')

    allocate (coils(1))
    call create_circular_coil(coils(1), 1.0_dp, 128)

    allocate (R(nR), Z(nZ), phi(nphi))
    call grid_from_bounding_box(Rmin, Rmax, nR, R, Zmin, Zmax, nZ, Z, nphi, phi)

    call vector_potential_biot_savart_fourier( &
        coils, nmax, min_distance, max_eccentricity, use_convex_wall, &
        Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)

    allocate (fd_dAnphi_dR(0:nmax, nR, nZ, size(coils)))
    allocate (fd_dAnphi_dZ(0:nmax, nR, nZ, size(coils)))
    fd_dAnphi_dR = (0.0_dp, 0.0_dp)
    fd_dAnphi_dZ = (0.0_dp, 0.0_dp)

    max_rel_err_R = 0.0_dp
    max_rel_err_Z = 0.0_dp
    mid_R = (nR + 1)/2
    mid_Z = (nZ + 1)/2
    allocate (dphi_error_vs_R(nR))
    allocate (dphi_error_vs_Z(nZ))
    dphi_error_vs_R = 0.0_dp
    dphi_error_vs_Z = 0.0_dp

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
                            dphi_error_vs_R(kR) = max(dphi_error_vs_R(kR), rel_err)
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
    call save_error_data(R, dphi_error_vs_R, trim(plot_directory)//'/dAphi_dR_error.csv')
    call save_error_data(Z, dphi_error_vs_Z, trim(plot_directory)//'/dAphi_dZ_error.csv')

    call validate_cartesian_gradients(R, Z, plot_directory)

    call generate_plots(plot_directory)

    if (max_rel_err_R > rel_tol .or. max_rel_err_Z > rel_tol) then
        call print_fail
        write (output_unit, '(a, 1x, es12.5)') 'max rel err dAphi/dR:', max_rel_err_R
        write (output_unit, '(a, 1x, es12.5)') 'max rel err dAphi/dZ:', max_rel_err_Z
        call coil_deinit(coils(1))
        error stop 'Vector potential derivative comparison failed beyond tolerance'
    end if

    call print_ok
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

    subroutine save_error_data(axis_vals, error_vals, filename)
        real(dp), intent(in) :: axis_vals(:)
        real(dp), intent(in) :: error_vals(:)
        character(len=*), intent(in) :: filename
        integer :: fid, k

        open(newunit=fid, file=filename, status='replace', action='write')
        write(fid, '(a)') 'coord,error'
        do k = 1, size(axis_vals)
            write(fid, '(es24.16e3, ",", es24.16e3)') axis_vals(k), error_vals(k)
        end do
        close(fid)
    end subroutine save_error_data

    subroutine generate_plots(plot_dir)
        character(len=*), intent(in) :: plot_dir
        character(len=1024) :: command
        integer :: istat

        command = 'python3 "'//trim(script_directory)// &
                  '/plot_vector_potential_derivs.py" "'//trim(plot_dir)//'"'
        call execute_command_line(trim(command), exitstat=istat)
        if (istat /= 0) then
            write (output_unit, '(a)') 'Warning: failed to generate plots with Python'
        end if
    end subroutine generate_plots

    subroutine validate_cartesian_gradients(R, Z, plot_dir)
        real(dp), intent(in) :: R(:), Z(:)
        character(len=*), intent(in) :: plot_dir
        real(dp), allocatable :: grad_AX(:, :, :, :), grad_AY(:, :, :, :), grad_AZ(:, :, :, :)
        real(dp), allocatable :: grad_profile(:)
        integer :: fid, kR, kZ, icomp

        allocate(grad_AX(3, nR, nZ, 1))
        allocate(grad_AY(3, nR, nZ, 1))
        allocate(grad_AZ(3, nR, nZ, 1))
        allocate(grad_profile(nR))

        open (newunit=fid, file='debug_cartesian_gradients.dat', status='old', action='read', form='unformatted')
        read (fid) grad_AX
        read (fid) grad_AY
        read (fid) grad_AZ
        close (fid)

        write (output_unit, '(a)') 'Cartesian gradients loaded successfully'
        write (output_unit, '(a, 3es12.4)') '  Sample grad_AX at midpoint:', grad_AX(:, nR/2, nZ/2, 1)
        write (output_unit, '(a, 3es12.4)') '  Sample grad_AY at midpoint:', grad_AY(:, nR/2, nZ/2, 1)
        write (output_unit, '(a, 3es12.4)') '  Sample grad_AZ at midpoint:', grad_AZ(:, nR/2, nZ/2, 1)
        write (output_unit, '(a, 2es12.4)') '  grad_AX range:', minval(grad_AX), maxval(grad_AX)
        write (output_unit, '(a, 2es12.4)') '  grad_AY range:', minval(grad_AY), maxval(grad_AY)
        write (output_unit, '(a, 2es12.4)') '  grad_AZ range:', minval(grad_AZ), maxval(grad_AZ)

        ! Save grad_AX components along R at Z=0
        do icomp = 1, 3
            do kR = 1, nR
                grad_profile(kR) = grad_AX(icomp, kR, nZ/2, 1)
            end do
            call save_error_data(R, grad_profile, trim(plot_dir)//'/grad_AX_comp'//char(ichar('0')+icomp)//'_vs_R.csv')
        end do

        ! Save grad_AY components along R at Z=0
        do icomp = 1, 3
            do kR = 1, nR
                grad_profile(kR) = grad_AY(icomp, kR, nZ/2, 1)
            end do
            call save_error_data(R, grad_profile, trim(plot_dir)//'/grad_AY_comp'//char(ichar('0')+icomp)//'_vs_R.csv')
        end do

        ! Save grad_AZ components along R at Z=0
        do icomp = 1, 3
            do kR = 1, nR
                grad_profile(kR) = grad_AZ(icomp, kR, nZ/2, 1)
            end do
            call save_error_data(R, grad_profile, trim(plot_dir)//'/grad_AZ_comp'//char(ichar('0')+icomp)//'_vs_R.csv')
        end do
    end subroutine validate_cartesian_gradients

end program test_coil_tools_vector_potential_derivs
