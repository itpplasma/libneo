program test_boozer_chartmap_asymmetric_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap, only: export_boozer_chartmap
    use boozer_coordinates_mod, only: n_theta_B, n_phi_B, use_B_r
    use boozer_sub, only: get_boozer_coordinates, sfl_quality_config_t
    use netcdf, only: nf90_open, nf90_nowrite, nf90_inq_varid, nf90_get_var, &
                      nf90_close, nf90_noerr, nf90_strerror

    implicit none

    character(*), parameter :: wout = 'wout_updown_tokamak.nc'
    character(*), parameter :: chartmap = 'updown_axis.chartmap.nc'
    type(sfl_quality_config_t) :: quality
    integer :: file, var
    real(dp) :: zaxis_m(1)
    real(dp), allocatable :: axis_z_cm(:, :, :)

    call check(nf90_open(wout, nf90_nowrite, file))
    call check(nf90_inq_varid(file, 'zaxis_cc', var))
    call check(nf90_get_var(file, var, zaxis_m))
    call check(nf90_close(file))

    use_B_r = .false.
    quality%allow_poor_equilibrium = .true.
    quality%override_reason = 'coarse public tokamak fixture'
    call get_boozer_coordinates(wout, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3, &
                                quality=quality)
    call export_boozer_chartmap(chartmap)

    allocate (axis_z_cm(1, n_theta_B - 1, n_phi_B - 1))
    call check(nf90_open(chartmap, nf90_nowrite, file))
    call check(nf90_inq_varid(file, 'z', var))
    call check(nf90_get_var(file, var, axis_z_cm, count=shape(axis_z_cm)))
    call check(nf90_close(file))
    if (maxval(abs(axis_z_cm - 100.0_dp*zaxis_m(1))) > 1.0e-9_dp) &
        error stop 'Boozer chartmap did not preserve the asymmetric axis offset'

contains

    subroutine check(status)
        integer, intent(in) :: status

        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine check

end program test_boozer_chartmap_asymmetric_axis
