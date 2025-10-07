program geoflux_coord_dump
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use math_constants, only : pi
    use geoflux_field, only : spline_geoflux_data
    use geoflux_coordinates, only : geoflux_to_cyl

    implicit none

    integer, parameter :: max_len = 512
    integer, parameter :: n_surface = 5
    integer, parameter :: n_surface_samples = 180
    integer, parameter :: n_theta_lines = 6
    integer, parameter :: n_radial_samples = 120
    real(dp), parameter :: two_pi = 2.0_dp * pi
    real(dp), parameter :: s_min = 5.0d-3
    real(dp), parameter :: s_max = 9.5d-1

    character(len=max_len) :: geqdsk_file
    character(len=max_len) :: output_file
    integer :: unit
    real(dp) :: surface_values(n_surface)
    real(dp) :: theta_values(n_theta_lines)
    real(dp) :: x_geo(3)
    real(dp) :: x_cyl(3)
    real(dp) :: theta, s_value
    integer :: i_surf, i_sample, i_theta

    call init_argument(1, geqdsk_file)
    call init_argument(2, output_file)

    call spline_geoflux_data(trim(geqdsk_file), 64, 128)

    surface_values = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta_values = [0.0_dp, pi/3.0_dp, 2.0_dp*pi/3.0_dp, pi, &
        4.0_dp*pi/3.0_dp, 5.0_dp*pi/3.0_dp]

    open(newunit=unit, file=trim(output_file), status='replace', &
        action='write', form='formatted')

    do i_surf = 1, n_surface
        do i_sample = 0, n_surface_samples
            theta = two_pi * real(i_sample, dp) / real(n_surface_samples, dp)
            x_geo = [surface_values(i_surf), theta, 0.0_dp]
            call geoflux_to_cyl(x_geo, x_cyl)
            write(unit, '(i1,1x,i3,1x,f12.6,1x,f12.6)') 0, i_surf, x_cyl(1), x_cyl(3)
        end do
    end do

    do i_theta = 1, n_theta_lines
        do i_sample = 0, n_radial_samples
            s_value = s_min + (s_max - s_min) * real(i_sample, dp) &
                / real(n_radial_samples, dp)
            x_geo = [s_value, theta_values(i_theta), 0.0_dp]
            call geoflux_to_cyl(x_geo, x_cyl)
            write(unit, '(i1,1x,i3,1x,f12.6,1x,f12.6)') 1, i_theta, x_cyl(1), x_cyl(3)
        end do
    end do

    close(unit)

contains

    subroutine init_argument(index, value)
        integer, intent(in) :: index
        character(len=*), intent(out) :: value

        integer :: status, length

        value = ''
        call get_command_argument(index, value=value, length=length, status=status)
        if (status /= 0 .or. length <= 0) then
            write(*, '("Usage: geoflux_coord_dump <geqdsk> <output>")')
            error stop
        end if
        value = adjustl(value)
    end subroutine init_argument

end program geoflux_coord_dump
