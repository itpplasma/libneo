program gframe_from_geoflux_dump
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use math_constants, only: pi
    use geoflux_field, only: spline_geoflux_data
    use geoflux_coordinates, only: geoflux_to_cyl
    use gframe_boundary, only: gframe_boundary_t, allocate_boundary
    use libneo_coordinates, only: coordinate_system_t, gframe_coordinate_system_t, &
        make_gframe_coordinate_system
    implicit none

    integer, parameter :: max_len = 512
    integer, parameter :: n_surface = 5
    integer, parameter :: n_surface_samples = 180
    integer, parameter :: n_theta_lines = 6
    integer, parameter :: n_radial_samples = 120
    integer, parameter :: ntheta_boundary = 65
    integer, parameter :: nzeta_boundary = 1
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
    type(gframe_boundary_t) :: boundary
    class(coordinate_system_t), allocatable :: cs
    type(gframe_coordinate_system_t), pointer :: gcs
    real(dp) :: u_g(3), x_g(3)

    call init_argument(1, geqdsk_file)
    call init_argument(2, output_file)

    call spline_geoflux_data(trim(geqdsk_file), 64, 128)
    call build_boundary_from_geoflux(boundary)
    call make_gframe_coordinate_system(cs, boundary, nrho = 65)

    select type (gcs => cs)
    type is (gframe_coordinate_system_t)
        continue
    class default
        write(*, '(a)') "gframe_from_geoflux_dump: allocation failed"
        error stop
    end select

    surface_values = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta_values = [0.0_dp, pi/3.0_dp, 2.0_dp*pi/3.0_dp, pi, &
        4.0_dp*pi/3.0_dp, 5.0_dp*pi/3.0_dp]

    open(newunit=unit, file=trim(output_file), status='replace', &
        action='write', form='formatted')

    do i_surf = 1, n_surface
        do i_sample = 0, n_surface_samples
            theta = two_pi * real(i_sample, dp) / real(n_surface_samples, dp)
            u_g = [surface_values(i_surf), theta, 0.0_dp]
            call gcs%evaluate_point(u_g, x_g)
            write(unit, '(i1,1x,i3,1x,f12.6,1x,f12.6)') 0, i_surf, x_g(1), x_g(3)
        end do
    end do

    do i_theta = 1, n_theta_lines
        do i_sample = 0, n_radial_samples
            s_value = s_min + (s_max - s_min) * real(i_sample, dp) &
                / real(n_radial_samples, dp)
            u_g = [s_value, theta_values(i_theta), 0.0_dp]
            call gcs%evaluate_point(u_g, x_g)
            write(unit, '(i1,1x,i3,1x,f12.6,1x,f12.6)') 1, i_theta, x_g(1), x_g(3)
        end do
    end do

    close(unit)

contains

    subroutine build_boundary_from_geoflux(boundary)
        type(gframe_boundary_t), intent(out) :: boundary

        integer :: jt
        real(dp) :: theta_val
        real(dp) :: x_geo(3), x_cyl(3)

        call allocate_boundary(boundary, ntheta_boundary, nzeta_boundary)
        boundary%zeta = 0.0_dp

        do jt = 1, ntheta_boundary
            theta_val = two_pi * real(jt - 1, dp) / real(ntheta_boundary - 1, dp)
            boundary%theta(jt) = theta_val
            x_geo = [1.0_dp, theta_val, 0.0_dp]
            call geoflux_to_cyl(x_geo, x_cyl)
            boundary%Rb(jt, 1) = x_cyl(1)
            boundary%Zb(jt, 1) = x_cyl(3)
        end do
    end subroutine build_boundary_from_geoflux

    subroutine init_argument(index, value)
        integer, intent(in) :: index
        character(len=*), intent(out) :: value

        integer :: status, length

        value = ''
        call get_command_argument(index, value=value, length=length, status=status)
        if (status /= 0 .or. length <= 0) then
            write(*, '("Usage: gframe_from_geoflux_dump <geqdsk> <output>")')
            error stop
        end if
        value = adjustl(value)
    end subroutine init_argument

end program gframe_from_geoflux_dump
