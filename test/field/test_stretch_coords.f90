program test_stretch_coords
    use util_for_test, only: print_test, print_ok, print_fail
    use libneo_kinds, only: dp
    implicit none

    call test_stretch_coords_large_file
    call test_stretch_coords_empty_file

contains

    subroutine test_stretch_coords_large_file
        use field_sub, only: stretch_coords
        use input_files, only: convexfile

        real(dp), parameter :: tolerance = 1.0e-9_dp
        real(dp), parameter :: pi_value = 2.0_dp * asin(1.0_dp)
        real(dp), parameter :: major_radius = 10.0_dp
        real(dp), parameter :: radius_small = 0.4_dp
        real(dp), parameter :: radius_large = 0.8_dp
        real(dp), parameter :: angle_probe = 5.5_dp
        real(dp), parameter :: rho_probe = 1.0_dp
        real(dp), parameter :: R0 = major_radius + 0.5_dp * &
            (radius_large - radius_small)
        character(len=*), parameter :: full_file = "test_convexwall_full.dat"
        character(len=*), parameter :: truncated_file = "test_convexwall_truncated.dat"

        integer :: unit_full, unit_trunc
        integer :: i
        real(dp) :: angle
        real(dp) :: rm_full, zm_full, rm_trunc, zm_trunc
        real(dp) :: r_input, z_input

        call print_test("test_stretch_coords_large_file")

        open(newunit=unit_full, file=full_file, status='replace', action='write')
        open(newunit=unit_trunc, file=truncated_file, status='replace', action='write')
        do i = 1, 150
            angle = (real(i, dp) - 0.5_dp) * 2.0_dp * pi_value / 150.0_dp
            call write_point(unit_full, major_radius, radius_small, radius_large, angle)
            if (i <= 100) then
                call write_point(unit_trunc, major_radius, radius_small, &
                                 radius_large, angle)
            end if
        end do
        close(unit_full)
        close(unit_trunc)

        r_input = R0 + rho_probe * cos(angle_probe)
        z_input = rho_probe * sin(angle_probe)

        convexfile = truncated_file
        call stretch_coords(r_input, z_input, rm_trunc, zm_trunc)

        convexfile = full_file
        call stretch_coords(r_input, z_input, rm_full, zm_full)

        call cleanup_file(full_file)
        call cleanup_file(truncated_file)

        if (abs(rm_full - rm_trunc) <= tolerance .and. &
            abs(zm_full - zm_trunc) <= tolerance) then
            write(*, '(A)') &
                'ERROR: stretch_coords results identical with truncated data.'
            call print_fail
            error stop
        end if

        call print_ok
    end subroutine test_stretch_coords_large_file

    subroutine write_point(unit, major_radius, radius_small, radius_large, angle)
        integer, intent(in) :: unit
        real(dp), intent(in) :: major_radius, radius_small, radius_large, angle
        real(dp) :: cos_angle, sin_angle, radius_profile

        cos_angle = cos(angle)
        sin_angle = sin(angle)
        if (cos_angle > 0.0_dp) then
            radius_profile = radius_large
        else
            radius_profile = radius_small
        end if
        write(unit, '(2es24.16)') major_radius + radius_profile * cos_angle, &
                                   radius_profile * sin_angle
    end subroutine write_point

    subroutine cleanup_file(filename)
        character(len=*), intent(in) :: filename
        logical :: exists
        integer :: unit

        inquire(file=filename, exist=exists)
        if (exists) then
            open(newunit=unit, file=filename, status='old')
            close(unit, status='delete')
        end if
    end subroutine cleanup_file

    subroutine test_stretch_coords_empty_file
        use input_files, only: convexfile

        character(len=*), parameter :: empty_file = "empty_convexwall.dat"
        integer :: unit

        call print_test("test_stretch_coords_empty_file")

        open(newunit=unit, file=empty_file, status='replace', action='write')
        close(unit)

        convexfile = empty_file

        call cleanup_file(empty_file)

        call print_ok
    end subroutine test_stretch_coords_empty_file

end program test_stretch_coords
