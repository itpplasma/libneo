program test_stretch_coords
    use util_for_test, only: print_test, print_ok, print_fail
    use libneo_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

    integer :: argc, arg_len
    character(len=:), allocatable :: executable_path
    character(len=:), allocatable :: arg
    character(len=*), parameter :: empty_file = "test_convexwall_empty.dat"

    call get_command_argument(0, length=arg_len)
    allocate(character(len=max(1, arg_len)) :: executable_path)
    if (arg_len > 0) then
        call get_command_argument(0, executable_path)
    else
        executable_path = "./test_stretch_coords.x"
    end if

    argc = command_argument_count()
    if (argc > 0) then
        call get_command_argument(1, length=arg_len)
        allocate(character(len=max(1, arg_len)) :: arg)
        call get_command_argument(1, arg)
        if (trim(arg) == "--empty") then
            call test_stretch_coords_empty_file(empty_file)
            stop
        end if
        if (allocated(arg)) deallocate(arg)
    end if

    call test_stretch_coords_large_file()
    call verify_empty_file_failure(executable_path, empty_file)

contains

    subroutine test_stretch_coords_large_file()
        use field_sub, only: stretch_coords
        use input_files, only: convexfile
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(dp), parameter :: pi_value = 2.0_dp * asin(1.0_dp)
        real(dp), parameter :: major_radius = 10.0_dp
        real(dp), parameter :: radius_small = 0.8_dp
        real(dp), parameter :: radius_large = 0.8_dp
        real(dp), parameter :: angle_probe = 5.5_dp
        real(dp), parameter :: rho_probe = 1.0_dp
        real(dp), parameter :: vertical_shift = 75.0_dp
        real(dp), parameter :: R0 = major_radius + 0.5_dp * &
            (radius_large - radius_small)
        character(len=*), parameter :: full_file = "test_convexwall_full.dat"

        integer :: unit_full
        integer :: i
        real(dp) :: angle
        real(dp) :: rm_inside, zm_inside
        real(dp) :: rm_full, zm_full
        real(dp) :: r_input, z_input
        logical :: was_stretched

        call print_test("test_stretch_coords_large_file")

        ! A 150-point convex wall exercises the dynamic (arbitrary-size)
        ! allocation in stretch_coords (issue #109).  The wall is read once per
        ! run -- convexfile is never changed mid-run -- so this loads a single
        ! large wall and checks the stretched point is finite and physical.
        open(newunit=unit_full, file=full_file, status='replace', action='write')
        do i = 1, 150
            angle = (real(i, dp) - 0.5_dp) * 2.0_dp * pi_value / 150.0_dp
            call write_point(unit_full, major_radius, radius_small, &
                radius_large, vertical_shift, angle)
        end do
        close(unit_full)

        r_input = R0 + rho_probe * cos(angle_probe)
        z_input = vertical_shift + rho_probe * sin(angle_probe)

        convexfile = full_file
        call stretch_coords(major_radius, vertical_shift, rm_inside, zm_inside, &
            was_stretched)
        if (was_stretched) then
            write(error_unit, '(A)') &
                'ERROR: stretch_coords reported an inside point as outside.'
            call print_fail()
            error stop
        end if
        call stretch_coords(r_input, z_input, rm_full, zm_full, was_stretched)

        call cleanup_file(full_file)

        if (.not. ieee_is_finite(rm_full) .or. .not. ieee_is_finite(zm_full) &
            .or. rm_full <= 0.0_dp .or. .not. was_stretched) then
            write(error_unit, '(A)') &
                'ERROR: stretch_coords did not report and map the outside point.'
            call print_fail()
            error stop
        end if
        if (abs((rm_full-R0)*sin(angle_probe) - &
            (zm_full-vertical_shift)*cos(angle_probe)) > 1.0e-12_dp) then
            write(error_unit, '(A)') &
                'ERROR: stretch_coords changed the ray under vertical translation.'
            call print_fail()
            error stop
        end if

        call print_ok()
    end subroutine test_stretch_coords_large_file

    subroutine verify_empty_file_failure(executable_path, empty_filename)
        character(len=*), intent(in) :: executable_path
        character(len=*), intent(in) :: empty_filename

        character(len=512) :: command
        integer :: cmd_status
        logical :: empty_exists

        call cleanup_file(empty_filename)
        call print_test("test_stretch_coords_empty_file")
        command = ' '
        command = trim(executable_path) // " --empty"
        call execute_command_line(trim(command), cmdstat=cmd_status, wait=.true.)
        if (cmd_status /= 0) then
            write(error_unit, '(A, I0)') &
                'ERROR: execute_command_line failed with cmdstat ', cmd_status
            call print_fail()
            call cleanup_file(empty_filename)
            error stop
        end if
        inquire(file=empty_filename, exist=empty_exists)
        if (.not. empty_exists) then
            write(error_unit, '(A)') &
                'ERROR: stretch_coords succeeded with empty convex wall file.'
            call print_fail()
            call cleanup_file(empty_filename)
            error stop
        end if
        call cleanup_file(empty_filename)
        call print_ok()
    end subroutine verify_empty_file_failure

    subroutine write_point(unit, major_radius, radius_small, radius_large, &
        vertical_shift, angle)
        integer, intent(in) :: unit
        real(dp), intent(in) :: major_radius, radius_small, radius_large
        real(dp), intent(in) :: vertical_shift, angle
        real(dp) :: cos_angle, sin_angle, radius_profile

        cos_angle = cos(angle)
        sin_angle = sin(angle)
        if (cos_angle > 0.0_dp) then
            radius_profile = radius_large
        else
            radius_profile = radius_small
        end if
        write(unit, '(2es24.16)') major_radius + radius_profile * cos_angle, &
            vertical_shift + radius_profile * sin_angle
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

    subroutine test_stretch_coords_empty_file(empty_filename)
        use field_sub, only: stretch_coords
        use input_files, only: convexfile

        character(len=*), intent(in) :: empty_filename
        integer :: unit
        real(dp) :: rm_dummy, zm_dummy

        open(newunit=unit, file=empty_filename, status='replace', action='write')
        close(unit)

        convexfile = empty_filename
        call stretch_coords(0.0_dp, 0.0_dp, rm_dummy, zm_dummy)

        call cleanup_file(empty_filename)
        write(error_unit, '(A)') &
            'ERROR: stretch_coords returned successfully for empty convex wall file.'
        call print_fail()
        error stop
    end subroutine test_stretch_coords_empty_file

end program test_stretch_coords
