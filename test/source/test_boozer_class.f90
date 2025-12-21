program test_boozer_class

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use io, only: boozer_data_type

    implicit none

    type(boozer_data_type) :: boozer_file_orig, boozer_file_read_again

    real(dp), parameter :: tolerance = 1.0d-12
    logical :: error_found
    integer :: clock_count
    character(len=256) :: out_path

    error_found = .false.
    call system_clock(count=clock_count)
    write (out_path, '(a,i0,a)') '/tmp/libneo_test_boozer_', clock_count, '.dat'

    call boozer_file_orig%read_data('../../test/resources/input_boozer_file.dat')

    call execute_command_line('rm -f ' // trim(out_path), wait=.true.)
    call boozer_file_orig%write_data(trim(out_path))

    call boozer_file_read_again%read_data(trim(out_path))

    if (boozer_file_orig%get_m0b() .le. 0) then
        write (*, *) 'STOP size m0b equal or smaller than zero.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_n0b() .lt. 0) then
        write (*, *) 'STOP size n0b smaller than zero.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_nsurf() .le. 0) then
        write (*, *) 'STOP size nsurf equal or smaller than zero.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_nper() .le. 0) then
        write (*, *) 'STOP size nper equal or smaller than zero.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_m0b() /= boozer_file_read_again%get_m0b()) then
        write (*, *) 'STOP size m0b does not match.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_n0b() /= boozer_file_read_again%get_n0b()) then
        write (*, *) 'STOP size n0b does not match.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_nsurf() /= boozer_file_read_again%get_nsurf()) then
        write (*, *) 'STOP size nsurf does not match.'
        error_found = .true.
    end if

    if (boozer_file_orig%get_nper() /= boozer_file_read_again%get_nper()) then
        write (*, *) 'STOP size nper does not match.'
        error_found = .true.
    end if

    if (abs(boozer_file_orig%get_flux() - boozer_file_read_again%get_flux()) > &
        tolerance) then
        write (*, *) 'STOP flux does not match.'
        error_found = .true.
    end if

    if (abs(boozer_file_orig%get_a() - boozer_file_read_again%get_a()) > tolerance) then
        write (*, *) 'STOP a does not match.'
        error_found = .true.
    end if

    if (abs(boozer_file_orig%get_R() - boozer_file_read_again%get_R()) > tolerance) then
        write (*, *) 'STOP R does not match.'
        error_found = .true.
    end if

    call execute_command_line('rm -f ' // trim(out_path), wait=.true.)

    if (error_found) then
        stop 'STOP error occurred'
    end if

end program test_boozer_class
