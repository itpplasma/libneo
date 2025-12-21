program test_efit_class

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use io, only: efit_data_type

    implicit none

    type(efit_data_type) :: efit_file_orig, efit_file_read_again

    real(dp), parameter :: tolerance = 1.0d-12
    logical :: error_found
    integer :: clock_count
    character(len=256) :: out_path

    error_found = .false.
    call system_clock(count=clock_count)
    write (out_path, '(a,i0,a)') '/tmp/libneo_test_efit_', clock_count, '.dat'

    call efit_file_orig%read_data('../../test/resources/input_efit_file.dat')

    call execute_command_line('rm -f ' // trim(out_path), wait=.true.)
    call efit_file_orig%write_data(trim(out_path))

    call efit_file_read_again%read_data(trim(out_path))

    if (efit_file_orig%get_nwEQD() .le. 0) then
        write (*, *) 'STOP size nw equal or smaller than zero.'
        error_found = .true.
    end if

    if (efit_file_orig%get_nhEQD() .le. 0) then
        write (*, *) 'STOP size nh equal or smaller than zero.'
        error_found = .true.
    end if

    if (efit_file_orig%get_nwEQD() /= efit_file_read_again%get_nwEQD()) then
        write (*, *) 'STOP size nw does not match.'
        error_found = .true.
    end if

    if (efit_file_orig%get_nhEQD() /= efit_file_read_again%get_nhEQD()) then
        write (*, *) 'STOP size nh does not match.'
        error_found = .true.
    end if

    if (abs(efit_file_orig%get_psiSep() - efit_file_read_again%get_psiSep()) > &
        tolerance) then
        write (*, *) 'STOP psiSep does not match.'
        error_found = .true.
    end if

    if (abs(efit_file_orig%get_bt0() - efit_file_read_again%get_bt0()) > tolerance) then
        write (*, *) 'STOP bt0 does not match.'
        error_found = .true.
    end if

    if (abs(efit_file_orig%get_rzero() - efit_file_read_again%get_rzero()) > &
        tolerance) then
        write (*, *) 'STOP rzero does not match.'
        error_found = .true.
    end if

    call execute_command_line('rm -f ' // trim(out_path), wait=.true.)

    if (error_found) then
        stop 'STOP error occurred'
    end if

end program test_efit_class
