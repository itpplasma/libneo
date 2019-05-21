program test_efit_class

  use io, only: efit_data_type

  ! Read file, write file

  type(efit_data_type) :: efit_file_orig, efit_file_read_again

  real, parameter :: tolerance = 1.0e-12
  logical :: error_found = .false.

  call efit_file_orig%read_data('../tests/resources/input_efit_file.dat')

  call efit_file_orig%write_data('output_efit_file.dat')

  call efit_file_read_again%read_data('output_efit_file.dat')

  if (efit_file_orig%get_nwEQD() .le. 0) then
    write(*,*) 'STOP size nw equal or smaller than zero.'
    error_found = .true.
  end if

  if (efit_file_orig%get_nhEQD() .le. 0) then
    write(*,*) 'STOP size nh equal or smaller than zero.'
    error_found = .true.
  end if

  if (efit_file_orig%get_nwEQD() /= efit_file_read_again%get_nwEQD()) then
    write(*,*) 'STOP size nw does not match.'
    error_found = .true.
  end if

  if (efit_file_orig%get_nhEQD() /= efit_file_read_again%get_nhEQD()) then
    write(*,*) 'STOP size nh does not match.'
    error_found = .true.
  end if

  if (abs(efit_file_orig%get_psiSep() - efit_file_read_again%get_psiSep()) > tolerance) then
    write(*,*) 'STOP psiSep does not match.'
    error_found = .true.
  end if

  if (abs(efit_file_orig%get_bt0() - efit_file_read_again%get_bt0()) > tolerance) then
    write(*,*) 'STOP bt0 does not match.'
    error_found = .true.
  end if

  if (abs(efit_file_orig%get_rzero() - efit_file_read_again%get_rzero()) > tolerance) then
    write(*,*) 'STOP rzero does not match.'
    error_found = .true.
  end if

!~   if (abs(efit_file_orig%get_rad() - efit_file_read_again%get_rad()) > tolerance) then
!~     write(*,*) 'STOP rad does not match.'
!~     error_found = .true.
!~   end if

!~   if (abs(efit_file_orig%get_zet() - efit_file_read_again%get_zet()) > tolerance) then
!~     write(*,*) 'STOP zet does not match.'
!~     error_found = .true.
!~   end if

!~   if (abs(efit_file_orig%get_psiRZ() - efit_file_read_again%get_psiRZ()) > tolerance) then
!~     write(*,*) 'STOP psiRZ does not match.'
!~     error_found = .true.
!~   end if

  if (error_found) then
    stop 'STOP error occured'
  end if

end program test_efit_class
