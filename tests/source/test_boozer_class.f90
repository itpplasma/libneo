program test_boozer_class

  use io, only: boozer_data_type

  ! Read file, write file

  type(boozer_data_type) :: boozer_file_orig, boozer_file_read_again

  real, parameter :: tolerance = 1.0e-12
  logical :: error_found = .false.

  call boozer_file_orig%read_data('../tests/resources/input_boozer_file.dat')

  call boozer_file_orig%write_data('output_boozer_file.dat')

  call boozer_file_read_again%read_data('output_boozer_file.dat')

  if (boozer_file_orig%get_m0b() .le. 0) then
    write(*,*) 'STOP size m0b equal or smaller than zero.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_n0b() .lt. 0) then
    write(*,*) 'STOP size n0b smaller than zero.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_nsurf() .le. 0) then
    write(*,*) 'STOP size nsurf equal or smaller than zero.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_nper() .le. 0) then
    write(*,*) 'STOP size nper equal or smaller than zero.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_m0b() /= boozer_file_read_again%get_m0b()) then
    write(*,*) 'STOP size m0b does not match.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_n0b() /= boozer_file_read_again%get_n0b()) then
    write(*,*) 'STOP size n0b does not match.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_nsurf() /= boozer_file_read_again%get_nsurf()) then
    write(*,*) 'STOP size nsurf does not match.'
    error_found = .true.
  end if

  if (boozer_file_orig%get_nper() /= boozer_file_read_again%get_nper()) then
    write(*,*) 'STOP size nper does not match.'
    error_found = .true.
  end if

  if (abs(boozer_file_orig%get_flux() - boozer_file_read_again%get_flux()) > tolerance) then
    write(*,*) 'STOP flux does not match.'
    error_found = .true.
  end if

  if (abs(boozer_file_orig%get_a() - boozer_file_read_again%get_a()) > tolerance) then
    write(*,*) 'STOP a does not match.'
    error_found = .true.
  end if

  if (abs(boozer_file_orig%get_R() - boozer_file_read_again%get_R()) > tolerance) then
    write(*,*) 'STOP R does not match.'
    error_found = .true.
  end if

  if (error_found) then
    stop 'STOP error occured'
  end if

end program test_boozer_class
