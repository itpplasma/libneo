program tilted_coil_axis_field
  use iso_fortran_env, only: dp => real64
  use neo_biotsavart_field, only: biotsavart_field_t

  implicit none

  character(len=1024) :: coil_file, axis_file, output_file
  integer :: argc, npts, iu_axis, iu_out, i, ios
  real(dp), allocatable :: points(:, :)
  real(dp) :: x(3), B(3)
  type(biotsavart_field_t) :: field

  argc = command_argument_count()
  if (argc /= 3) then
    write (*, '(a)') 'Usage: tilted_coil_axis_field.x <coil_file> <axis_file> <output_file>'
    stop 1
  end if

  call get_command_argument(1, coil_file)
  call get_command_argument(2, axis_file)
  call get_command_argument(3, output_file)

  call field%biotsavart_field_init(trim(coil_file))

  open(newunit = iu_axis, file = trim(axis_file), status = 'old', action = 'read', iostat = ios)
  if (ios /= 0) then
    write (*, '(a,1x,a)') 'Failed to open axis file:', trim(axis_file)
    stop 2
  end if

  read(iu_axis, *, iostat = ios) npts
  if (ios /= 0 .or. npts < 1) then
    write (*, '(a)') 'Invalid axis file header'
    stop 3
  end if

  allocate(points(3, npts))
  do i = 1, npts
    read(iu_axis, *, iostat = ios) points(:, i)
    if (ios /= 0) then
      write (*, '(a,i0)') 'Failed to read axis point', i
      stop 4
    end if
  end do
  close(iu_axis)

  open(newunit = iu_out, file = trim(output_file), status = 'replace', action = 'write', iostat = ios)
  if (ios /= 0) then
    write (*, '(a,1x,a)') 'Failed to open output file:', trim(output_file)
    stop 5
  end if

  do i = 1, npts
    x = points(:, i)
    call field%compute_bfield(x, B)
    write(iu_out, '(3(es24.16e3,1x))') B
  end do
  close(iu_out)

end program tilted_coil_axis_field
