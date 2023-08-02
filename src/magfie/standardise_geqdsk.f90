program standardise_geqdsk
  use geqdsk_tools, only : geqdsk_t, geqdsk_read, geqdsk_classify, &
    geqdsk_standardise, geqdsk_write, geqdsk_deinit

  implicit none

  character(len = 1024) :: unprocessed_file, standardised_file
  type(geqdsk_t) :: geqdsk

  if (command_argument_count() == 2) then
    call get_command_argument(1, unprocessed_file)
    call get_command_argument(2, standardised_file)
  else
    write(*,*) 'Error: unexpected number of arguments. Usage:'
    write(*,*) ''
    write(*,*) '  standardise_geqdsk.x infilename outfilename'
    write(*,*) ''
    write(*,*) 'infilename: string with path+name to unprocessed GEQDSK file which to standardise'
    write(*,*) 'outfilename: string with path+name of GEQDSK file where to write the output, i.e. the standardised file.'
    error stop 'unexpected number of parameters'
  end if

  call geqdsk_read(geqdsk, trim(unprocessed_file))
  call geqdsk_classify(geqdsk)
  call geqdsk_standardise(geqdsk)
  call geqdsk_write(geqdsk, trim(standardised_file))
  call geqdsk_deinit(geqdsk)
end program standardise_geqdsk
