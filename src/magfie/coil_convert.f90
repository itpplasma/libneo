program coil_convert
  use iso_fortran_env, only: error_unit
  use coil_tools, only: coil_t, coil_init, coil_deinit, coils_append, &
    process_fixed_number_of_args, check_number_of_args, &
    coils_write_AUG, coils_read_AUG, coils_writenemov, coils_readnemov, &
    coils_write_GPEC, coils_read_GPEC

  implicit none

  character(len = *), parameter :: incommensurable_fmt = &
    '("Cannot distribute ", i0, " coil(s) to ", i0, " output file(s).")'

  integer :: argc, in_ncoil, out_ncoil, kc
  character(len = 1024) :: arg, in_type, out_type
  character(len = 1024), dimension(:), allocatable :: in_files, out_files
  type(coil_t), dimension(:), allocatable :: coils, more_coils

  argc = command_argument_count()
  if (argc < 1) then
    call usage()
    stop
  end if
  if (argc == 1) then
    call get_command_argument(1, arg)
    if (arg == '-h' .or. arg == '--help' .or. arg == '-?') then
      call usage()
      stop
    else
      write (error_unit, '("Unrecognized argument ", a, ".")') trim(arg)
      call usage()
      error stop
    end if
  endif

  call check_number_of_args(1)
  call get_command_argument(1, in_type)
  call process_fixed_number_of_args(2, in_ncoil, in_files)
  call check_number_of_args(3 + in_ncoil)
  call get_command_argument(3 + in_ncoil, out_type)
  call process_fixed_number_of_args(4 + in_ncoil, out_ncoil, out_files)

  if (in_type == 'AUG') then
    allocate(coils(in_ncoil))
    do kc = 1, in_ncoil
      call coils_read_AUG(trim(in_files(kc)), coils(kc))
    end do
  else if (in_type == 'GPEC') then
    call coils_read_GPEC(trim(in_files(1)), coils)
    do kc = 2, in_ncoil
      call coils_read_GPEC(trim(in_files(kc)), more_coils)
      call coils_append(coils, more_coils)
    end do
  else if (in_type == 'Nemov') then
    call coils_readnemov(trim(in_files(1)), coils)
    do kc = 2, in_ncoil
      call coils_readnemov(trim(in_files(kc)), more_coils)
      call coils_append(coils, more_coils)
    end do
  else
    write (error_unit, '("Unknown input type ", a, ".")') trim(in_type)
    error stop
  end if
  ! not actually necessary for format AUG
  if (.not. allocated(coils)) then
    write (error_unit, '("Unexpected error: coils array not allocated.")')
    error stop
  end if
  in_ncoil = size(coils)

  if (out_type == 'AUG') then
    if (in_ncoil /= out_ncoil) then
      write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
      error stop
    end if
    do kc = 1, out_ncoil
      call coils_write_AUG(trim(out_files(kc)), coils(kc))
    end do
  else if (out_type == 'GPEC') then
    if (out_ncoil > 1) then
      if (in_ncoil /= out_ncoil) then
        write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
        error stop
      end if
      do kc = 1, out_ncoil
        call coils_write_GPEC(trim(out_files(kc)), coils(kc:kc))
      end do
    else
      call coils_write_GPEC(trim(out_files(1)), coils)
    end if
  else if (out_type == 'Nemov') then
    if (out_ncoil > 1) then
      if (in_ncoil /= out_ncoil) then
        write (error_unit, incommensurable_fmt) in_ncoil, out_ncoil
        error stop
      end if
      do kc = 1, out_ncoil
        call coils_writenemov(trim(out_files(kc)), coils(kc:kc))
      end do
    else
      call coils_writenemov(trim(out_files(1)), coils)
    end if
  else
    write (error_unit, '("Unknown output type ", a, ".")') trim(out_type)
    error stop
  end if

  if (allocated(coils)) then
    do kc = 1, size(coils)
      call coil_deinit(coils(kc))
    end do
    deallocate(coils)
  end if
  if (allocated(in_files)) deallocate(in_files)
  if (allocated(out_files)) deallocate(out_files)

contains

  subroutine usage()
    write (*, '("Usage:")')
    write (*, '("  coil_convert.x <coil_format> <n_files> <file_1> [ <file_2> ... ]")')
    write (*, '("                 <coil_format> <n_files> <file_1> [ <file_2> ... ]")')
    write (*, '()')
    write (*, '("  coil_format   supportedd formats are AUG, GPEC, and Nemov")')
    write (*, '("  n_files       integer indicating the number of coil file names that follow")')
    write (*, '("  file_#        coil file, to be processed in the given order")')
    write (*, '()')
    write (*, '("The first set of input arguments specify the input files,")')
    write (*, '("the second set of input arguments specify the output files.")')
  end subroutine usage

end program coil_convert
