program vacfield

  use iso_fortran_env, only: dp => real64, error_unit
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use coil_tools, only: coil_t, coil_init, coil_deinit, coils_append, &
       process_fixed_number_of_args, check_number_of_args, &
       AUG_coils_read, AUG_coils_read_Nemov, AUG_coils_read_GPEC, &
       AUG_coils_write_Fourier, read_currents_Nemov, &
       Biot_Savart_sum_coils, write_Bvac_Nemov

  implicit none

  ! use default values for now
  real(dp), parameter :: Rmin = 75.0d0, Rmax = 267.0d0, Zmin = -154.0d0, Zmax = 150.4d0
  integer, parameter :: nmax = 64, nR = 150, nZ = 300, nphi = 512
  integer :: argc, ncoil, kc
  character(len = 1024) :: arg, coil_type, field_type, field_file, currents_file
  character(len = :), dimension(:), allocatable :: coil_files
  type(coil_t), dimension(:), allocatable :: coils, more_coils
  real(dp), allocatable :: Ic(:), Bvac(:, :, :, :)

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
  call get_command_argument(1, coil_type)
  call process_fixed_number_of_args(2, ncoil, coil_files)
  call check_number_of_args(3 + ncoil)
  call get_command_argument(3 + ncoil, field_type)

  if (coil_type == 'AUG') then
     allocate(coils(ncoil))
     do kc = 1, ncoil
        call AUG_coils_read(trim(coil_files(kc)), coils(kc))
     end do
  else if (coil_type == 'GPEC') then
     call AUG_coils_read_GPEC(trim(coil_files(1)), coils)
     do kc = 2, ncoil
        call AUG_coils_read_GPEC(trim(coil_files(kc)), more_coils)
        call coils_append(coils, more_coils)
     end do
  else if (coil_type == 'Nemov') then
     call AUG_coils_read_Nemov(trim(coil_files(1)), coils)
     do kc = 2, ncoil
        call AUG_coils_read_Nemov(trim(coil_files(kc)), more_coils)
        call coils_append(coils, more_coils)
     end do
  else
     write (error_unit, '("Unknown input type ", a, ".")') trim(coil_type)
     error stop
  end if
  ! not actually necessary for format AUG
  if (.not. allocated(coils)) then
     write (error_unit, '("Unexpected error: coils array not allocated.")')
     error stop
  end if

  if (field_type == 'Fourier') then
     call check_number_of_args(4 + ncoil)
     call get_command_argument(4 + ncoil, field_file)
     call h5_init
     h5overwrite = .true.
     call AUG_coils_write_Fourier(trim(field_file), coils, &
          nmax, Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ)
     call h5_deinit
  else if (field_type == 'sum') then
     call check_number_of_args(4 + ncoil)
     call get_command_argument(4 + ncoil, field_file)
     call check_number_of_args(5 + ncoil)
     call get_command_argument(5 + ncoil, currents_file)
     call read_currents_Nemov(trim(currents_file), Ic)
     call Biot_Savart_sum_coils(coils, Ic, &
          Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
     call write_Bvac_Nemov(trim(field_file), Rmin, Rmax, Zmin, Zmax, Bvac)
     if (allocated(Ic)) deallocate(Ic)
     if (allocated(Bvac)) deallocate(Bvac)
  else
     write (error_unit, '("unknown output type ", a)') trim(field_type)
     error stop
  end if

  if (allocated(coils)) then
     do kc = 1, size(coils)
        call coil_deinit(coils(kc))
     end do
     deallocate(coils)
  end if

contains

  subroutine usage()
    write (*, '("Usage:")')
    write (*, '("  vacfield.x <coil_format> <n_files> <file_1> [ <file_2> ... ]")')
    write (*, '("             <field_format> <field_file> [ <current_file> ]")')
    write (*, '()')
    write (*, '("  coil_format   supportedd formats are AUG, GPEC, and Nemov")')
    write (*, '("  n_files       integer indicating the number of coil file names that follow")')
    write (*, '("  file_#        coil file, to be processed in the given order")')
    write (*, '("  field_format  ''sum'' sums over coils and currents, ")')
    write (*, '("                ''Fourier'' computes Fourier modes per coil with unit current")')
    write (*, '("  field_file    output file, plaintext for ''sum'' and HDF5 for ''Fourier''")')
    write (*, '("  current_file  plaintext file containing currents for use with ''sum''")')
  end subroutine usage

end program vacfield
