program vacfield

  use iso_fortran_env, only: dp => real64, error_unit
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use math_constants, only: current_si_to_cgs, C
  use coil_tools, only: coil_t, coil_init, coil_deinit, coils_append, &
    process_fixed_number_of_args, check_number_of_args, &
    coils_read_AUG, coils_read_Nemov, coils_read_GPEC, &
    read_currents, Biot_Savart_sum_coils, write_Bvac_Nemov, &
    Biot_Savart_Fourier, write_Bnvac_Fourier, &
    Vector_Potential_Biot_Savart_Fourier, write_An_arrays

  implicit none

  real(dp) :: Rmin, Rmax, Zmin, Zmax, prefactor
  integer :: nR, nZ, nphi, nmax
  namelist /coil_field/ Rmin, Rmax, Zmin, Zmax, nR, nZ, nphi, nmax, prefactor
  integer :: argc, fid, status, num_coilfiles, kc, ncoil
  character(len = 1024) :: arg, err_msg, &
    coil_type, field_type, grid_file, field_file, currents_file
  character(len = :), dimension(:), allocatable :: coil_files
  type(coil_t), dimension(:), allocatable :: coils, more_coils
  real(dp), allocatable :: Ic(:), Bvac(:, :, :, :)
  complex(dp), allocatable :: Bnvac(:, :, :, :, :)
  complex(dp), dimension(:, :, :, :), allocatable :: AnR_array, Anphi_array, AnZ_array, dAnphi_dR_array, dAnphi_dZ_array

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
  call process_fixed_number_of_args(2, num_coilfiles, coil_files)
  call check_number_of_args(3 + num_coilfiles)
  call get_command_argument(3 + num_coilfiles, field_type)
  call check_number_of_args(4 + num_coilfiles)
  call get_command_argument(4 + num_coilfiles, grid_file)

  ! set default values for ASDEX Upgrade
  ! Rmin = 75.0d0
  ! Rmax = 267.0d0
  ! Zmin = -154.0d0
  ! Zmax = 150.4d0
  ! nR = 150
  ! nZ = 300
  ! nphi = 512
  ! nmax = 64
  ! prefactor = 5 * current_si_to_cgs / C

  ! read namelist input
  open(newunit = fid, file = trim(grid_file), status = 'old', action = 'read')
  read (fid, nml = coil_field, iostat = status, iomsg = err_msg)
  close(fid)
  if (status /= 0) then
    write (error_unit, '("Error ", i0, " occurred while reading from ", a, ": ", a, ".")') &
      status, trim(grid_file), trim(err_msg)
    error stop
  end if

  ! read coil files
  if (coil_type == 'AUG') then
    allocate(coils(num_coilfiles))
    do kc = 1, num_coilfiles
      call coils_read_AUG(trim(coil_files(kc)), coils(kc))
    end do
  else if (coil_type == 'GPEC') then
    call coils_read_GPEC(trim(coil_files(1)), coils)
    do kc = 2, num_coilfiles
      call coils_read_GPEC(trim(coil_files(kc)), more_coils)
      call coils_append(coils, more_coils)
    end do
  else if (coil_type == 'Nemov') then
    call coils_read_Nemov(trim(coil_files(1)), coils)
    do kc = 2, num_coilfiles
      call coils_read_Nemov(trim(coil_files(kc)), more_coils)
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
    call check_number_of_args(5 + num_coilfiles)
    call get_command_argument(5 + num_coilfiles, field_file)
    call h5_init
    h5overwrite = .true.
    call Biot_Savart_Fourier(coils, nmax, &
      Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bnvac)
    call write_Bnvac_Fourier(trim(field_file), Bnvac, Rmin, Rmax, Zmin, Zmax)
    deallocate(Bnvac)
    call h5_deinit
  else if (field_type == 'sum') then
    call check_number_of_args(5 + num_coilfiles)
    call get_command_argument(5 + num_coilfiles, field_file)
    call check_number_of_args(6 + num_coilfiles)
    call get_command_argument(6 + num_coilfiles, currents_file)
    call read_currents(trim(currents_file), Ic)
    call Biot_Savart_sum_coils(coils, prefactor * Ic, &
      Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
    call write_Bvac_Nemov(trim(field_file), Rmin, Rmax, Zmin, Zmax, Bvac)
    deallocate(Ic)
    deallocate(Bvac)
  else if (field_type == 'vector_potential') then
    print*, 'hello'
    call check_number_of_args(5 + num_coilfiles)
    call get_command_argument(5 + num_coilfiles, field_file)
    call Vector_Potential_Biot_Savart_Fourier(coils, nmax, &
      Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR_array, Anphi_array, AnZ_array, dAnphi_dR_array, dAnphi_dZ_array, ncoil)
    call write_An_arrays(trim(field_file), ncoil, nmax, &
      Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, AnR_array, Anphi_array, AnZ_array, dAnphi_dR_array, dAnphi_dZ_array)
    !deallocate(Bnvac)
    print*, 'hello'
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
    write (*, '("             <field_format> <grid_file> <field_file> [ <current_file> ]")')
    write (*, '()')
    write (*, '("  coil_format   supportedd formats are AUG, GPEC, and Nemov")')
    write (*, '("  n_files       integer indicating the number of coil file names that follow")')
    write (*, '("  file_#        coil file, to be processed in the given order")')
    write (*, '("  field_format  ''sum'' sums over coils and currents, ")')
    write (*, '("                ''Fourier'' computes Fourier modes per coil with unit current")')
    write (*, '("  grid_file     input file for namelist specifying grid parameters")')
    write (*, '("  field_file    output file, plaintext for ''sum'' and HDF5 for ''Fourier''")')
    write (*, '("  current_file  plaintext file containing currents for use with ''sum''")')
  end subroutine usage

end program vacfield
