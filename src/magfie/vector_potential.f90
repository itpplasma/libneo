program vector_potential

  use iso_fortran_env, only: dp => real64, error_unit
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use math_constants, only: current_si_to_cgs, C
  use coil_tools, only: coil_t, coil_init, coil_deinit, coils_append, &
    process_fixed_number_of_args, check_number_of_args, &
    coils_read_AUG, coils_read_Nemov, coils_read_GPEC, &
    read_currents, Biot_Savart_sum_coils, write_Bvac_Nemov, &
    Biot_Savart_Fourier, write_Bnvac_Fourier

  implicit none

  real(dp) :: Rmin, Rmax, Zmin, Zmax, prefactor
  integer :: nR, nZ, nphi, nmax
  integer :: argc, fid, status, ncoil, kc
  character(len = 1024) :: arg, err_msg, &
    coil_type, field_type, grid_file, field_file, currents_file
  character(len = :), dimension(:), allocatable :: coil_files
  type(coil_t), dimension(:), allocatable :: coils, more_coils
  real(dp), allocatable :: Ic(:), Bvac(:, :, :, :)
  complex(dp), allocatable :: Bnvac(:, :, :, :, :)

  !procedure 
  !1) read biotsavart.inp
  !2) read coildata.m
  !3) compute vector potential
  !4) write netcdf file

  ! read biotsavart.inp
  open(newunit = fid, file = 'biotsavart.inp', status = 'old', action = 'read')
  read (fid, nml = coil_field, iostat = status, iomsg = err_msg)
  close(fid)
  if (status /= 0) then
    write (error_unit, '("Error ", i0, " occurred while reading from ", a, ": ", a, ".")') &
      status, trim(grid_file), trim(err_msg)
    error stop
  end if

  do ks = 1, nseg
    read (fid, *) R_Z_phi(ks, :)
  end do

  ! read coil files
    allocate(coils(ncoil))
    do kc = 1, ncoil
      call coils_read_AUG(trim(coil_files(kc)), coils(kc))
    end do

  ! compute vector potential

  ! write netcdf file

end program vector_potential
