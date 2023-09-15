program vacfield

  use iso_fortran_env, only: dp => real64, error_unit
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use coil_tools, only: AUG_coils_read, AUG_coils_write_Nemov, AUG_coils_read_Nemov, &
       AUG_coils_write_GPEC, AUG_coils_read_GPEC, AUG_coils_write_Fourier, read_currents_Nemov, &
       Biot_Savart_sum_coils, write_Bvac_Nemov

  implicit none

  ! use default values for now
  real(dp), parameter :: Rmin = 75.0d0, Rmax = 267.0d0, Zmin = -154.0d0, Zmax = 150.4d0
  integer, parameter :: nmax = 64, nR = 150, nZ = 300, nphi = 512
  integer :: argc, nseg, ncoil, nwind
  character(len = 1024) :: in_type, out_type, in_dir, out_dir
  real(dp), allocatable :: XYZ(:, :, :), Ic(:), Bvac(:, :, :, :)

  argc = command_argument_count()
  if (argc >= 4) then
     call get_command_argument(1, in_type)
     call get_command_argument(2, out_type)
     call get_command_argument(3, in_dir)
     call get_command_argument(4, out_dir)
  else
     error stop 'Expected 4 command line arguments, see documentation'
  endif

  if (in_type == 'AUG') then
     call AUG_coils_read(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else if (in_type == 'GPEC') then
     call AUG_coils_read_GPEC(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else if (in_type == 'Nemov') then
     call AUG_coils_read_Nemov(trim(in_dir), ncoil, nseg, nwind, XYZ)
  else
     write (error_unit, '("unknown input type ", a)') trim(in_type)
     error stop
  end if
  if (out_type == 'Fourier') then
     call h5_init
     h5overwrite = .true.
     call AUG_coils_write_Fourier(trim(out_dir), ncoil, nseg, nwind, XYZ, &
          nmax, Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ)
     call h5_deinit
  else if (out_type == 'GPEC') then
     call AUG_coils_write_GPEC(trim(out_dir), ncoil, nseg, nwind, XYZ)
  else if (out_type == 'Nemov') then
     call AUG_coils_write_Nemov(trim(out_dir), ncoil, nseg, XYZ)
  else if (out_type == 'sum') then
     call read_currents_Nemov(trim(in_dir), Ic)
     call Biot_Savart_sum_coils(ncoil, nseg, nwind, XYZ, Ic, &
          Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
     call write_Bvac_Nemov(trim(out_dir), Rmin, Rmax, Zmin, Zmax, Bvac)
     if (allocated(Ic)) deallocate(Ic)
     if (allocated(Bvac)) deallocate(Bvac)
  else
     write (error_unit, '("unknown output type ", a)') trim(out_type)
     error stop
  end if
  if (allocated(XYZ)) deallocate(XYZ)

end program vacfield
