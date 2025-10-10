program tilted_coil_fourier_modes
  use util_for_test, only: print_test, print_ok, print_fail
  use coil_tools, only: coil_t, coil_init, coil_deinit, &
                       vector_potential_biot_savart_fourier, write_Anvac_fourier, &
                       biot_savart_fourier, write_Bnvac_fourier, &
                       biot_savart_sum_coils, write_Bvac_nemov
  use libneo_kinds, only: dp
  use math_constants, only: length_si_to_cgs
  use hdf5_tools, only: h5_init, h5_deinit
  implicit none

  type(coil_t), allocatable :: coils(:)
  complex(dp), allocatable :: AnR(:,:,:,:), Anphi(:,:,:,:), AnZ(:,:,:,:)
  complex(dp), allocatable :: dAnphi_dR(:,:,:,:), dAnphi_dZ(:,:,:,:)
  complex(dp), allocatable :: Bn(:,:,:,:,:)
  real(dp), allocatable :: Bvac(:,:,:,:)
  real(dp), allocatable :: Ic(:)
  real(dp) :: Rmin, Rmax, Zmin, Zmax
  integer :: nR, nZ, nphi, nmax, nseg, k, iostat
  real(dp) :: min_distance, max_eccentricity
  logical :: use_convex_wall
  real(dp) :: x, y, z
  real(dp) :: scale

  call print_test("tilted coil Fourier field data generation")

  scale = length_si_to_cgs

  ! Grid parameters (converted to cm)
  Rmin = scale * 0.5_dp
  Rmax = scale * 4.0_dp
  Zmin = scale * (-1.5_dp)
  Zmax = scale * 2.5_dp
  nR = 50
  nZ = 50
  nphi = 192
  nmax = 16  ! Include modes n=0,1,2,...,16

  ! Biot-Savart parameters
  min_distance = 1.0e-6_dp
  max_eccentricity = 0.999_dp
  use_convex_wall = .false.

  ! Load coil geometry
  allocate(coils(1))

  ! Read number of segments from file
  open(unit=10, file='tilted_coil.dat', status='old', action='read')
  read(10, *) nseg

  ! Initialize coil
  call coil_init(coils(1), nseg, 1)

  ! Read coil coordinates
  do k = 1, nseg
    read(10, *, iostat=iostat) x, y, z
    if (iostat /= 0) then
      call print_fail
      error stop "Error reading coil geometry"
    end if
    coils(1)%XYZ(1, k) = scale * x
    coils(1)%XYZ(2, k) = scale * y
    coils(1)%XYZ(3, k) = scale * z
  end do
  close(10)

  ! Compute vector potential Fourier modes
  call vector_potential_biot_savart_fourier(coils, nmax, min_distance, max_eccentricity, use_convex_wall, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, &
    AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)

  ! Compute magnetic field Fourier modes (per-coil)
  call biot_savart_fourier(coils, nmax, Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bn)

  call h5_init()

  ! Write to NetCDF and HDF5 files
  call write_Anvac_fourier('tilted_coil_Anvac.nc', 1, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, &
    AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)
  call write_Bnvac_fourier('tilted_coil_Bnvac.h5', Bn, Rmin, Rmax, Zmin, Zmax)

  allocate(Ic(1))
  Ic(1) = 1.0_dp
  call biot_savart_sum_coils(coils, Ic, Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, Bvac)
  call write_Bvac_nemov('tilted_coil_Bvac.dat', Rmin, Rmax, Zmin, Zmax, Bvac)

  call h5_deinit()

  deallocate(Ic, Bvac, Bn)
  call coil_deinit(coils(1))
  deallocate(coils)

  call print_ok
end program tilted_coil_fourier_modes
