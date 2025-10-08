program test_ntor0_tilted_coil
  use util_for_test, only: print_test, print_ok, print_fail
  use coil_tools, only: coil_t, coil_init, coil_deinit, &
                       vector_potential_biot_savart_fourier, write_Anvac_fourier
  use libneo_kinds, only: dp
  implicit none

  type(coil_t), allocatable :: coils(:)
  complex(dp), allocatable :: AnR(:,:,:,:), Anphi(:,:,:,:), AnZ(:,:,:,:)
  complex(dp), allocatable :: dAnphi_dR(:,:,:,:), dAnphi_dZ(:,:,:,:)
  real(dp) :: Rmin, Rmax, Zmin, Zmax
  integer :: nR, nZ, nphi, nmax, nseg, k, iostat
  real(dp) :: min_distance, max_eccentricity
  logical :: use_convex_wall
  real(dp) :: x, y, z

  call print_test("ntor=0 field reconstruction for tilted coil")

  ! Grid parameters
  Rmin = 0.5_dp
  Rmax = 4.0_dp
  Zmin = -1.5_dp
  Zmax = 2.5_dp
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
    coils(1)%XYZ(1, k) = x
    coils(1)%XYZ(2, k) = y
    coils(1)%XYZ(3, k) = z
  end do
  close(10)

  ! Compute vector potential Fourier modes
  call vector_potential_biot_savart_fourier(coils, nmax, min_distance, max_eccentricity, use_convex_wall, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, &
    AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)

  ! Write to NetCDF file
  call write_Anvac_fourier('tilted_coil_Anvac.nc', 1, nmax, &
    Rmin, Rmax, Zmin, Zmax, nR, nphi, nZ, &
    AnR, Anphi, AnZ, dAnphi_dR, dAnphi_dZ)

  call coil_deinit(coils(1))
  deallocate(coils)

  call print_ok
end program test_ntor0_tilted_coil
