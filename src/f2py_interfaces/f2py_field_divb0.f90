subroutine field_divb0_initialize_from_grid(nr_in, np_in, nz_in, ntor_in, &
    rmin_in, rmax_in, pmin_in, pmax_in, zmin_in, zmax_in, br_in, bp_in, bz_in)
  use libneo_kinds, only : dp
  use field_c_mod, only : icall_c, ntor, nr, np, nz, icftype, rmin, rmax, pmin, pmax, zmin, zmax
  implicit none

  interface
    subroutine vector_potentials(nr_in,np_in,nz_in,ntor_in,      &
             rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in,  &
             br,bp,bz)
      use libneo_kinds, only : dp
      implicit none
      integer, intent(in) :: nr_in, np_in, nz_in, ntor_in
      real(dp), intent(in) :: rmin_in, rmax_in, pmin_in, pmax_in, zmin_in, zmax_in
      real(dp), intent(inout) :: br(nr_in, np_in, nz_in)
      real(dp), intent(inout) :: bp(nr_in, np_in, nz_in)
      real(dp), intent(inout) :: bz(nr_in, np_in, nz_in)
    end subroutine vector_potentials
  end interface

  integer, intent(in) :: nr_in, np_in, nz_in, ntor_in
  real(dp), intent(in) :: rmin_in, rmax_in, pmin_in, pmax_in, zmin_in, zmax_in
  real(dp), intent(in) :: br_in(nr_in, np_in, nz_in)
  real(dp), intent(in) :: bp_in(nr_in, np_in, nz_in)
  real(dp), intent(in) :: bz_in(nr_in, np_in, nz_in)

  real(dp), allocatable :: br(:,:,:), bp(:,:,:), bz(:,:,:)

  allocate(br(nr_in, np_in, nz_in))
  allocate(bp(nr_in, np_in, nz_in))
  allocate(bz(nr_in, np_in, nz_in))

  br = br_in
  bp = bp_in
  bz = bz_in

  call vector_potentials(nr_in, np_in, nz_in, ntor_in, &
       rmin_in, rmax_in, pmin_in, pmax_in, zmin_in, zmax_in, br, bp, bz)

  nr = nr_in
  np = np_in
  nz = nz_in
  ntor = ntor_in
  rmin = rmin_in
  rmax = rmax_in
  pmin = pmin_in
  pmax = pmax_in
  zmin = zmin_in
  zmax = zmax_in
  icftype = 4
  icall_c = 1

  deallocate(br, bp, bz)

end subroutine field_divb0_initialize_from_grid


subroutine field_divb0_eval(n_points, r, phi, z, br, bp, bz)
  use libneo_kinds, only : dp
  implicit none

  interface
    subroutine field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
      use libneo_kinds, only : dp
      implicit none
      real(dp), intent(in) :: r, phi, z
      real(dp), intent(out) :: Br, Bp, Bz
      real(dp), intent(out) :: dBrdR, dBrdp, dBrdZ
      real(dp), intent(out) :: dBpdR, dBpdp, dBpdZ
      real(dp), intent(out) :: dBzdR, dBzdp, dBzdZ
    end subroutine field_divfree
  end interface

  integer, intent(in) :: n_points
  real(dp), intent(in) :: r(n_points), phi(n_points), z(n_points)
  real(dp), intent(out) :: br(n_points), bp(n_points), bz(n_points)

  real(dp) :: dBrdR, dBrdp, dBrdZ
  real(dp) :: dBpdR, dBpdp, dBpdZ
  real(dp) :: dBzdR, dBzdp, dBzdZ
  integer :: i

  do i = 1, n_points
    call field_divfree(r(i), phi(i), z(i), br(i), bp(i), bz(i), &
         dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  end do

end subroutine field_divb0_eval
