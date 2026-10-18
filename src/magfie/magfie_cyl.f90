!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module magfie_cyl_mod
  implicit none
  private

  ! Values of the optional status argument of magfie_with_status.
  integer, parameter, public :: MAGFIE_OK = 0
  integer, parameter, public :: MAGFIE_OUTSIDE_LCFS = 1
  integer, parameter, public :: MAGFIE_OUTSIDE_BOUNDARY = 2

  public :: magfie_with_status

contains

subroutine magfie_with_status(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, status)
  ! Computes magnetic field module in units of the magnetic code  - bmod,
  ! square root of determinant of the metric tensor               - sqrtg,
  ! derivatives of the logarythm of the magnetic field module
  ! over coordinates                                              - bder,
  ! covariant componets of the unit vector of the magnetic
  ! field direction                                               - hcovar,
  ! contravariant components of this vector                       - hctrvr,
  ! contravariant component of the curl of this vector            - hcurl
  ! Order of coordinates is the following: x(1)=R (big radius),
  ! x(2)=phi (toroidal angle), x(3)=Z (altitude).
  !
  !  Input parameters:
  !            formal:  x                -    array of coordinates
  !  Output parameters:
  !            formal:  bmod
  !                     sqrtg
  !                     bder
  !                     hcovar
  !                     hctrvr
  !                     hcurl
  !
  !            optional status (intent(out)): MAGFIE_OK, MAGFIE_OUTSIDE_LCFS
  !                     (outside the LCFS while allow_sol is off) or
  !                     MAGFIE_OUTSIDE_BOUNDARY (outside the convex boundary).
  !                     Nothing is printed; ierrfield is set as before.
  !
  !  Called routines:  field

  use field_eq_mod, only : psi_axis,psi_sep,ierrfield,allow_sol
  use field_sub, only : psif
  use libneo_kinds, only : dp
  use field_sub,    only : field

  implicit none

  real(dp), intent(in) :: x(3)
  real(dp), intent(out) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
  integer, intent(out), optional :: status
  real(dp) :: hr,hf,hz

  real(dp) :: br,bf,bz, &
  BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ
  logical :: outside_boundary

  CALL field(x(1),x(2),x(3),br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ, &
      outside_boundary)

  ! The convex file bounds the interpolation domain.  SOL extension only
  ! relaxes the separatrix check; it must never make a point outside the
  ! interpolation boundary look valid after stretch_coords has clamped it.
  if (outside_boundary) then
    ierrfield=1
    if (present(status)) status = MAGFIE_OUTSIDE_BOUNDARY
  else if(allow_sol) then
    ierrfield=0
    if (present(status)) status = MAGFIE_OK
  else if((psif-psi_axis)/(psi_sep-psi_axis).gt.1d0) then
    ierrfield=1
    if (present(status)) status = MAGFIE_OUTSIDE_LCFS
  else
    ierrfield=0
    if (present(status)) status = MAGFIE_OK
  endif

  bmod = dsqrt(br**2 + bf**2 + bz**2)  ! B
  sqrtg = x(1)
  hr = br/bmod
  hf = bf/bmod
  hz = bz/bmod

  bder(1) = (brr*hr + bfr*hf + bzr*hz) / bmod
  bder(2) = (brf*hr + bff*hf + bzf*hz) / bmod
  bder(3) = (brz*hr + bfz*hf + bzz*hz) / bmod

  hcovar(1) = hr
  hcovar(2) = hf*x(1)
  hcovar(3) = hz

  hctrvr(1) = hr
  hctrvr(2) = hf/x(1)
  hctrvr(3) = hz

  ! hcurl = (curl vecB)/B and not curl vech = curl (vecB/B) (?)
  hcurl(1)=((bzf-x(1)*bfz)/bmod    + hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
  hcurl(2)=((brz-bzr)/bmod         + hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
  hcurl(3)=((bf+x(1)*bfr-brf)/bmod + hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
end subroutine magfie_with_status
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module magfie_cyl_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  ! Legacy external entry point (implicit-interface callers such as NEO-RT and
  ! POTATO); an optional argument cannot be added here without breaking them.
  ! Use magfie_with_status to obtain the outside-LCFS status.
  use libneo_kinds, only : dp
  use magfie_cyl_mod, only : magfie_with_status

  implicit none

  real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

  call magfie_with_status(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
end subroutine magfie
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
