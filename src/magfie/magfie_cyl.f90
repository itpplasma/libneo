!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
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
  !  Called routines:  field

  use field_eq_mod, only : psi_axis,psi_sep,psif,ierrfield
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind) :: x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  real(kind=real_kind) :: hr,hf,hz

  real(kind=real_kind) :: br,bf,bz, &
  BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ

  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)

  CALL field(x(1),x(2),x(3),br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

  if((psif-psi_axis)/(psi_sep-psi_axis).gt.1d0) then
    print *,'magfie: point is outside separatrix, (R,Z) = ',x(1),x(3)
    ierrfield=1
  else
    ierrfield=0
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
  return
end subroutine magfie
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
