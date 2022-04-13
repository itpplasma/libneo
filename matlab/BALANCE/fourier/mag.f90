!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!
! Computes magnetic field module normalized to axial value  - bmod,
! square root of determinant of the metric tensor           - sqrtg,
! derivatives of the logarythm of the magnetic field module
! over coordinates                                          - bder,
! covariant componets of the unit vector of the magnetic
! field direction                                           - hcovar,
! contravariant components of this vector                   - hctrvr,
! derivatives of the covariant components of this vector    - hcoder(i,j)
! here hcoder(i,j) is the derivative of hcovar(j) over x(i)
! for given set of coordinates x(i).
! Oder of coordinates is the following: x(1)=R (big radius), 
! x(2)=phi (toroidal angle), x(3)=z (altitude).
!
!  Input parameters:
!            formal:  x                -    array of coordinates
!  Output parameters:
!            formal:  bmod
!                     sqrtg
!                     bder
!                     hcovar
!                     hctrvr
!                     hcoder
!
!  Called routines:  field_eq
!
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
      double precision hr,hf,hz
!
      double precision ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ, &
      BRK,BZK,BRRK,BRZK,BZRK,BZZK
!
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!
      rbig=max(x(1),1d-12)
!
!cccccc computation of gb in cylindrical co-ordinates cccccccc
      ri=rbig
      fii=x(2)
      zi=x(3)

      CALL field_eq(ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)

!ccccc end of gb computation cccccccccc
      bmod=dsqrt(br**2+bf**2+bz**2)
      sqrtg=rbig
      hr=br/bmod
      hf=bf/bmod
      hz=bz/bmod
!
      bder(1)=(brr*hr+bfr*hf+bzr*hz)/bmod
      bder(2)=(brf*hr+bff*hf+bzf*hz)/bmod
      bder(3)=(brz*hr+bfz*hf+bzz*hz)/bmod
!
      hcovar(1)=hr
      hcovar(2)=hf*rbig
      hcovar(3)=hz
!
      hctrvr(1)=hr
      hctrvr(2)=hf/rbig
      hctrvr(3)=hz
!
      hcoder(1,1)=brr/bmod-hcovar(1)*bder(1)
      hcoder(2,1)=brf/bmod-hcovar(1)*bder(2)
      hcoder(3,1)=brz/bmod-hcovar(1)*bder(3)
      hcoder(1,2)=(rbig*bfr+bf)/bmod-hcovar(2)*bder(1)
      hcoder(2,2)=rbig*bff/bmod-hcovar(2)*bder(2)
      hcoder(3,2)=rbig*bfz/bmod-hcovar(2)*bder(3)
      hcoder(1,3)=bzr/bmod-hcovar(3)*bder(1)
      hcoder(2,3)=bzf/bmod-hcovar(3)*bder(2)
      hcoder(3,3)=bzz/bmod-hcovar(3)*bder(3)
!
      hctder(1,1)=brr/bmod-hctrvr(1)*bder(1)
      hctder(2,1)=brf/bmod-hctrvr(1)*bder(2)
      hctder(3,1)=brz/bmod-hctrvr(1)*bder(3)
      hctder(1,2)=(bfr-bf/rbig)/(rbig*bmod)-hctrvr(2)*bder(1)
      hctder(2,2)=bff/(rbig*bmod)-hctrvr(2)*bder(2)
      hctder(3,2)=bfz/(rbig*bmod)-hctrvr(2)*bder(3)
      hctder(1,3)=bzr/bmod-hctrvr(3)*bder(1)
      hctder(2,3)=bzf/bmod-hctrvr(3)*bder(2)
      hctder(3,3)=bzz/bmod-hctrvr(3)*bder(3)
!
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
