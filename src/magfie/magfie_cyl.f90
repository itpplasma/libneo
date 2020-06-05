!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
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
!  Called routines:  GBhs,GBRZd 
!
      use magfield_mod, only : ierrfield
!
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision hr,hf,hz
!
      double precision ri,fii,zi,br,bf,bz, &
      BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ, &
      BRK,BZK,BRRK,BRZK,BZRK,BZZK
!
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      rbig=max(x(1),1d-12)
!
!cccccc computation of gb in cylindrical co-ordinates cccccccc
      ri=rbig
      fii=x(2)
      zi=x(3)

! 19.03.2010 CALL GBhs(ri,fii,zi,br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
      CALL gbtj2(ri,fii,zi,br,bf,bz,BRR,BRF,BRZ,BFR,BFF,BFZ,BZR,BZF,BZZ)
!
      if(ierrfield.eq.1) return
! 18.03.2010      CALL GBRZd(ri,zi,BRK,BZK,BRRK,BRZK,BZRK,BZZK)
!
! Here we add the vertical field
! 18.03.2010      br=br+BRK
! 18.03.2010      bz=bz+BZK
!
! 18.03.2010      BRR=BRR+BRRK
! 18.03.2010      BRZ=BRZ+BRZK
! 18.03.2010      BZR=BZR+BZRK
! 18.03.2010      BZZ=BZZ+BZZK
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
      hcurl(1)=((bzf-rbig*bfz)/bmod                        &
     +          hcovar(2)*bder(3)-hcovar(3)*bder(2))/sqrtg
      hcurl(2)=((brz-bzr)/bmod                             &
     +          hcovar(3)*bder(1)-hcovar(1)*bder(3))/sqrtg
      hcurl(3)=((bf+rbig*bfr-brf)/bmod                     &
     +          hcovar(1)*bder(2)-hcovar(2)*bder(1))/sqrtg
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
