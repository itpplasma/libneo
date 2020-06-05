!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! For testing: zero electric field
!
subroutine elefie(x, derphi)
  double precision, dimension(3) :: x,derphi
  derphi = 0d0
end subroutine elefie


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine elefie_mtprofile(x,derphi)
!
! Computes the derivatives of the electrostatic potential over
! coordinates (covariant electric field). Potential is normalized
! to reference temperature : phinor = e phi / T.
!
!   Input parameters:
!             formal:    x       -   array of coordinates
!
!   Output parameters:
!             formal:    derphi  -   array of derivatives
!

  USE polylag_3, only : mp,indef, plag1d
  use elefie_mod, only: Mtprofile, v0, Z1, am1, rbig, escale, bscale
  use neo_input, only : flux
  use constants, only: e_mass, e_charge, p_mass, c, pi
        
  implicit none
  
  integer :: ierr
  double precision, dimension(3) :: x,derphi
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: bder,hcovar,hctrvr,hcurl
  double precision :: r,phi,z,psi,phi_el,phi_el_pr,phi_el_prpr
  double precision :: qsafety

  integer,          dimension(mp) :: indu
  double precision, dimension(mp) :: xp,fp
  double precision :: s, der, dxm1
  double precision :: Mt, vth

  !   
  derphi=0.d0
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  qsafety = hctrvr(2)/hctrvr(3)
  
  s = x(1)
  
  dxm1=1.d0/(Mtprofile(2,1)-Mtprofile(1,1))
  call indef(s,Mtprofile(1,1),dxm1,size(Mtprofile,1),indu)
  
  xp=Mtprofile(indu,1)
  fp=Mtprofile(indu,2)
  call plag1d(s,fp,dxm1,xp,Mt,der)
  
  fp=Mtprofile(indu,3)
  call plag1d(s,fp,dxm1,xp,vth,der)

  Mt = Mt/bscale
  
  derphi(1) = -(1.0d8*flux*bscale/(2.0d0*pi))*Mt*vth/(qsafety*c*rbig)
  derphi(1) = derphi(1)*Z1*e_charge/(am1*p_mass*v0**2/2.d0) ! normalization

  derphi(1) = -escale*derphi(1) ! TODO check if this is correct sign
!  derphi(1) = 0.1d0 ! TODO remove this test
  return
end subroutine elefie_mtprofile
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine velo(tau,z,vz)
!
!
!  Computes the components of the 5-D velocity          -  vz
!  for given set of phase space coordinates             -  z.
!
!  Warning !!!  The dimensionless time is used (see below)
!
!  Order of the coordinates is the following:
!  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
!                                 dimension of the general covariant
!                                 space coordinate system
!  z(4) = p                     - momentum  module normalized to
!                                 thermal momentum and sqrt(2);
!  z(5) = alambd                - cosine of the pitch-angle
!
!  Input parameters:
!            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
!                     z      -   see above
!            common:  rmu    -   inverse relativistic temperature
!                     ro0    -   Larmor radius for the reference
!                                magnetic field and temperature:
!                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
!                     p0     -   dimensionless momentum module in the
!                                initial point
!                     alamb0 -   cos of pitch-angle in the initial point
!  Output parameters:
!            formal:  vz     -   see above
!
!  Called routines: magfie, elefie
!
      use parmot_mod, only : rmu,ro0
      use magfield_mod, only : ierrfield
!
      implicit none
!
      integer :: i
!
      double precision, intent(in)  :: tau, z
      double precision, intent(out) :: vz
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision derphi
      double precision p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
      double precision rmumag,rovsqg,rosqgb,rovbm
      double precision a_phi,a_b,a_c,hstar
      double precision s_hc,hpstar,phidot,blodot,bra
      double precision pardeb
!
      dimension z(5), vz(5)
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
      dimension derphi(3)
      dimension a_phi(3),a_b(3),a_c(3),hstar(3)
!
      do 1 i=1,3
        x(i)=z(i)
 1    continue
!
! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
!            bmod   - dimensionless magnetic field module: bmod=B/B_ref
!            sqrtg  - Jacobian of space coordinates (square root of
!                     metric tensor
!            bder   - derivatives of logarithm of bmod over space coords 
!                     (covariant vector)
!            hcovar - covariant components of the unit vector along
!                     the magnetic field
!            hctrvr - contravariant components of the unit vector along
!                     the magnetic field
!            hcurl  - contravariant components of the curl of this vector
!
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      if(ierrfield.ne.0) then
        vz=0.d0
        return
      endif
! in elefie: x(i)   - space coords (input, see above)
!            derphi - derivatives of the dimensionless electric potential
!                     phihat=e*phi/T over space coords (covar. vector)
!
      call elefie(x,derphi)
!
      p=z(4)
      alambd=z(5)
!
      p2=p*p
      ovmu=2.d0/rmu
      gamma2=p2*ovmu+1.d0
      gamma=dsqrt(gamma2)
      ppar=p*alambd
! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
      vpa=ppar/gamma
      coala=(1.d0-alambd**2)
! rmumag - magnetic moment
      rmumag=.5d0*p2*coala/bmod
!
      rovsqg=ro0/sqrtg
      rosqgb=.5d0*rovsqg/bmod
      rovbm=ro0/bmod
!
      a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
      a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
      a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
      a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
      a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
      a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg
!
      s_hc=0.d0
      do i=1,3
        a_c(i)=hcurl(i)*rovbm
        s_hc=s_hc+a_c(i)*hcovar(i)
        hstar(i)=hctrvr(i)+ppar*a_c(i)
      enddo
      hpstar=1.d0+ppar*s_hc
!
! velocities in the coordinate space
!
! phidot - derivative of the dmls el. potential over dmls time
! blodot - derivative of the logarith of the mag. field module over dmls time
      phidot=0.d0
      blodot=0.d0
      do i=1,3
        bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
        vz(i)=bra/hpstar
        phidot=phidot+vz(i)*derphi(i)
        blodot=blodot+vz(i)*bder(i)
      enddo
!
! velocities in the phase space
!
      vz(4)=-0.5d0*gamma*phidot/p
!      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
      vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
            + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine orbit_timestep(z,dtau,dtaumin,ierr)
!
! collisions
      use collis_alp, only : swcoll,iswmod
! end collisions
      use magfield_mod, only : ierrfield
!
      implicit none
!
      integer, parameter          :: ndim=5, nstepmax=100000000
      double precision, parameter :: relerr=1d-6
!
      integer :: ierr,j
      double precision :: dtau,dtaumin,phi,tau1,tau2
!
      double precision, dimension(2)    :: y
      double precision, dimension(ndim) :: z
! 24.08.2011
      double precision, parameter :: vdr_dv=0.03d0
! 24.08.2011 end
! collisions
      integer :: ierrcol
      double precision :: dtauc
! end collisions
!
      external velo
!
      if(dtaumin*nstepmax.le.dtau) then
        ierr=2
        print *,'orbit_timestep: number of steps exceeds nstepmax'
        return
      endif
!
      ierr=0
      y(1)=z(1)
      phi=z(2)
      y(2)=z(3)
!
      call chamb(y,phi,ierr)
!
      if(ierr.ne.0) return
      tau1=0.d0
! 24.08.2011      tau2=dtaumin
      tau2=dtaumin/(dabs(z(5))+vdr_dv)
!
      do while(tau2.lt.dtau)
!
        call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)
!
! collisions
        if(swcoll) then
          dtauc=tau2-tau1
!
          call stost(z,dtauc,iswmod,ierrcol)
!
        endif
! end collisions
!
        y(1)=z(1)
        phi=z(2)
        y(2)=z(3)
!
        call chamb(y,phi,ierr)
        ierr=max(ierr,ierrfield)
!
        if(ierr.ne.0) return
        tau1=tau2
! 24.08.2011        tau2=tau2+dtaumin
        tau2=tau2+dtaumin/(dabs(z(5))+vdr_dv)
      enddo
!
      tau2=dtau
!
      call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)
!
! collisions
      if(swcoll) then
        dtauc=tau2-tau1
!
        call stost(z,dtauc,iswmod,ierrcol)
!
      endif
! end collisions
!
      y(1)=z(1)
      phi=z(2)
      y(2)=z(3)
!
      call chamb(y,phi,ierr)
      ierr=max(ierr,ierrfield)
!
      if(ierr.ne.0) return
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhs_mflint(phi,y,dery)
!
! Computes the right hand side of the magnetic field line equation for
! the integration over the toroidal angle, subintegrand for the flux tube 
! volume $1/B^\varphi$ and subintegrants for Boozer $B_{00}$ computation
! $B^2/B^\varphi$ and $B^3/B^\varphi$                       -   dery
!
! Oder of coordinates in magfie is the following: x(1)=R (big radius),
! x(2)=phi (toroidal angle), x(3)=z (altitude).
!
!  Input parameters:
!            formal:
!                    y(1:2) - coordinates in the poloidal plane (phi=const)
!                    y(3:5) - integrals
!  Output parameters:
!            formal:
!                 dery(1:5) - vector of the right-hand side
!  Called routines:  magfie
!
      use magfield_mod, only : ierrfield
!
      double precision :: phi
      double precision, dimension(5) :: y,dery
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      x(1)=y(1)
      x(2)=phi
      x(3)=y(2)
!
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      if(ierrfield.ne.0) then
        dery=0.d0
        return
      endif
      dery(1)=hctrvr(1)/hctrvr(2)
      dery(2)=hctrvr(3)/hctrvr(2)
      dery(3)=1.d0/(bmod*hctrvr(2))
      dery(4)=bmod/hctrvr(2)
      dery(5)=bmod**2/hctrvr(2)
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 01.09.2011      subroutine integrate_mfl(npoi,dphi,rbeg,phibeg,zbeg,         &
      subroutine integrate_mfl(npoi,dphii,rbeg,phibeg,zbeg,         &
                               xstart,bstart,volstart,bmod00,ierr)
!
      use magfield_mod, only : ierrfield
!
      implicit none
!
      integer, parameter          :: ndim=5
      double precision, parameter :: relerr=1d-8
      integer :: npoi,i,ierr
      double precision :: dphi,rbeg,phibeg,zbeg,bmod00,phi,phiold
      double precision, dimension(3,npoi) :: xstart
      double precision, dimension(npoi)   :: bstart,volstart
      double precision, dimension(ndim)   :: y
!
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
! 01.09.2011
! 03.09.2011      integer, parameter          :: ndphi=10
! 08.09.2011      integer, parameter          :: ndphi=4
! 08.09.2011      integer, parameter          :: ndphi=3
! 10.09.2011      integer, parameter          :: ndphi=5
      integer, parameter          :: ndphi=13
      integer          :: i2
      double precision :: dphii
! 01.09.2011 end
!
      external :: rhs_mflint
! 01.09.2011
      dphi=dphii/ndphi
! 01.09.2011 end!
      ierr=0
!
      phi=phibeg
      y(1)=rbeg
      y(2)=zbeg
      y(3:5)=0.d0
!
      call chamb(y(1:2),phi,ierr)
!
      x(1)=y(1)
      x(2)=phi
      x(3)=y(2)
!
      call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      ierr=max(ierr,ierrfield)
      if(ierr.ne.0) return
      xstart(:,1)=x
      bstart(1)=bmod
      volstart(1)=y(3)
!
      do i=2,npoi
! 01.09.2011 
      do i2=1,ndphi
! 01.09.2011 end
        phiold=phi
        phi=phiold+dphi
!
        call odeint_allroutines(y,ndim,phiold,phi,relerr,rhs_mflint)
        if(ierrfield.ne.0) exit
! 01.09.2011 
      enddo
! 01.09.2011 end
!
        call chamb(y(1:2),phi,ierr)
!
        x(1)=y(1)
        x(2)=phi
        x(3)=y(2)
!
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        ierr=max(ierr,ierrfield)
        if(ierr.ne.0) return
!
        xstart(:,i)=x
        bstart(i)=bmod
        volstart(i)=y(3)
      enddo
!
      volstart=volstart/volstart(npoi)
      bmod00=y(5)/y(4)
!
      end subroutine integrate_mfl
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine binsrc(p,nmin,nmax,xi,i)
!
! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.
!
  implicit none
!
  integer                                :: n,nmin,nmax,i,imin,imax,k
  double precision                       :: xi
  double precision, dimension(nmin:nmax) :: p
!
  imin=nmin
  imax=nmax
  n=nmax-nmin
!
  do k=1,n
    i=(imax-imin)/2+imin
    if(p(i).gt.xi) then
      imax=i
    else
      imin=i
    endif
    if(imax.eq.imin+1) exit
  enddo
!
  i=imax
!
  return
  end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine regst(z,dtau,ierr)
!
      implicit none
!
      integer, parameter          :: ndim=5
      double precision, parameter :: relerr=1d-6
!
      integer :: ierr,j
      double precision :: dtau,tau1,tau2
!
      double precision, dimension(ndim) :: z
!
      external velo
!
      ierr=0
!
      tau1=0.d0
      tau2=dtau
!
      call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)
!
      if(z(1).gt.1.d0) ierr=1
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
