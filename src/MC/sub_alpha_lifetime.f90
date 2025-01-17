!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! For testing: zero electric field
subroutine elefie(x, derphi)
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind), dimension(3) :: x,derphi
  derphi = 0d0
end subroutine elefie


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine velo(tau,z,vz)
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

  use libneo_kinds, only : real_kind
      use parmot_mod, only : rmu,ro0

      implicit none

      integer :: i

  real(kind=real_kind), intent(in)  :: tau, z
  real(kind=real_kind), intent(out) :: vz
  real(kind=real_kind) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  real(kind=real_kind) derphi
  real(kind=real_kind) p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
  real(kind=real_kind) rmumag,rovsqg,rosqgb,rovbm
  real(kind=real_kind) a_phi,a_b,a_c,hstar
  real(kind=real_kind) s_hc,hpstar,phidot,blodot,bra
  real(kind=real_kind) pardeb

      dimension z(5), vz(5)
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
      dimension derphi(3)
      dimension a_phi(3),a_b(3),a_c(3),hstar(3)

      do 1 i=1,3
        x(i)=z(i)
 1    continue

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

  ! TODO: error handling magfie
  !      if(ierrfield.ne.0) then
  !        vz=0.d0
  !        return
  !      endif
  ! in elefie: x(i)   - space coords (input, see above)
  !            derphi - derivatives of the dimensionless electric potential
  !                     phihat=e*phi/T over space coords (covar. vector)

      call elefie(x,derphi)

      p=z(4)
      alambd=z(5)

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

      rovsqg=ro0/sqrtg
      rosqgb=.5d0*rovsqg/bmod
      rovbm=ro0/bmod

      a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
      a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
      a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
      a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
      a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
      a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg

      s_hc=0.d0
      do i=1,3
        a_c(i)=hcurl(i)*rovbm
        s_hc=s_hc+a_c(i)*hcovar(i)
        hstar(i)=hctrvr(i)+ppar*a_c(i)
      enddo
      hpstar=1.d0+ppar*s_hc

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

  ! velocities in the phase space

      vz(4)=-0.5d0*gamma*phidot/p
      vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
            + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))

end subroutine velo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine orbit_timestep(z,dtau,dtaumin,ierr)
  ! collisions
  use collis_alp, only : swcoll,iswmod
  ! end collisions
  use libneo_kinds, only : real_kind
  use odeint_sub, only : odeint_allroutines

  implicit none

  integer, parameter          :: ndim=5, nstepmax=100000000
  real(kind=real_kind), parameter :: relerr=1d-6
  real(kind=real_kind), parameter :: vdr_dv=0.03d0

  integer :: ierr
  real(kind=real_kind) :: dtau,dtaumin,phi,tau1,tau2

  real(kind=real_kind), dimension(2)    :: y
  real(kind=real_kind), dimension(ndim) :: z

  ! collisions
  integer :: ierrcol
  real(kind=real_kind) :: dtauc
  ! end collisions

  external velo

  if(dtaumin*nstepmax.le.dtau) then
    ierr=2
    print *,'orbit_timestep: number of steps exceeds nstepmax'
    return
  endif

  ierr=0
  y(1)=z(1)
  phi=z(2)
  y(2)=z(3)

  call chamb(y,phi,ierr)

  if(ierr.ne.0) return
  tau1=0.d0
  tau2=dtaumin/(dabs(z(5))+vdr_dv)

  do while(tau2.lt.dtau)

    call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)

    ! collisions
    if (swcoll) then
      dtauc=tau2-tau1

      call stost(z,dtauc,iswmod,ierrcol)

    endif
    ! end collisions

    y(1)=z(1)
    phi=z(2)
    y(2)=z(3)

    call chamb(y,phi,ierr)
    ! TODO: error handling magfie
    if(ierr.ne.0) return
    tau1=tau2
    tau2=tau2+dtaumin/(dabs(z(5))+vdr_dv)
  end do

  tau2=dtau

  call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)

  ! collisions
  if (swcoll) then
    dtauc=tau2-tau1

    call stost(z,dtauc,iswmod,ierrcol)

  end if
  ! end collisions

  y(1)=z(1)
  phi=z(2)
  y(2)=z(3)

  call chamb(y,phi,ierr)
  ! TODO: error handling magfie
  if(ierr.ne.0) return

end subroutine orbit_timestep

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rhs_mflint(phi,y,dery)
  ! Computes the right hand side of the magnetic field line equation for
  ! the integration over the toroidal angle, subintegrand for the flux tube
  ! volume $1/B^\varphi$ and subintegrants for Boozer $B_{00}$ computation
  ! $B^2/B^\varphi$ and $B^3/B^\varphi$                       -   dery
  !
  ! Order of coordinates in magfie is the following: x(1)=R (big radius),
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

  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind) :: phi
  real(kind=real_kind), dimension(5) :: y,dery
  real(kind=real_kind) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)

  x(1)=y(1)
  x(2)=phi
  x(3)=y(2)

  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)


  ! TODO: error handling magfie
  dery(1)=hctrvr(1)/hctrvr(2)
  dery(2)=hctrvr(3)/hctrvr(2)
  dery(3)=1.d0/(bmod*hctrvr(2))
  dery(4)=bmod/hctrvr(2)
  dery(5)=bmod**2/hctrvr(2)

end subroutine rhs_mflint

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine integrate_mfl(npoi,dphii,rbeg,phibeg,zbeg,         &
                               xstart,bstart,volstart,bmod00,ierr)
  use libneo_kinds, only : real_kind
  use odeint_sub, only : odeint_allroutines

  implicit none

  integer, parameter :: ndim=5
  integer, parameter :: ndphi=13
  real(kind=real_kind), parameter :: relerr=1d-8

  integer :: npoi,i,ierr
  real(kind=real_kind) :: dphi,rbeg,phibeg,zbeg,bmod00,phi,phiold
  real(kind=real_kind), dimension(3,npoi) :: xstart
  real(kind=real_kind), dimension(npoi)   :: bstart,volstart
  real(kind=real_kind), dimension(ndim)   :: y

  real(kind=real_kind) x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
  integer          :: i2
  real(kind=real_kind) :: dphii

  external :: rhs_mflint

  dphi=dphii/ndphi

  ierr=0

  phi=phibeg
  y(1)=rbeg
  y(2)=zbeg
  y(3:5)=0.d0

  call chamb(y(1:2),phi,ierr)

  x(1)=y(1)
  x(2)=phi
  x(3)=y(2)

  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)

  ! TODO: error handling magfie
  if(ierr.ne.0) return
  xstart(:,1)=x
  bstart(1)=bmod
  volstart(1)=y(3)

  do i=2,npoi

    do i2=1,ndphi
      phiold=phi
      phi=phiold+dphi

      call odeint_allroutines(y,ndim,phiold,phi,relerr,rhs_mflint)
      ! TODO: error handling magfie

    end do

    call chamb(y(1:2),phi,ierr)

    x(1)=y(1)
    x(2)=phi
    x(3)=y(2)

    call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    ! TODO: error handling magfie
    if(ierr.ne.0) return

    xstart(:,i)=x
    bstart(i)=bmod
    volstart(i)=y(3)
  end do

  volstart=volstart/volstart(npoi)
  bmod00=y(5)/y(4)

end subroutine integrate_mfl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine binsrc(p,nmin,nmax,xi,i)
  ! Finds the index  i  of the array of increasing numbers   p  with dimension  n
  ! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.

  use libneo_kinds, only : real_kind

  implicit none

  integer                                :: n,nmin,nmax,i,imin,imax,k
  real(kind=real_kind)                       :: xi
  real(kind=real_kind), dimension(nmin:nmax) :: p

  imin=nmin
  imax=nmax
  n=nmax-nmin

  do k=1,n
    i=(imax-imin)/2+imin
    if(p(i).gt.xi) then
      imax=i
    else
      imin=i
    endif
    if(imax.eq.imin+1) exit
  end do

  i=imax

end subroutine binsrc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine regst(z,dtau,ierr)
  use libneo_kinds, only : real_kind
  use odeint_sub, only : odeint_allroutines

  implicit none

  integer, parameter          :: ndim=5
  real(kind=real_kind), parameter :: relerr=1d-6

  integer :: ierr
  real(kind=real_kind) :: dtau,tau1,tau2

  real(kind=real_kind), dimension(ndim) :: z

  external velo

  ierr=0

  tau1=0.d0
  tau2=dtau

  call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo)

  if(z(1).gt.1.d0) ierr=1

end subroutine regst
