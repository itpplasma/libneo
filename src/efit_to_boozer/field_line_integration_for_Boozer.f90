!
  subroutine flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta,           &
                              rmn,rmx,zmn,zmx,raxis,zaxis,sigma,      &
                              rbeg,rsmall,qsaf,psi_pol,               &
                              psi_tor_vac,psi_tor_plas,C_const,       &
                              R_ts,Z_ts,bmod_ts,sqgnorm_ts,Gfunc_ts)
!
  use field_sub
  use odeint_allroutines_sub
  use rhs_boozer_sub, only: rhs_axis, rhs_surf
  use field_eq_mod,  only : icall_eq,nrad,nzet,rad,zet,rtf,btf
  use rhs_surf_mod, only : dz_dphi
  use efit_to_boozer_mod, only : psimax
!
  implicit none
!
  integer, parameter :: neq_ax=4,neq=5
  integer, parameter :: niter_axis=20  !number of iterations for finding axis
  integer, parameter :: niter=10       !number of iterations for Newton method
  integer, parameter :: nstep_min=10   !minimum number of steps
!
  integer :: nstep,nsurfmax,nlabel,ntheta
  integer :: i,j,nsurf,nmap,isurf,iter
!
  double precision, parameter :: pi = 3.14159265358979d0
  double precision, dimension(4), parameter :: win   = (/-1.d0, 13.d0, 13.d0, -1.d0 /) / 24.d0
  double precision, dimension(4), parameter :: wleft = (/ 9.d0, 19.d0, -5.d0,  1.d0 /) / 24.d0
  double precision :: rmn,rmx,zmn,zmx,raxis,zaxis
  double precision :: relerr,phiout
  double precision :: phi,rrr,ppp,zzz
  double precision :: aiota,hbr
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: psi_axis,h,sig,phi_sep,sigma
!
  double precision, dimension(neq_ax)        :: yax
  double precision, dimension(neq)           :: ymet
  double precision, dimension(nlabel)        :: rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor_plas,C_const
  double precision, dimension(ntheta,nlabel) :: R_ts,Z_ts,bmod_ts,sqgnorm_ts,Gfunc_ts
!
  double precision, dimension(:), allocatable :: dpsitor_dR

  associate(dummy => nstep)
  end associate
!
!
  nmap=10         !number of maps for finding magnetic axis
!
! Initialization of the field:
!
  icall_eq=-1
  rrr=1.d0
  ppp=0.d0
  zzz=0.d0
!
  call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! End of initialization
!
! Computation box:
  rmn=rad(1)
  rmx=rad(nrad)
  zmn=zet(1)
  zmx=zet(nzet)
!
  rrr=0.5d0*(rmn+rmx)
  ppp=0.d0
  zzz=0.5d0*(zmn+zmx)
!
! Search for the magnetic axis
!
  relerr = 1.d-9
  phi=0.d0
  phiout=2.d0*pi
  yax=0.d0
  yax(1)=rrr
  yax(2)=zzz
  do iter=1,niter_axis
    yax(3:4)=0.d0
    do i=1,nmap
      call odeint_allroutines(yax,neq_ax,phi,phiout,relerr,rhs_axis)
    enddo
    yax(1:2)=yax(3:4)/(phiout-phi)/dfloat(nmap)
  enddo
  raxis=yax(1)
  zaxis=yax(2)
!
  call field_eq(raxis,ppp,zaxis,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  psi_axis=psif
  print *,'toroidal field = ',btf
  open(1,file='btor_rbig.dat')
  write (1,*) btf,rtf
  close(1)
!
! End of search for the magnetic axis
!
  hbr=(rmx-raxis)/nsurfmax
!
  rrr=raxis+hbr
  zzz=zaxis
  call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Direction of poloidal field:
  sigma=sign(1.d0,Bz*Bp)
!
!-------------------------------------------------------------------------------
!
! Scan of flux surfaces
!
  h=2.d0*pi/nstep_min
!
  nsurf = 0
  surf: do isurf=1,nsurfmax
    phi=0.d0
    phiout=h
    phi_sep=phi
    ymet(1)=raxis+hbr*isurf
    ymet(2)=zaxis
    ymet(3)=0.d0
    ymet(4)=0.d0
    ymet(5)=0.d0
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
    phi_sep=phi_sep+phiout
    sig=ymet(2)-zaxis
    do while(sig*(ymet(2)-zaxis).gt.0.d0)
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
      if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
          ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
        nsurf=isurf-1
        exit surf
      endif
    enddo
    sig=ymet(2)-zaxis
    do while(sig*(ymet(2)-zaxis).gt.0.d0)
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
      if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
          ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
        nsurf=isurf-1
        exit surf
      endif
      if (psif > psimax) then
        nsurf=isurf-1
        exit surf
      endif
    enddo
!
! Newton method
    do iter=1,niter
      phiout=(zaxis-ymet(2))/dz_dphi
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
    enddo
!
  enddo surf
!
  nsurf=nsurf-1  !last point is bad, remove it
!
  print *,'Field line integration: separatrix found'
!
!------------------------------------------------------------------------------

! Re-define start points step size in R for data storage
  hbr=hbr*dfloat(nsurf)/dfloat(nlabel)
!
!------------------------------------------------------------------------------
!
! Computation of flux functions: effective radius, safety factor, poloidal and toroidal fluxes
!
  allocate(dpsitor_dR(0:nlabel))
  dpsitor_dR(0)=0.d0
!
  do isurf=1,nlabel
    phi=0.d0
    phiout=h
    phi_sep=phi
    rbeg(isurf)=raxis+hbr*dfloat(isurf)
    ymet(1)=rbeg(isurf)
    ymet(2)=zaxis
    ymet(3)=0.d0
    ymet(4)=0.d0
    ymet(5)=0.d0
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
    phi_sep=phi_sep+phiout
    sig=ymet(2)-zaxis
    do while(sig*(ymet(2)-zaxis).gt.0.d0)
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
    enddo
    sig=ymet(2)-zaxis
    do while(sig*(ymet(2)-zaxis).gt.0.d0)
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
    enddo
! Newton method
    do iter=1,niter
      phiout=(zaxis-ymet(2))/dz_dphi
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
      phi_sep=phi_sep+phiout
    enddo
!
    aiota=2.d0*pi/phi_sep
    rrr=ymet(1)
    zzz=ymet(2)
!
    call  field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    rsmall(isurf)=sqrt(abs(ymet(3))/pi)
    qsaf(isurf)=1.d0/aiota
    psi_pol(isurf)=psif-psi_axis
    psi_tor_vac(isurf)=-ymet(4)/(2.d0*pi)
    C_const(isurf)=ymet(5)*aiota/(2.d0*pi)
!
    dpsitor_dR(isurf)=qsaf(isurf)*rrr*Bz
  enddo
!
  psi_tor_plas(1)=sum(dpsitor_dR(0:3)*wleft)*hbr
  do isurf=2,nlabel-1
    psi_tor_plas(isurf)=psi_tor_plas(isurf-1)+sum(dpsitor_dR(isurf-2:isurf+1)*win)*hbr
  enddo
  psi_tor_plas(nlabel)=psi_tor_plas(nlabel-1)+sum(dpsitor_dR(nlabel:nlabel-3:-1)*wleft)*hbr
!
  deallocate(dpsitor_dR)
!
  print *,'Normalized poloidal flux psi at the LCMS = ',psi_pol(nlabel),'(psi on axis = 0)'
!
  print *,'Field line integration: flux functions done'
!
!------------------------------------------------------------------------------
!
! Compute 2D functions:
!
  do isurf=1,nlabel
    phi = 0.d0
    phiout = 2.d0*pi*qsaf(isurf)/ntheta
!
    ymet(1) = rbeg(isurf)
    ymet(2) = zaxis
    ymet(3) = 0.d0
    ymet(4) = 0.d0
    ymet(5) = 0.d0
!
    do j=1,ntheta
!
      call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_surf)
!
      rrr=ymet(1)
      zzz=ymet(2)
!
      call field_eq(rrr,phi,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
                   ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
      R_ts(j,isurf) = rrr
      Z_ts(j,isurf) = zzz
      bmod_ts(j,isurf) = sqrt(Br**2+Bp**2+Bz**2)
      sqgnorm_ts(j,isurf) = rrr/abs(Bp)
      Gfunc_ts(j,isurf) = ymet(5)/C_const(isurf)-phiout*dfloat(j)
    enddo
!
  enddo
!
  print *,'Field line integration: 2D functions done'
!
!-------------------------------------------------------------------------------
!
  end subroutine flint_for_Boozer
! -----------------------------------------------------------------
