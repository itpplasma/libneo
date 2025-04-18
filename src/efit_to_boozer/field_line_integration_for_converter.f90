!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module field_line_integration_for_coverter_mod
    logical :: prop = .true.
    double precision :: rmn,rmx,zmn,zmx
  end module

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_line_integration_for_coverter(rbeg,zaxis,ntheta,phimax,F_B,F_H,dF_B_dphi,ierr)
!
! Computes transformation functions $F_B$ and $F_H$, Eq.(58) and derivative $\rd F_B / \rd \varphi$, Eq.(60).
! All functions are normalized by end values - output is:
! F_B = $\hat F_B(\varphi) = F_B(\varphi) / F_B(\varphi_{max})$
! F_H = $\hat F_H(\varphi) = F_H(\varphi) / F_H(\varphi_{max})$
! dF_B_dphi = $\hat \rd F_B(\varphi)/ \rd \varphi =(\rd F_B(\varphi) / \rd \varphi) / F_B(\varphi_{max})$
! which corresponds to Eq.(62) in case $\vartheta_{B0} = \vartheta_{H0} = 0$.
!
  use field_line_integration_for_coverter_mod, only : prop,rmn,rmx,zmn,zmx
  use rhs_converter_mod, only : dz_dphi
  use field_eq_mod, only : nrad,nzet,rad,zet,icall_eq
!
  implicit none
!
  integer, parameter :: neq=4
  integer, parameter :: niter=10       !number of iterations for Newton method
  integer, parameter :: nstep_min=10   !minimum number of steps
!
  integer          :: ntheta,ierr,iter,j
  double precision :: rbeg,zaxis,phimax
!
  double precision :: phi,rrr,ppp,zzz
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: relerr,h,phiout,phi_sep,sig
!
  double precision, dimension(ntheta) :: F_B,F_H,dF_B_dphi
  double precision, dimension(neq)    :: ymet
!
  external :: rhs_converter
!
  relerr = 1.d-9
!
! Initialization of the field:
!
  if(prop) then
    prop = .false.
    icall_eq=-1
    rrr=1.d0
    ppp=0.d0 
    zzz=0.d0 
! 
    call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
              ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
! 
! Computation box:
    rmn=rad(1)
    rmx=rad(nrad)
    zmn=zet(1)
    zmx=zet(nzet)
  endif
! 
! End of initialization
!
! First field line integration to find ploidal period:
!
  ierr = 0
  h=6.28d0/nstep_min
!
  phi=0.d0
  phiout=h
  phi_sep=phi
  ymet(1)=rbeg
  ymet(2)=zaxis
  ymet(3)=0.d0
  ymet(4)=0.d0
  call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_converter)
  phi_sep=phi_sep+phiout
  sig=ymet(2)-zaxis
  do while(sig*(ymet(2)-zaxis).gt.0.d0)
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_converter)
    phi_sep=phi_sep+phiout
    if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
        ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
      print *,'Error: open field line'
      ierr = 1
      return
    endif
  enddo
  sig=ymet(2)-zaxis
  do while(sig*(ymet(2)-zaxis).gt.0.d0)
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_converter)
    phi_sep=phi_sep+phiout
    if( ymet(1).lt.rmn .or. ymet(1).gt.rmx .or.            &
        ymet(2).lt.zmn .or. ymet(2).gt.zmx ) then
      print *,'Error: open field line'
      ierr = 1
      return
    endif
  enddo
!
! Newton method
  do iter=1,niter
    phiout=(zaxis-ymet(2))/dz_dphi
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_converter)
    phi_sep=phi_sep+phiout
  enddo
!
  phimax = phi_sep
!
! Second field line integration to compute transformation funtions
!
  phi = 0.d0
  phiout = phimax/ntheta
!
  ymet(1)=rbeg
  ymet(2)=zaxis
  ymet(3)=0.d0
  ymet(4)=0.d0
!
  do j=1,ntheta
!
    call odeint_allroutines(ymet,neq,phi,phiout,relerr,rhs_converter)
!
    rrr=ymet(1)
    zzz=ymet(2)
!
    call field_eq(rrr,phi,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ  &
                 ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    F_H(j) = ymet(3)
    F_B(j) = ymet(4)
    dF_B_dphi(j) = (Br**2+Bp**2+Bz**2)*rrr/Bp
  enddo
!
  dF_B_dphi = dF_B_dphi/F_B(ntheta)
  F_H = F_H/F_H(ntheta)
  F_B = F_B/F_B(ntheta)
!
  end subroutine field_line_integration_for_coverter
!
