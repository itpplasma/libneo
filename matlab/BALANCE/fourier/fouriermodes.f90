!
  use rhs1_mod
  use rhs2_mod
  use bdivfree_mod, only : ntor
  use field_eq_mod, only : icall_eq,rtf,btf,nrad,nzet,rad,zet             &
                         , psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use field_c_mod,  only : icall_c
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: ndim,ndimc,nstep,nmap,ntotstep,niter,iter,i,j,ind,m,n,i1
  integer :: nsurf,isurf,nsurfmax,ndim_fc,ntheta,nsqpsi,nsubstep
  integer :: k1,k2,k3,k4,numbig,npoisep,nlabel,k
  double precision :: h,phi,phibeg,rbeg,zbeg,raxis,zaxis,hr,sig,eps_four
  double precision :: hh,hh0,sig0,rmn,rmx,zmn,zmx,psi_axis,dphidpsi
  double precision :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                      ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: sigma,themin,themax,theqtmin,theqtmax
  double precision :: rindmin,rindmax,sqpsimin,sqpsimax,x
  double precision :: theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz
  double precision :: rsmall_max,phitor_max,psisurf_max,hpsif,weight
  double precision :: flabel_min,flabel_max,hsqpsi,hsqpsim1
  double precision, dimension(4) :: xp,coef
  double precision, dimension(:), allocatable :: y,theta_qt,flabel
  double precision, dimension(:), allocatable :: rsmall,qsaf,psisurf,phitor
  double precision, dimension(:), allocatable :: sqpsi,startrind,phinorm_arr
  double precision, dimension(:), allocatable :: startrind_lab,sqpsi_lab
  double precision, dimension(:), allocatable :: qsaf_lab,volume
  double precision, dimension(:), allocatable :: rbigmax,rbigmin
  double precision, dimension(:,:), allocatable :: theta_of_theta_qt
  double precision, dimension(:,:), allocatable :: separ,dummy
  double complex, dimension(:),     allocatable :: cy
  double complex, dimension(:,:,:), allocatable :: armn,azmn,apsimn,athetmn
!
  external :: rhs1,rhs2
!
  open(1,file='fouriermodes.inp')
  read (1,*) mpol     !number of poloidal modes
  read (1,*) nstep    !number of integration steps
  read (1,*) nsqpsi   !grid size over radial variable for 1D quantities
  read (1,*) nlabel   !grid size over radial variable for 2D and 3D quantities
  read (1,*) ntheta   !grid size over poloidal angle
  read (1,*) nsurfmax !number of starting points between the 
                      !magnetic axis and right box boundary
                      !when searching for the separatrix
  close(1)
!
  numbig=1000000
  nmap=10         !number of maps for finding magnetic axis
  niter=10        !number of iterations for Newton method
  eps_four=0.1d0  !integration step in x for sin(x) 
!
! Initialization of the field:
!
  icall_eq=-1
  icall_c=-1
  rrr=0.d0
  ppp=0.d0
  zzz=0.d0
!
  call field(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! End of initialization
!
! Computation box:
  rmn=rad(1) 
  rmx=rad(nrad)
  zmn=zet(1)
  zmx=zet(nzet)
  open(1,file='box_size.dat')
  write(1,*) 'rmn = ',rmn
  write(1,*) 'rmx = ',rmx
  write(1,*) 'zmn = ',zmn
  write(1,*) 'zmx = ',zmx
  close(1)
!
  rbeg=0.5d0*(rmn+rmx)
  zbeg=0.5d0*(zmn+zmx)
!
!
  rrr=rbeg
  ppp=0.d0
  zzz=zbeg
!
  ntotstep=nstep*nmap
  h=2.d0*pi/nstep
!
! Search for the magnetic axis
!
  ndim=5
  isw_rhs1=1
!
  allocate(y(ndim))
  phi=0.d0
  y=0.d0
  y(1)=rbeg
  y(2)=zbeg
  do iter=1,niter
    phibeg=phi
    y(4:5)=0.d0
    do i=1,ntotstep
      call RK4D(y,ndim,phi,h,rhs1)
    enddo
    y(1:2)=y(4:5)/(phi-phibeg)
  enddo
  raxis=y(1)
  zaxis=y(2)
  open(1,file='axis.dat')
  write(1,*) 'raxis = ',raxis
  write(1,*) 'zaxis = ',zaxis
  close(1)
!
  call  field_eq(raxis,ppp,zaxis,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
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
!  deallocate(y)
!  ndim=4
  isw_rhs1=2
!  allocate(y(ndim))
!
  ndim_fc=1+4*ntor*(2*mpol+1)
  allocate(cy(ndim_fc))
!
!
  hr=(rmx-raxis)/nsurfmax
!
  rrr=raxis+hr
  zzz=zaxis
  call field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Direction of poloidal field:
  sigma=sign(1.d0,Bzet*Bphi)
!
  allocate(psisurf(nsurfmax))
  allocate(phinorm_arr(nsurfmax))
!
! Scan of flux surfaces
!
  allocate(dummy(2,numbig))
!
  surf: do isurf=1,nsurfmax
    i1=0
    phi=0.d0
    y(1)=raxis+hr*isurf
    y(2)=zaxis
    y(3)=0.d0
    y(4)=0.d0
    y(5)=0.d0
    call RK4D(y,ndim,phi,h,rhs1)
    sig=y(2)-zaxis
    i1=i1+1
    dummy(:,i1)=y(1:2)
    do while(sig*(y(2)-zaxis).gt.0.d0)
      call RK4D(y,ndim,phi,h,rhs1)
      i1=i1+1
      dummy(:,i1)=y(1:2)
      if( y(1).lt.rmn .or. y(1).gt.rmx .or.            &
          y(2).lt.zmn .or. y(2).gt.zmx ) then
        nsurf=isurf-1
        exit surf
      endif
    enddo
    sig=y(2)-zaxis
    do while(sig*(y(2)-zaxis).gt.0.d0)
      call RK4D(y,ndim,phi,h,rhs1)
      i1=i1+1
      dummy(:,i1)=y(1:2)
      if( y(1).lt.rmn .or. y(1).gt.rmx .or.            &
          y(2).lt.zmn .or. y(2).gt.zmx ) then
        nsurf=isurf-1
        exit surf
      endif
    enddo
!
! Newton method
    do iter=1,niter
      hh=(zaxis-y(2))/dz_dphi
      sig0=sig
      call RK4D(y,ndim,phi,hh,rhs1)
    enddo
!    i1=i1+1
    dummy(:,i1)=y(1:2)
!
    rrr=y(1)
    zzz=y(2)
!
    call  field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    psisurf(isurf)=psif-psi_axis
    rsmall_max=sqrt(abs(y(3))/pi)
    psisurf_max=psif-psi_axis
    phitor_max=y(4)/(2.d0*pi)
    phinorm_arr(isurf)=phitor_max
    if(allocated(separ)) deallocate(separ)
    npoisep=i1
    allocate(separ(2,npoisep))
    separ=dummy(:,1:npoisep)
  enddo surf
!
  nsurf=nsurf-1  !last point is bad, remove it
!print *,nsurf
  open(111,file="separ.dat")
  write(111,*) psisurf(nsurf)+psi_axis,psi_axis
  do i1=1,npoisep
    write(111,*) separ(:,i1)
  enddo
  close(111)
  deallocate(dummy,separ)
!
! Write the normalized toroidal flux on the equidistant grid 
! of the poloidal flux:
!
  phinorm_arr(1:nsurf)=phinorm_arr(1:nsurf)/phitor_max
  hpsif=(psisurf(nsurf)-psisurf(1))/(nsurf-1)
  open(375,file='phinorm_arr.dat')
  write (375,*) nsurf,psi_axis,hpsif,psisurf_max,-phitor_max
  i1=1
  do i=1,nsurf
    psif=psisurf(1)+hpsif*(i-1)
!    do while(psisurf(i1).le.psif)
    do while((psisurf(i1)-psif)/hpsif.le.0.d0)
      i1=i1+1
    enddo
    k3=min(nsurf-1,max(3,i1))
    k4=k3+1
    k2=k3-1
    k1=k3-2
    write (375,*)  phinorm_arr(k1)                              &
                 *(psif-psisurf(k2))/(psisurf(k1)-psisurf(k2))  &
                 *(psif-psisurf(k3))/(psisurf(k1)-psisurf(k3))  &
                 *(psif-psisurf(k4))/(psisurf(k1)-psisurf(k4))  &
                 + phinorm_arr(k2)                              &
                 *(psif-psisurf(k1))/(psisurf(k2)-psisurf(k1))  &
                 *(psif-psisurf(k3))/(psisurf(k2)-psisurf(k3))  &
                 *(psif-psisurf(k4))/(psisurf(k2)-psisurf(k4))  &
                 + phinorm_arr(k3)                              &
                 *(psif-psisurf(k1))/(psisurf(k3)-psisurf(k1))  &
                 *(psif-psisurf(k2))/(psisurf(k3)-psisurf(k2))  &
                 *(psif-psisurf(k4))/(psisurf(k3)-psisurf(k4))  &
                 + phinorm_arr(k4)                              &
                 *(psif-psisurf(k1))/(psisurf(k4)-psisurf(k1))  &
                 *(psif-psisurf(k2))/(psisurf(k4)-psisurf(k2))  &
                 *(psif-psisurf(k3))/(psisurf(k4)-psisurf(k3)) 
  enddo
  close(375)
!
!
!
!
! Re-interpolate for the equi-distant grid over sqrt(psi)
!
  allocate(sqpsi(0:nsurf),startrind(0:nsqpsi))
!
  sqpsi(0)=0.d0
  sqpsi(1:nsurf)=sqrt(abs(psisurf(1:nsurf)))
  rindmin=0
  rindmax=nsurf
!
  call invert_mono_reg(nsurf,sqpsi,rindmin,rindmax,nsqpsi, &
                       startrind,sqpsimin,sqpsimax)
  sqpsimin=sqpsimax/nsqpsi
!
  deallocate(psisurf)
  allocate(rsmall(nsqpsi),qsaf(nsqpsi),psisurf(nsqpsi),phitor(nsqpsi))
  allocate(volume(nsqpsi),rbigmin(nsqpsi),rbigmax(nsqpsi))
  rbigmin=raxis
  rbigmax=raxis
!
! Computation of iota, flux label, surface area
!
  do isurf=1,nsqpsi
print *,isurf,nsqpsi
    phi=0.d0
    y(1)=raxis+hr*startrind(isurf)
!write (111,*) startrind(isurf)
    y(2)=zaxis
    y(3)=0.d0
    y(4)=0.d0
    y(5)=0.d0
    rbigmin(isurf)=min(rbigmin(isurf),y(1))
    rbigmax(isurf)=max(rbigmax(isurf),y(1))
    call RK4D(y,ndim,phi,h,rhs1)
    sig=y(2)-zaxis
    do while(sig*(y(2)-zaxis).gt.0.d0)
      call RK4D(y,ndim,phi,h,rhs1)
      rbigmin(isurf)=min(rbigmin(isurf),y(1))
      rbigmax(isurf)=max(rbigmax(isurf),y(1))
    enddo
    sig=y(2)-zaxis
    do while(sig*(y(2)-zaxis).gt.0.d0)
      call RK4D(y,ndim,phi,h,rhs1)
      rbigmin(isurf)=min(rbigmin(isurf),y(1))
      rbigmax(isurf)=max(rbigmax(isurf),y(1))
    enddo
! Newton method
    do iter=1,niter
      hh=(zaxis-y(2))/dz_dphi
      sig0=sig
      call RK4D(y,ndim,phi,hh,rhs1)
      rbigmin(isurf)=min(rbigmin(isurf),y(1))
      rbigmax(isurf)=max(rbigmax(isurf),y(1))
    enddo
!
    aiota=2.d0*pi/phi
    rrr=y(1)
    zzz=y(2)
!
    call  field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    rsmall(isurf)=sqrt(abs(y(3))/pi)
    qsaf(isurf)=1.d0/aiota
    psisurf(isurf)=psif-psi_axis
    phitor(isurf)=y(4)/(2.d0*pi)
    volume(isurf)=abs(y(5))*pi
  enddo
!
  open(1,file='equil_r_q_psi.dat')
  write (1,*) '# ntor = ',ntor,' mpol = ',-mpol,':',mpol, 'nrad = ',nsqpsi
  write (1,*) '# rmax',rsmall_max,'   psimax',psisurf_max,'   phimax',phitor_max
  write (1,*) '# raidus r,  safety factor q, poloidal flux psi, toroidal flux phi,  d phi / d psi, &
                 geom. radius r, volume, R_beg, Z_beg, R_min, R_max'
  do i=1,nsqpsi
    i1=min(nsqpsi-1,max(2,i))
    dphidpsi=(phitor(i1+1)-phitor(i1-1))/(psisurf(i1+1)-psisurf(i1-1))
    write (1,*) sqrt(2.d0*abs(phitor(i)/btf)),qsaf(i),psisurf(i),phitor(i), &
                dphidpsi,rsmall(i),volume(i),raxis+hr*startrind(i),zaxis,   &
                rbigmin(i),rbigmax(i)
  enddo
  close(1)
!
!
! Re-interpolate sqrt(psi) to the equi-distant grid over r*q
!
  allocate(flabel(0:nsqpsi),sqpsi_lab(0:nlabel),startrind_lab(0:nlabel))
  allocate(qsaf_lab(0:nlabel))
!
  flabel(0)=0.d0
  flabel(1:nsqpsi)=rsmall(1:nsqpsi)*qsaf(1:nsqpsi)
  sqpsimin=0
!
  call invert_mono_reg(nsqpsi,flabel,sqpsimin,sqpsimax,nlabel, &
                       sqpsi_lab,flabel_min,flabel_max)
  flabel_min=flabel_max/nlabel
  sqpsimin=sqpsimax/nsqpsi
!
! sqpsi_lab - sqrt(psi) on the equidistant grid of flabel
! 
! Re-calculate startpoints for the equidistant grid of flabel
!
  hsqpsi=sqpsimax/nsqpsi
  hsqpsim1=1.d0/hsqpsi
  startrind_lab(0)=0.d0
  do i=1,nlabel
    k=int(sqpsi_lab(i)*hsqpsim1)
    k=min(max(0,k-2),nsqpsi-4)
    do j=1,4
      xp(j)=(k+j)*hsqpsi
    enddo
    call coefs_bdf(sqpsi_lab(i),xp,hsqpsim1,coef)
    startrind_lab(i)=sum(coef*startrind(k+1:k+4))
    qsaf_lab(i)=flabel_min*i/sum(coef*rsmall(k+1:k+4))
  enddo
!
  allocate(theta_qt(0:nstep))
  allocate(theta_of_theta_qt(nlabel,0:ntheta))
  themin=0.d0
  themax=2.d0*pi
!
  allocate(armn(-mpol:mpol,ntor,nlabel))
  allocate(azmn(-mpol:mpol,ntor,nlabel))
  allocate(apsimn(-mpol:mpol,ntor,nlabel))
  allocate(athetmn(-mpol:mpol,ntor,nlabel))
!
! Fourier expansion of the perturbation field over the poloidal angle
!
  do isurf=1,nlabel
!
    aiota=1.d0/qsaf_lab(isurf)
    h=2.d0*pi/aiota/nstep
    hh=eps_four/max(1,mpol)
    nsubstep=ceiling(h/hh)
    hh=h/nsubstep
    phi=0.d0
    cy=dcmplx(0.d0,0.d0)
    cy(1)=dcmplx(raxis+startrind_lab(isurf)*hr,zaxis)
    theta_qt(0)=0.d0
!
    do i=1,nstep
      do j=1,nsubstep
        call RK4C(cy,ndim_fc,phi,hh,rhs2)
      enddo
      rrr=dble(cy(1))
      zzz=aimag(cy(1))
      ppp=phi
      theta_qt(i)=mod(sigma*atan2(zzz-zaxis,rrr-raxis)+2*pi,2*pi)
    enddo
    if(theta_qt(nstep).lt.pi) theta_qt(nstep)=theta_qt(nstep)+2.d0*pi
    theta_qt=theta_qt*2.d0*pi/theta_qt(nstep)
!
    call invert_mono_per(nstep,theta_qt,themin,themax,ntheta,         &
                         theta_of_theta_qt(isurf,:),theqtmin,theqtmax)
!
    do i=0,ntheta
      theta_of_theta_qt(isurf,i)=theta_of_theta_qt(isurf,i)-themax*i/ntheta
    enddo
!
    ind=1
    do n=1,ntor
      do m=-mpol,mpol
        ind=ind+1
        armn(m,n,isurf)=cy(ind)
        ind=ind+1
        azmn(m,n,isurf)=cy(ind)
        ind=ind+1
        apsimn(m,n,isurf)=cy(ind)
        ind=ind+1
        athetmn(m,n,isurf)=cy(ind)
      enddo
    enddo
print *,isurf,nlabel,', ',nstep*nsubstep,' integration steps'
!
  enddo
!
  open(1,form='unformatted',file='theta_of_theta_qt_flabel.dat')
  write (1) nsqpsi,nlabel,ntheta,sqpsimin,sqpsimax,flabel_min,flabel_max &
           ,raxis,zaxis,psi_axis,sigma
  write (1) theta_of_theta_qt
  write (1) flabel
  close(1)
!
  open(11,file='thetabooz.dat')
  do isurf=1,nlabel
    write (11,*) theta_of_theta_qt(isurf,:)
  enddo
  close(11)
!
  open(1,form='unformatted',file='amn.dat')
  write (1) ntor,mpol,nlabel,flabel_min,flabel_max
  write (1) apsimn,athetmn
  close(1)
!
  stop
!
  do n=1,3
    do m=-mpol,mpol
      do i=1,nlabel
        write (100*n+m,*) rsmall(i),real(apsimn(m,n,i)),real(aimag(apsimn(m,n,i)))
        write (1000+100*n+m,*) rsmall(i),real(athetmn(m,n,i)),real(aimag(athetmn(m,n,i)))
        write (2000+100*n+m,*) rsmall(i),real(armn(m,n,i)),real(aimag(armn(m,n,i)))
        write (3000+100*n+m,*) rsmall(i),real(azmn(m,n,i)),real(aimag(azmn(m,n,i)))
      enddo
    enddo
  enddo
!
!
  stop
  end
!
