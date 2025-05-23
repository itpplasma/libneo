module field_sub

implicit none

integer, parameter :: dp = kind(1.0d0)

real(dp) :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2

! Make temporary variables threadprivate
!$omp threadprivate(psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2)

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine read_field_input(input_file)
  use input_files, only : iunit, gfile, pfile, convexfile, fluxdatapath, ieqfile
  use field_c_mod, only : ntor, icftype
  use field_eq_mod, only : nwindow_r, nwindow_z
  use field_mod, only : ipert, iequil, ampl
  use inthecore_mod, only : cutoff

  character(*), intent(in), optional :: input_file

  if (present(input_file)) then
    open(iunit, file=input_file, status='old', action='read')
  else
    open(iunit, file='field_divB0.inp', status='old', action='read')
  end if

  read(iunit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
                             ! 3=plas+vac with derivatives
  read(iunit,*) iequil       ! 0=perturbation alone, 1=with equilibrium
  read(iunit,*) ampl         ! amplitude of perturbation, a.u.
  read(iunit,*) ntor         ! number of toroidal harmonics
  read(iunit,*) cutoff       ! inner cutoff in psi/psi_a units
  read(iunit,*) icftype      ! type of coil file
  read(iunit,*) gfile        ! equilibrium file
  read(iunit,*) pfile        ! coil        file
  read(iunit,*) convexfile   ! convex file for stretchcoords
  read(iunit,*) fluxdatapath ! directory with data in flux coord.
  read(iunit,*) nwindow_r    ! widow size for filtering of psi array over R
  read(iunit,*) nwindow_z    ! widow size for filtering of psi array over Z
  read(iunit,*,err=1) ieqfile! equilibrium file type (0 - old, 1 - EFIT)
1 close(iunit)
end subroutine read_field_input

subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

  use field_c_mod, only : icall_c
  use field_mod, only : icall, ipert, iequil, ampl
  use inthecore_mod, only : incore
  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(in) :: r, z
  real(dp) :: p,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  real(dp) :: rm,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc &
                     ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc

  if(icall .eq. 0) then
    icall = 1
    call read_field_input
    print *, 'Perturbation field',ipert,'Ampl',ampl
    if(icall_c.eq.-1) ipert=1
  end if

  call stretch_coords(r,z,rm,zm)

  if(iequil.eq.0) then
    call set_zero(Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,   &
                  dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  else
    call field_eq(rm,p,zm,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  end if

  if(ipert.gt.0) then

    if(ipert.gt.1) then
      call inthecore(rm,zm)
    else
      incore=-1
    endif

    ! vacuum perturbation coil field:

    call field_c(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc   &
                ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)

    call add_scaled(Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ  &
      ,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc,ampl)

    if(incore.gt.-1) then
      ! perturbation coil field with plasma shielding:

      if(ipert.eq.2) then

        call field_fourier(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc           &
                          ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)

        Br = Br + Brc*ampl
        Bp = Bp + Bpc*ampl
        Bz = Bz + Bzc*ampl

      else if(ipert.eq.3) then

        call field_fourier_derivs(rm,p,zm,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc    &
                                 ,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc)

        call add_scaled( &
          Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ  &
         ,Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,dBpdRc,dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZ &
         ,ampl)

      end if

    end if

   end if

end subroutine field

! ========================================================================
subroutine field_eq(r,ppp,z,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                   ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

  use input_files, only : ieqfile
  use field_eq_mod, only : use_fpol,skip_read,icall_eq,nrad,nzet,icp,nwindow_r,&
    nwindow_z,psib,btf,rtf,hrad,hzet,psi_axis,psi_sep,&
    psi,psi0,splfpol,splpsi,rad,zet,imi,ima,jmi,jma,ipoint
  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(in) :: r, ppp, z
  real(dp), intent(out) :: Brad, Bphi, Bzet, dBrdR, dBrdp, dBrdZ
  real(dp), intent(out) :: dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

  integer :: ierr,i

  real(dp) :: rrr,zzz
  real(dp) :: psihat,fpol,fpol_prime

  associate(dummy => ppp)
  end associate

  !-------first call: read data from disk-------------------------------
  if (icall_eq < 1) then
    if (.not. skip_read) then
      if (ieqfile == 0) then
        call read_dimeq0(nrad, nzet)
      elseif (ieqfile == 2) then
        call read_dimeq_west(nrad, nzet)
      else
        call read_dimeq1(nrad, nzet)
      end if
      allocate(rad(nrad), zet(nzet))
      allocate(psi0(nrad, nzet), psi(nrad, nzet))
    end if
    if (use_fpol) then
      if (.not. skip_read) then
        allocate(splfpol(0:5, nrad))
        call read_eqfile2(nrad, nzet, psi_axis, psi_sep, btf, rtf, &
                          splfpol(0, :), rad, zet, psi)
      end if
      psib = -psi_axis
      psi_sep = (psi_sep - psi_axis) * 1.d8
      splfpol(0, :) = splfpol(0, :) * 1.d6
      call spline_fpol
    else
      if (.not. skip_read) then
        if (ieqfile == 0) then
          call read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
        elseif (ieqfile == 2) then
          call read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)
        else
          call read_eqfile1(nrad, nzet, psib, btf, rtf, rad, zet, psi)
        end if
      end if
    end if

    ! Filtering:
    do i=1,nzet
      call window_filter(nrad,nwindow_r,psi(:,i),psi0(:,i))
    end do

    do i=1,nrad
      call window_filter(nzet,nwindow_z,psi0(i,:),psi(i,:))
    end do
    ! End filtering

    rad = rad*1.d2 ! cm
    zet = zet*1.d2 ! cm
    rtf = rtf*1.d2 ! cm
    psi = psi*1.d8
    psib= psib*1.d8
    btf = btf*1.d4

    psi=psi+psib

    hrad = rad(2) - rad(1)
    hzet = zet(2) - zet(1)

    ! rectangular domain:
    allocate( imi(nzet),ima(nzet),jmi(nrad),jma(nrad) )
    imi = 1
    ima = nrad
    jmi = 1
    jma = nzet

    !  Computation of the number of data in splpsi
    icp = 0
    do i=1,nzet
      if ( imi(i) .gt. 0 .and. ima(i) .gt. 0 ) then
         icp = icp + ima(i) - imi(i) + 1
      end if
    end do
    write(6,*) 'number of points in the table:  ',icp

    allocate( splpsi(6,6,icp), ipoint(nrad,nzet) )

    call s2dcut(nrad,nzet,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)

    if(icall_eq.eq.-1) then
      ! Quit after initialization with zero field
      call set_zero(Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ, &
                    dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
      icall_eq = 1
      return
    end if
    icall_eq = 1
  end if
  ! ------- end first call ----------------------------------------------

  rrr=max(rad(1),min(rad(nrad),r))
  zzz=max(zet(1),min(zet(nzet),z))

  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
                psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

  Brad = -dpsidz/rrr
  Bzet =  dpsidr/rrr

  ! axisymmetric case:
  dBrdp = 0.
  dBpdp = 0.
  dBzdp = 0.

  dBrdR = -d2psidrdz/rrr+dpsidz/rrr**2
  dBzdZ =  d2psidrdz/rrr
  dBrdZ = -d2psidz2/rrr
  dBzdR =  d2psidr2/rrr-dpsidr/rrr**2

  if(use_fpol) then
    psihat=psif/psi_sep
    if(psihat.gt.1.d0) then
      fpol=splfpol(0,nrad)
      fpol_prime=0.d0
    else
      call splint_fpol(psihat,fpol,fpol_prime)
    end if
    Bphi = fpol/rrr
    dBpdR = fpol_prime*dpsidr/(psi_sep*rrr)-fpol/rrr**2
    dBpdZ = fpol_prime*dpsidz/(psi_sep*rrr)
  else
    Bphi = btf*rtf/rrr
    dBpdR = -btf*rtf/rrr**2
    dBpdZ = 0.
  end if

end subroutine field_eq

! ----------- Runov Original Method --------------------------------
subroutine read_dimeq0(nrad,nzet)
  use input_files, only : eqfile

  implicit none

  integer, intent(out) :: nrad, nzet

  open(11,file=eqfile)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)

  read(11,111) nrad
  read(11,111) nzet
111  format(12x,i3)

  close(11)
  return
end subroutine read_dimeq0

subroutine read_eqfile0(nrad, nzet, psib, btf, rtf, rad, zet, psi)
  use input_files, only : eqfile
  use libneo_kinds, only : dp

  implicit none

  integer :: i,j,k

  integer, intent(in) :: nrad, nzet
  real(dp), intent(out) :: psib, btf, rtf
  real(dp), intent(out) :: rad(nrad), zet(nzet)
  real(dp), intent(out) :: psi(nrad,nzet)

  integer :: dum

  open(11,file=eqfile)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)

  read(11,111) dum !nrad
  read(11,111) dum !nzet
  read(11,112) psib
  read(11,112) btf
  read(11,112) rtf
  read(11,*)
  read(11,*)

  print *, nrad, nzet, psib, btf, rtf

  read(11,113)(rad(i),i=1,nrad)
  read(11,*)
  read(11,*)
  read(11,113)(zet(i),i=1,nzet)
  read(11,*)
  read(11,*)
  read(11,113)((psi(j,k),j=1,nrad),k=1,nzet)

  close(11)
  return

111  format(12x,i3)
112  format(12x,f21.2)
113  format(5(e17.4))
end subroutine read_eqfile0


! ----------- Read gfile directly --------------------------------
subroutine read_dimeq1(nwEQD,nhEQD)
  use input_files, only : iunit,gfile
  implicit none

  integer, intent(out) :: nwEQD, nhEQD

  integer :: i, idum
  character(len=10) :: dummy(6)

  open(unit=iunit,file=trim(gfile),status='old',action='read')
  read(iunit,2000)(dummy(i),i=1,6),idum,nwEQD,nhEQD
  close(iunit)
  return

2000  format(6a8,3i4)
end subroutine read_dimeq1


subroutine read_eqfile1(nwEQD,nhEQD,psiSep, bt0, rzero, rad, zet, psiRZ)
  use input_files, only : iunit, gfile
  use libneo_kinds, only : dp

  implicit none

  integer, intent(inout) :: nwEQD, nhEQD
  real(dp), intent(out) :: bt0, rzero, psiSep
  real(dp), intent(out) :: rad(nwEQD), zet(nhEQD)
  real(dp), dimension(nwEQD,nhEQD), intent(out) :: psiRZ

  integer :: gunit, idum
  character(len=10) :: dummy(6)
  integer :: i,j
  real(dp) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  real(dp) :: plas_cur, psiAxis
  real(dp), dimension(nwEQD) :: fpol,pres,ffprim,pprime,qpsi

  integer :: n_bndyxy,nlimEQD
  real(dp), dimension(:), allocatable :: LCFS, limEQD

  gunit=iunit

  open(unit=gunit,file=trim(gfile),status='old',action='read')

  ! Equilibrium Parameters
  read(gunit,2000)(dummy(i),i=1,6),idum,nwEQD,nhEQD
  write(*,*) 'READ_EQFILE1: ',trim(gfile),nwEQD,nhEQD
  read(gunit,2010,end=55,err=250)xdim,zdim,rzero,r1,zmid
  write(*,*) xdim, zdim, rzero, r1, zmid
  read(gunit,2010,end=55,err=250)rmaxis,zmaxis,psiAxis,psiSep,bt0
  write(*,*) rmaxis,zmaxis,psiAxis,psiSep,bt0
  read(gunit,2010,end=55,err=250)plas_cur,psiAxis,xdum,rmaxis,xdum
  write(*,*) plas_cur,psiAxis,xdum,rmaxis,xdum
  read(gunit,2010,end=55,err=250)zmaxis,xdum,psiSep,xdum,xdum
  write(*,*) zmaxis,xdum,psiSep,xdum,xdum
  read(gunit,2010,end=55,err=250)(fpol(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(pres(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(ffprim(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(pprime(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
  read(gunit,2010,end=55,err=250)(qpsi(i),i=1,nwEQD)
  print *, 'Equilibrium Done.', trim(gfile)
  ! Boundary Data
  read(gunit,*,end=55,err=250)n_bndyxy,nlimEQD

  if (n_bndyxy > 0) then
    allocate(LCFS(2*n_bndyxy))
    read(gunit,2010,end=55,err=250)(LCFS(i),i=1,2*n_bndyxy)
  end if

  if (nlimEQD > 0) then
    allocate(limEQD(2*nlimEQD))
    read(gunit,2010,end=55,err=250)(limEQD(i),i=1,2*nlimEQD)
  end if

  close(gunit)

  call set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  return

2000  format(6a8,3i4)
2010  format(5(e16.9))
55    print *, 'READ_EQFILE1: Early EOF in',trim(gfile); STOP
250   print *, 'READ_EQFILE1: Error reading ',trim(gfile); STOP

end subroutine read_eqfile1


subroutine read_eqfile2(nwEQD,nhEQD,psiAxis,psiSep,bt0,rzero,fpol,rad,zet,psiRZ)
  use input_files, only : iunit, gfile
  use libneo_kinds, only : dp
  implicit none

  integer, intent(inout) :: nwEQD, nhEQD
  real(dp), intent(out) :: bt0, rzero, psiAxis, psiSep
  real(dp), dimension(nwEQD), intent(out) :: fpol
  real(dp), intent(out) :: rad(nwEQD), zet(nhEQD)
  real(dp), dimension(nwEQD,nhEQD), intent(out) :: psiRZ

  integer :: gunit, idum
  character(len=10) :: dummy(6)
  integer :: i,j
  real(dp) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  real(dp) :: plas_cur
  real(dp), dimension(nwEQD) :: pres,ffprim,pprime,qpsi

  integer :: n_bndyxy,nlimEQD
  real(dp), dimension(:), allocatable :: LCFS, limEQD

  gunit=iunit

  open(unit=gunit,file=trim(gfile),status='old',action='read')

  ! Equilibrium Parameters
  read(gunit,2000)(dummy(i),i=1,6),idum,nwEQD,nhEQD
  write(*,*) 'READ_EQFILE1: ',trim(gfile),nwEQD,nhEQD
  read(gunit,2010,end=55,err=250)xdim,zdim,rzero,r1,zmid
  write(*,*) xdim, zdim, rzero, r1, zmid
  read(gunit,2010,end=55,err=250)rmaxis,zmaxis,psiAxis,psiSep,bt0
  write(*,*) rmaxis,zmaxis,psiAxis,psiSep,bt0
  read(gunit,2010,end=55,err=250)plas_cur,psiAxis,xdum,rmaxis,xdum
  write(*,*) plas_cur,psiAxis,xdum,rmaxis,xdum
  read(gunit,2010,end=55,err=250)zmaxis,xdum,psiSep,xdum,xdum
  write(*,*) zmaxis,xdum,psiSep,xdum,xdum
  read(gunit,2010,end=55,err=250)(fpol(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(pres(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(ffprim(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)(pprime(i),i=1,nwEQD)
  read(gunit,2010,end=55,err=250)((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
  read(gunit,2010,end=55,err=250)(qpsi(i),i=1,nwEQD)
  print *, 'Equilibrium Done.', trim(gfile)
  ! Boundary Data
  read(gunit,*,end=55,err=250)n_bndyxy,nlimEQD
  allocate(LCFS(2*n_bndyxy))
  allocate(limEQD(2*nlimEQD))
  read(gunit,2010,end=55,err=250)(LCFS(i),i=1,2*n_bndyxy)
  read(gunit,2010,end=55,err=250)(limEQD(i),i=1,2*nlimEQD)

  close(gunit)

  call set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  return

2000  format(6a8,3i4)
2010  format(5(e16.9))
55    print *, 'READ_EQFILE2: Early EOF in',trim(gfile); STOP
250   print *, 'READ_EQFILE2: Error reading ',trim(gfile); STOP

end subroutine read_eqfile2


subroutine set_eqcoords(nwEQD,nhEQD,xdim,zdim,r1,zmid,rad,zet)
  use libneo_kinds, only : dp

  implicit none

  integer, intent(in) :: nwEQD, nhEQD
  real(dp), intent(in) :: xdim, zdim, r1, zmid
  real(dp), intent(out) :: rad(nwEQD), zet(nhEQD)

  integer :: j,k
  real(dp) :: z1

  do j=1,nwEQD
    rad(j) = r1 + (j-1)*(xdim/(nwEQD-1))
  end do

  z1 = zmid - zdim/2.0 ! check this definition wrt zmid
  do k=1,nhEQD ! runov chooses lower, probe chooses upper
    zet(k) = z1 + (k-1)*(zdim/(nhEQD-1))
  end do

end subroutine set_eqcoords


! ===========================================================================
! Input of axisymmetric equilibrium for WEST tokamak
subroutine read_dimeq_west(nrad, nzet)

  use input_files, only : iunit, gfile

  implicit none

  integer, intent(out) :: nrad, nzet

  open(unit = iunit, file = trim(gfile), status = 'old', action = 'read')
  read(iunit, *) nrad, nzet
  close(iunit)

end subroutine read_dimeq_west


!----------------------------------------------------------------------
subroutine read_eqfile_west(nrad, nzet, psib, btf, rtf, rad, zet, psi)

  use input_files, only : iunit, gfile
  use libneo_kinds, only : dp

  implicit none

  integer, intent(inout) :: nrad, nzet
  real(dp), intent(out) :: psib, btf, rtf
  real(dp), intent(out) :: rad(nrad), zet(nzet)
  real(dp), intent(out) :: psi(nrad, nzet)
  integer :: ir

  psib = 0.0d0

  open(unit = iunit, file = trim(gfile), status = 'old', action = 'read')
  read (iunit, *) nrad, nzet
  read (iunit, *) btf
  read (iunit, *) rad
  read (iunit, *) zet
  do ir = 1, nrad
     read (iunit, *) psi(ir, :)
  end do
  close(iunit)
  rtf = 0.5d0 * (rad(1) + rad(nrad))
  btf = btf / rtf

end subroutine read_eqfile_west


! ===========================================================================
subroutine field_c(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

  use field_c_mod, only : icall_c,ntor,nr,np,nz,icftype,rmin,pmin,zmin,rmax,pmax,zmax
  use libneo_kinds, only : dp

  implicit none

  real(dp) :: rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ  &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  real(dp), dimension(:,:,:), allocatable :: Br,Bp,Bz

  !-------first call: read data from disk-------------------------------
  if(icall_c .lt. 1) then
    print *,'coils: file type = ',icftype
    if(icftype.eq.1) then
      nr=129
      np=37
      nz=129
    elseif(icftype.eq.2) then
      nr=129
      np=33
      nz=129
    elseif(icftype.eq.3) then
      nr=129
      np=37
      nz=131
      icftype=1
    elseif(icftype.eq.4) then
      call read_sizes(nr,np,nz)
    else
      print *,'field_c: unknown coil file type'
      stop
    end if
    allocate(Br(nr,np,nz),Bp(nr,np,nz),Bz(nr,np,nz))

    if(icftype.lt.4) then
      call read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    else
      call read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)
    end if

    print *,'coils: nr,np,nz = ',nr,np,nz
    print *,'coils: rmin,rmax = ',rmin,rmax
    print *,'coils: zmin,zmax = ',zmin,zmax
    print *,'coils: pmin,pmax = ',pmin,pmax

    call vector_potentials(nr,np,nz,ntor,rmin,rmax,pmin,pmax,zmin,zmax,br,bp,bz)

    deallocate(Br,Bp,Bz)

    if(icall_c.eq.-1) then
      ! Quit after initialization with zero field
      call set_zero(Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ, &
                    dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
      icall_c = 1
      return
    end if
    icall_c = 1
  end if
  !------- end first call ----------------------------------------------

  call field_divfree(rrr,ppp,zzz,Brad,Bphi,Bzet,dBrdR,dBrdp,dBrdZ    &
                    ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

end subroutine field_c


! ===========================================================================
subroutine read_field0(rad,phi,zet,rmin,pmin,zmin,hrm1,hpm1,hzm1,Br,Bp,Bz)

  use input_files, only : cfile
  use math_constants, only : pi

  implicit real(8) (a-h, o-z)

  integer, parameter :: nr=64, np=37, nz=64
  integer :: i, j, k, icall

  dimension Bz(nr,np,nz)
  dimension Br(nr,np,nz),Bp(nr,np,nz)
  dimension rad(nr), phi(np), zet(nz)
  data icall/0/
  save

  associate(dummy => hrm1)
  end associate
  associate(dummy => hpm1)
  end associate
  associate(dummy => hzm1)
  end associate

  !-------first call: read data from disk-------------------------------
  open(1,file=cfile,status='old',action='read')
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)

  !---Input B      -->T = V*s/m/m
  do j=1,np-1  !only npmax-1 points are given
    do k=nz,1,-1  !reverse order of probe data
      do i=1,nr
        read(1,*) Br(i,j,k), Bp(i,j,k), Bz(i,j,k)

        Br(i,j,k) = Br(i,j,k)*1.d4
        Bp(i,j,k) = Bp(i,j,k)*1.d4
        Bz(i,j,k) = Bz(i,j,k)*1.d4

      end do
      read(1,*)
    end do
    read(1,*)
  end do
  close(1)

  rmin = 84.d0
  rmax = 254.d0
  zmin = -160.d0
  zmax = 160.d0
  pmin = 0.d0
  pmax = 2.d0*pi

  hrad = (rmax - rmin)/(nr-1)
  hphi = (pmax - pmin)/(np-1)
  hzet = (zmax - zmin)/(nz-1)

  do i=1,nr
    rad(i) = rmin + hrad*(i-1)
  end do
  do i=1,np
    phi(i) = pmin + hphi*(i-1)
  end do
  do i=1,nz
    zet(i) = zmin + hzet*(i-1)
  end do

  do i=1,nr
    do k=1,nz
      Br(i,np,k) = Br(i,1,k)
      Bp(i,np,k) = Bp(i,1,k)
      Bz(i,np,k) = Bz(i,1,k)
    end do
  end do
end subroutine read_field0


subroutine read_field1(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)

  use input_files, only : iunit,pfile
  use libneo_kinds, only : dp
  use math_constants, only : pi

  implicit none

  integer, intent(in) :: icftype, nr, np, nz
  real(dp), intent(out) :: rmin, rmax, pmin, pmax, zmin, zmax
  real(dp), dimension(nr,np,nz), intent(out) :: Br, Bp, Bz

  integer :: i,j,k
  real(dp) :: xdim,zdim,zmid,dum

  open(iunit,file=trim(pfile),status='old',action='read')
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)  !PROBE
  read(iunit,*)
  read(iunit,*)
  if(icftype.eq.2) then
    read(iunit,*)    !New Format
    read(iunit,*)    !New Format
  endif

  !---Input B      -->T = V*s/m/m
  do j=1,np-1   !only npmax-1 points are given
    do k=nz,1,-1  !reverse order of probe data
      do i=1,nr
        if(icftype.eq.1) then
          ! Old Format
          read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
        elseif(icftype.eq.2) then
          ! New Format
          read(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
        end if

        !Convert to CGS
        Br(i,j,k) = Br(i,j,k)*1.d4
        Bp(i,j,k) = Bp(i,j,k)*1.d4
        Bz(i,j,k) = Bz(i,j,k)*1.d4
      end do
      read(iunit,*)
    end do
    read(iunit,*)
  end do
  close(iunit)

  xdim=170.d0
  rmin=84.d0
  rmax=rmin+xdim

  pmin = 0.
  pmax = 2.*pi

  zdim=320.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0

  do i=1,nr
     do k=1,nz
        Br(i,np,k) = Br(i,1,k)
        Bp(i,np,k) = Bp(i,1,k)
        Bz(i,np,k) = Bz(i,1,k)
     end do
  end do
end subroutine read_field1


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine stretch_coords(r,z,rm,zm)
  use input_files, only : iunit,convexfile
  use libneo_kinds, only : dp
  use math_constants, only : TWOPI

  implicit none

  real(dp), intent(in) :: r, z
  real(dp), intent(out) :: rm, zm

  integer icall, i, j, nrz ! number of points "convex wall" in input file
  integer, parameter :: nrzmx=100 ! possible max. of nrz
  integer, parameter :: nrhotht=360
  integer :: iflag
  real(dp) R0,Rw, Zw, htht, a, b, rho, tht, rho_c, delta, dummy
  real(dp), dimension(0:1000):: rad_w, zet_w ! points "convex wall"
  real(dp), dimension(:), allocatable :: rho_w, tht_w
  real(dp), dimension(nrhotht) :: rho_wall, tht_wall ! polar coords of CW
  data icall /0/, delta/1./
  save
  !----------- 1st call --------------------------------------------------------
  !$omp critical
  if(icall .eq. 0) then
    icall = 1
    nrz = 0
    rad_w = 0.
    zet_w = 0.
    open(iunit, file=trim(convexfile), status='old', action='read')
    do i=1,nrzmx
      read(iunit,*,END=10)rad_w(i),zet_w(i)
      nrz = nrz + 1
    end do
10  continue
    close(iunit)

    allocate(rho_w(0:nrz+1), tht_w(0:nrz+1))
    R0 = (maxval(rad_w(1:nrz)) +  minval(rad_w(1:nrz)))*0.5
    do i=1,nrz
      rho_w(i) = sqrt( (rad_w(i)-R0)**2 + zet_w(i)**2 )
      tht_w(i) = atan2(zet_w(i),(rad_w(i)-R0))
      if(tht_w(i) .lt. 0.) tht_w(i) = tht_w(i) + TWOPI
    end do

    ! make sure points are ordered according to tht_w.
    do
      iflag = 0
      do i=1,nrz-1
        if (tht_w(i) > tht_w(i+1)) then
          iflag = 1
          dummy = rad_w(i+1)
          rad_w(i+1) = rad_w(i)
          rad_w(i) = dummy
          dummy = zet_w(i+1)
          zet_w(i+1) = zet_w(i)
          zet_w(i) = dummy
          dummy = rho_w(i+1)
          rho_w(i+1) = rho_w(i)
          rho_w(i) = dummy
          dummy = tht_w(i+1)
          tht_w(i+1) = tht_w(i)
          tht_w(i) = dummy
        end if
      end do
      if (iflag == 0) exit
    end do
    rad_w(0) = rad_w(nrz)
    zet_w(0) = zet_w(nrz)
    tht_w(0) = tht_w(nrz) - TWOPI
    rho_w(0) = rho_w(nrz)
    rad_w(nrz+1) = rad_w(1)
    zet_w(nrz+1) = zet_w(1)
    tht_w(nrz+1) = tht_w(1) + TWOPI
    rho_w(nrz+1) = rho_w(1)

    htht = TWOPI/(nrhotht-1)
    do i=2,nrhotht
      tht_wall(i) = htht*(i-1)
      do j=0,nrz
        if(tht_wall(i).ge.tht_w(j) .and. tht_wall(i).le.tht_w(j+1)) then
          if( abs((rad_w(j+1) - rad_w(j))/rad_w(j)) .gt. 1.e-3) then
            a = (zet_w(j+1) - zet_w(j))/(rad_w(j+1) - rad_w(j))
            b = zet_w(j) - a*(rad_w(j) - R0)
            Rw = b/(tan(tht_wall(i)) - a) + R0
            Zw = a*(Rw - R0) + b
          else
            a = (rad_w(j+1) - rad_w(j))/(zet_w(j+1) - zet_w(j))
            b = rad_w(j)-R0 - a*zet_w(j)
            Zw = b/(1./tan(tht_wall(i)) - a)
            Rw = a*Zw + b + R0
          end if
        end if
      end do
      rho_wall(i) = sqrt((Rw-R0)**2 + Zw**2)
    end do
    tht_wall(1) = 0.
    rho_wall(1) = rho_wall(nrhotht)

  end if
  !$omp end critical
  !----------- end of the 1st call --------------------------------------------
  rm = r
  zm = z
  rho = sqrt((r-R0)**2 + z**2)
  tht = atan2(z,(r-R0))
  if(tht .lt. 0.) tht = tht + TWOPI
  i = modulo(int(tht/htht), nrhotht-1) + 1
  rho_c = (rho_wall(i+1) - rho_wall(i))/(tht_wall(i+1) - tht_wall(i))   &
       *(tht - tht_wall(i)) + rho_wall(i)

  if(rho .ge. rho_c) then
     rho = rho_c + delta*atan2((rho-rho_c), delta)
     rm = rho*cos(tht) + R0
     zm = rho*sin(tht)
  end if

end subroutine stretch_coords


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine inthecore(R,Z)

  use inthecore_mod, only : prop,npoi,ijumpb,ibeg,iend,epssep,rc,zc,twopi,sig,&
    psi_sep,psi_cut,sigpsi,cutoff,rho2i,theti,incore,vacf,dvacdr,dvacdz,&
    d2vacdr2,d2vacdrdz,d2vacdz2,plaf,dpladr,dpladz,d2pladr2,d2pladrdz,d2pladz2
  use input_files,  only : iunit,fluxdatapath
  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(in) :: R,Z

  integer :: i
  real(dp) :: rho2,thet,xx,yy
  real(dp) :: weight,dweight,ddweight
  real(dp), dimension(4) :: x,y
  real(dp), dimension(:), allocatable :: ri,zi

  if(prop) then
    prop=.false.
    open(iunit,file=trim(fluxdatapath)//'/separ.dat')
    read(iunit,*) x(1:2)
    psi_sep=x(2)+(x(1)-x(2))*(1.d0-epssep)
    psi_cut=x(2)+(x(1)-x(2))*cutoff
    sigpsi=sign(1.d0,psi_sep-psi_cut)
    npoi=0
    do while(npoi.ge.0)
      npoi=npoi+1
      read(iunit,*,end=1) rc
    end do
 1  allocate(ri(0:npoi),zi(0:npoi),rho2i(0:npoi),theti(0:npoi))
    ri=0.d0
    zi=0.d0
    close(iunit)
    open(iunit,file=trim(fluxdatapath)//'/separ.dat')
    read(iunit,*)
    do i=1,npoi-1
      read(iunit,*) ri(i),zi(i)
    end do
    close(iunit)
    rc=sum(ri(1:npoi-1))/(npoi-1)
    zc=sum(zi(1:npoi-1))/(npoi-1)
    rho2i=(ri-rc)**2+(zi-zc)**2
    theti=atan2(zi-zc,ri-rc)
    sig=theti(2)-theti(1)
    do i=2,npoi-2
      if((theti(i+1)-theti(i))*sig.lt.0.d0) then
        ijumpb=i
        exit
      end if
    end do
    twopi=8.d0*atan2(1.d0,1.d0)
    ri(1:npoi-1-ijumpb)=rho2i(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=rho2i(1:ijumpb)
    rho2i=ri
    ri(1:npoi-1-ijumpb)=theti(ijumpb+1:npoi-1)
    ri(npoi-ijumpb:npoi-1)=theti(1:ijumpb)
    theti=ri
    deallocate(ri,zi)
    sig=theti(2)-theti(1)
    rho2i(npoi)=rho2i(1)
    theti(npoi)=theti(1)+sign(twopi,sig)
    rho2i(0)=rho2i(npoi-1)
    theti(0)=theti(npoi-1)-sign(twopi,sig)
  end if

  rho2=(r-rc)**2+(z-zc)**2
  thet=atan2(z-zc,r-rc)

  ibeg=0
  iend=npoi
  do while(ibeg+1.lt.iend)
    i=(ibeg+iend)/2
    if((thet-theti(i))*sig.gt.0.d0) then
      ibeg=i
    else
      iend=i
    end if
  end do
  iend=min(iend,npoi-1)
  ibeg=iend-1
  x=theti(ibeg-1:iend+1)
  y=rho2i(ibeg-1:iend+1)

  xx=thet
  yy=y(1)*(xx-x(2))/(x(1)-x(2))*(xx-x(3))/(x(1)-x(3))*(xx-x(4))/(x(1)-x(4)) &
    +y(2)*(xx-x(3))/(x(2)-x(3))*(xx-x(4))/(x(2)-x(4))*(xx-x(1))/(x(2)-x(1)) &
    +y(3)*(xx-x(4))/(x(3)-x(4))*(xx-x(1))/(x(3)-x(1))*(xx-x(2))/(x(3)-x(2)) &
    +y(4)*(xx-x(1))/(x(4)-x(1))*(xx-x(2))/(x(4)-x(2))*(xx-x(3))/(x(4)-x(3))

  if(rho2.gt.yy) then
    incore=-1
    return
  elseif((psif-psi_cut)*sigpsi.lt.0.d0) then
    incore=1
    return
  end if

  incore=0

  call localizer(psi_cut,psi_sep,psif,weight,dweight,ddweight)

  plaf=weight
  dpladr=dweight*dpsidr
  dpladz=dweight*dpsidz
  d2pladr2=ddweight*dpsidr**2+dweight*d2psidr2
  d2pladrdz=ddweight*dpsidr*dpsidz+dweight*d2psidrdz
  d2pladz2=ddweight*dpsidz**2+dweight*d2psidz2

  vacf=1.d0-plaf
  dvacdr=-dpladr
  dvacdz=-dpladz
  d2vacdr2=-d2pladr2
  d2vacdrdz=-d2pladrdz
  d2vacdz2=-d2pladz2

end subroutine inthecore


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine localizer(x1,x2,x,weight,dweight,ddweight)
  use libneo_kinds, only : dp

  implicit none

  real(dp), parameter :: c1=-6.283185307179586d0,c2=-1.414213562373095d0

  real(dp), intent(in) :: x1,x2,x
  real(dp), intent(out) :: weight,dweight,ddweight
  real(dp) :: t,exin

  t=(x-x1)/(x2-x1)

  if(t.le.0.d0) then
    weight=1.d0
    dweight=0.d0
    ddweight=0.d0
  elseif(t.ge.1.d0) then
    weight=0.d0
    dweight=0.d0
    ddweight=0.d0
  else
    exin=exp(c2/t)
    weight=exp(c1/(1.d0-t)*exin)
    dweight=weight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t)
    ddweight=dweight*c1*(1.d0/(1.d0-t)-c2/t**2)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)**2+2.d0*c2/t**3)*exin/(1.d0-t) &
            +weight*c1*(1.d0/(1.d0-t)-c2/t**2)**2*exin/(1.d0-t)
  end if

  dweight=dweight/(x2-x1)
  ddweight=ddweight/(x2-x1)**2

end subroutine localizer


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine window_filter(n,nw,arr_in,arr_out)
  use libneo_kinds, only : dp

  implicit none

  integer, intent(in) :: n,nw
  real(dp), dimension(n), intent(in) :: arr_in
  real(dp), dimension(n), intent(out) :: arr_out
  integer :: nwa,i

  do i=1,n
    nwa=min(nw,i-1,n-i)
    arr_out(i)=sum(arr_in(i-nwa:i+nwa))/(2*nwa+1)
  end do

end subroutine window_filter


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine read_field2(icftype,nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)

  use input_files, only : iunit, pfile
  use libneo_kinds, only : dp
  use math_constants, only : TWOPI

  implicit none

  integer, intent(in) :: icftype, nr, np, nz
  real(dp), intent(out) :: rmin, rmax, pmin, pmax, zmin, zmax
  real(dp), dimension(nr,np,nz), intent(out) :: Br,Bp,Bz

  integer :: i,j,k
  real(dp) :: xdim,zdim,zmid,dum

  open(iunit,file=trim(pfile),status='old',action='read')

  !---Input B      -->T = V*s/m/m
  do j=1,np-1   !only npmax-1 points are given
    do k=nz,1,-1  !reverse order of probe data
      do i=1,nr
        if(icftype.eq.1) then
          ! Old Format
          read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
        elseif(icftype.eq.2) then
          ! New Format
          read(iunit,*) dum,dum,dum,Br(i,j,k),Bp(i,j,k),Bz(i,j,k),dum,dum
        end if

        ! Convert to CGS
        Br(i,j,k) = Br(i,j,k)*1.d4
        Bp(i,j,k) = Bp(i,j,k)*1.d4
        Bz(i,j,k) = Bz(i,j,k)*1.d4
      end do
    end do
  end do
  close(iunit)

  xdim=300.d0
  rmin=100.d0
  rmax=rmin+xdim

  pmin = 0.
  pmax = TWOPI

  zdim=400.d0
  zmid=0.d0
  zmin=zmid - zdim/2.d0
  zmax=zmid + zdim/2.d0

  do i=1,nr
    do k=1,nz
      Br(i,np,k) = Br(i,1,k)
      Bp(i,np,k) = Bp(i,1,k)
      Bz(i,np,k) = Bz(i,1,k)
    end do
  end do
end subroutine read_field2

subroutine read_sizes(nr,np,nz)

  use input_files, only : iunit,pfile

  implicit none
  integer, intent(out) :: nr,np,nz

  open(iunit,file=pfile)
  read(iunit,*) nr,np,nz
  close(iunit)

end subroutine read_sizes

subroutine read_field4(nr,np,nz,rmin,rmax,pmin,pmax,zmin,zmax,Br,Bp,Bz)

  use input_files, only : iunit,pfile
  use libneo_kinds, only : dp

  implicit none

  integer, intent(inout) :: nr,np,nz
  integer :: i,j,k
  real(dp), intent(out) :: rmin,rmax,pmin,pmax,zmin,zmax
  real(dp), dimension(nr,np,nz), intent(out) :: Br,Bp,Bz

  open(iunit,file=pfile)
  read(iunit,*) nr,np,nz
  read(iunit,*) rmin,rmax
  read(iunit,*) pmin,pmax
  read(iunit,*) zmin,zmax
  do i=1,nr
    do j=1,np
      do k=1,nz
        read(iunit,*) Br(i,j,k),Bp(i,j,k),Bz(i,j,k)
      end do
    end do
  end do
  close(iunit)

end subroutine read_field4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine spline_fpol

  use field_eq_mod, only : nrad,hfpol,splfpol
  use libneo_kinds, only : dp
  use spl_three_to_five_sub, only: spl_five_reg

  implicit none

  real(dp), dimension(:), allocatable :: b,c,d,e,f

  allocate(b(nrad),c(nrad),d(nrad),e(nrad),f(nrad))

  hfpol=1.d0/dble(nrad-1)

  call spl_five_reg(nrad,hfpol,splfpol(0,:),b,c,d,e,f)

  splfpol(1,:)=b
  splfpol(2,:)=c
  splfpol(3,:)=d
  splfpol(4,:)=e
  splfpol(5,:)=f

  deallocate(b,c,d,e,f)

end subroutine spline_fpol


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine splint_fpol(x,f,fp)

  use field_eq_mod, only : hfpol,splfpol
  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(out) :: f,fp

  integer :: j,k
  real(dp) :: dx

  k=max(0,int(x/hfpol))
  dx=x-k*hfpol
  k=k+1

  f=splfpol(5,k)
  fp=0.d0
  do j=4,0,-1
    f=f*dx+splfpol(j,k)
    fp=fp*dx+splfpol(j+1,k)*dble(j+1)
  end do

end subroutine splint_fpol

subroutine set_zero(Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, &
  dBzdp, dBzdZ)

  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(out) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, &
    dBpdZ, dBzdR, dBzdp, dBzdZ

  Br=0.d0
  Bp=0.d0
  Bz=0.d0
  dBrdR=0.d0
  dBrdp=0.d0
  dBrdZ=0.d0
  dBpdR=0.d0
  dBpdp=0.d0
  dBpdZ=0.d0
  dBzdR=0.d0
  dBzdp=0.d0
  dBzdZ=0.d0
end subroutine set_zero

subroutine add_scaled(Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, &
  dBzdp, dBzdZ, Br1, Bp1, Bz1, dBrdR1, dBrdp1, dBrdZ1, dBpdR1, dBpdp1, dBpdZ1, &
  dBzdR1, dBzdp1, dBzdZ1, scale)

  use libneo_kinds, only : dp

  implicit none

  real(dp), intent(inout) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ, dBpdR,&
    dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
  real(dp), intent(in) :: Br1, Bp1, Bz1, dBrdR1, dBrdp1, dBrdZ1, dBpdR1, &
    dBpdp1, dBpdZ1, dBzdR1, dBzdp1, dBzdZ1
  real(dp), intent(in) :: scale

  Br=Br+scale*Br1
  Bp=Bp+scale*Bp1
  Bz=Bz+scale*Bz1
  dBrdR=dBrdR+scale*dBrdR1
  dBrdp=dBrdp+scale*dBrdp1
  dBrdZ=dBrdZ+scale*dBrdZ1
  dBpdR=dBpdR+scale*dBpdR1
  dBpdp=dBpdp+scale*dBpdp1
  dBpdZ=dBpdZ+scale*dBpdZ1
  dBzdR=dBzdR+scale*dBzdR1
  dBzdp=dBzdp+scale*dBzdp1
  dBzdZ=dBzdZ+scale*dBzdZ1
end subroutine add_scaled

end module field_sub
