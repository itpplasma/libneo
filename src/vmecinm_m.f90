subroutine vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s, &
                    nsurfb,nstrb,kparb,flux)
  !-----
  ! Usage:
  !
  !    call vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s,    &
  !               nsurfm,nstrm,kpar,torflux)
  !
  !  where scalars are:
  !
  !  nstrm - (integer) number of harmonics
  !  kpar  - (integer) number of radial points
  !
  !  vectors are:
  !  aiota(0:kpar)             - (double precision) iota profile
  !  sps(0:kpar) = 0:kpar      - (double precision) radial index as dble number
  !  phi(0:kpar)               - (double precision) toroidal flux
  !  s(0:kpar)                 - (double precision) normalized toroidal flux
  !  axm(nstrm)                - (double precision) poloidal mode numbers
  !  axn(nstrm)                - (double precision) toroidal mode numbers
  !
  !  matrices are:
  !  rmnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of R for various cos(m*theta - n*phi) harmonics
  !  zmnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of Z for various cos(m*theta - n*phi) harmonics
  !  almnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of lambda for various cos(m*theta - n*phi) harmonics
  !  rmns(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of R for various sin(m*theta - n*phi) harmonics
  !  zmns(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of Z for various sin(m*theta - n*phi) harmonics
  !  almns(nstrm,0:kpar)     - (double precision) profiles of Fourier
  !  amplitudes of lambda for various sin(m*theta - n*phi) harmonics

  use libneo_kinds, only : dp
  use math_constants, only : PI
  use new_vmec_stuff_mod, only : netcdffile, vmec_B_scale, vmec_RZ_scale
  use nctools_module, only : nc_open, nc_close, nc_get

  implicit none

  real(dp), parameter :: fac_b0 = 1d4, fac_r0 = 1d2
  real(dp) :: fac_b, fac_r

  integer :: nsurfb,nstrb,kparb,ncid,i
  real(dp) :: flux
  real(dp), dimension(nstrb)         :: axm,axn
  real(dp), dimension(0:kparb)       :: sps,aiota,phi,s
  real(dp), dimension(nstrb,0:kparb) :: rmnc,zmnc,almnc,lmnc
  real(dp), dimension(nstrb,0:kparb) :: rmns,zmns,almns,lmns
  integer :: lasym_int
  logical :: lasym

  associate(dummy => nsurfb)
  end associate

  fac_b = fac_b0 * vmec_B_scale
  fac_r = fac_r0 * vmec_RZ_scale

  do i=0,kparb
    sps(i)=dble(i)
  enddo

  call nc_open(netcdffile, ncid)

  call nc_get(ncid, 'lasym__logical__', lasym_int)
  lasym = (lasym_int == 1)

  call nc_get(ncid, 'phi', phi)
  phi = phi/(2*PI)  ! added by Christopher Albert, 2019-09-16 for correct normalization

  flux=phi(kparb)
  flux=flux*fac_b*fac_r**2
  phi=phi*fac_b*fac_r**2
  s=phi/flux

  call nc_get(ncid, 'xm', axm)
  call nc_get(ncid, 'xn', axn)
  call nc_get(ncid, 'iotaf', aiota)
  call nc_get(ncid, 'rmnc', rmnc)
  call nc_get(ncid, 'zmns', zmns)
  call nc_get(ncid, 'lmns', lmns)
  if (lasym) then
     call nc_get(ncid, 'rmns', rmns)
     call nc_get(ncid, 'zmnc', zmnc)
     call nc_get(ncid, 'lmnc', lmnc)
  else
     rmns = 0d0
     zmnc = 0d0
     lmnc = 0d0
  end if
  ! Convert half-mesh to full mesh for lambda
  ! added by Christopher Albert, 2020-02-11
  almnc(:,0) = 0d0
  almns(:,0) = 0d0
  do i=1,kparb-1
    almnc(:,i) = 0.5d0*(lmnc(:,i+1) + lmnc(:,i))
    almns(:,i) = 0.5d0*(lmns(:,i+1) + lmns(:,i))
  enddo
  almnc(:,kparb) = lmnc(:,kparb) &
    & + 0.5d0*(lmnc(:,kparb) - lmnc(:,kparb-1))

  almns(:,kparb) = lmns(:,kparb) &
    & + 0.5d0*(lmns(:,kparb) - lmns(:,kparb-1))

  rmnc = rmnc*fac_r
  zmnc = zmnc*fac_r
  rmns = rmns*fac_r
  zmns = zmns*fac_r

  call nc_close(ncid)

end subroutine vmecin

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)

  use libneo_kinds, only : dp
  use new_vmec_stuff_mod, only : nper,rmajor

  implicit none

  integer :: L1i
  real(dp) :: RT0,R0i,cbfi,bz0i,bf0

  associate(dummy => R0i)
  end associate
  associate(dummy => cbfi)
  end associate
  associate(dummy => bz0i)
  end associate
  associate(dummy => bf0)
  end associate

  L1i=nper
  RT0=rmajor*1.d2

end subroutine stevvo
