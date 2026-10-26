program hamada_to_boozer
  use binsrc_sub, only: binsrc
  use plag_coeff_sub, only: plag_coeff
  use efit_to_boozer_mod, only: psimax
  use efit_sfl, only: efit_sfl_t, efit_sfl_config_t, efit_sfl_status_t, efit_sfl_surface_t, &
                      contour_t, efit_sfl_symflux, efit_sfl_x_of_s, efit_sfl_raise, &
                      EFIT_SFL_COVARIANT, EFIT_SFL_JACOBIAN_CONST
  use efit_sfl_magdata, only: efit_to_boozer_read_options, efit_sfl_load, efit_sfl_symflux_eval
!
  implicit none
!
  integer, parameter :: nplag=4, nder=0
  double complex, parameter :: imun = (0.d0,1.d0)
!
  integer :: mpolmax,nrhopol,iunit,iunit1,iunit2,irhopol,ntor
  integer :: mdum_half,mdum_full
  integer :: nstep,nlabel,ntheta,nsurf,idummy,ilabel
  integer :: is, ibeg, iend, m, mm, ierr, j
  integer :: construction, write_gs_map
  double precision :: theta, nu_B, nu_H, dnu_B, dummy4(4)
  type(efit_sfl_t) :: sfl_B, sfl_H
  type(efit_sfl_config_t) :: sfl_cfg
  type(efit_sfl_status_t) :: sfl_status
  type(efit_sfl_surface_t) :: sv_B, sv_H
  type(contour_t) :: ct
!
  double precision :: boexp,dummy,raxis,zaxis,fourier_series_factor
  double precision :: hs, rho_pol_s ,twopi, aiota, phimax, theta_B0, theta_H0, step_dft
!
  double precision, dimension(0:nder,nplag) :: coef
!
  double precision, dimension(:),   allocatable :: rho_pol,dummy1d
  double precision, dimension(:),   allocatable :: R_beg,psi_pol,psi_tor
  double precision, dimension(:),   allocatable :: s,R_beg_s
  double precision, dimension(:),   allocatable :: F_B,F_H,dF_B_dphi
  double complex,   dimension(:),   allocatable :: comfac_BH,b_mn_B
  double complex,   dimension(:,:), allocatable :: b_mn,b_mn_hs
  double complex,   dimension(:,:), allocatable :: polfac_B,polfac_H,c_HB
!
  iunit1=1
  iunit2=2
  twopi = 8.d0*atan(1.d0)
!
!--------------------------------------------
! Read MARS data
!
  open(iunit1,file='mars_data_dimensions_and_extra.in')
  read(iunit1,*) mpolmax
  read(iunit1,*) nrhopol
  read(iunit1,*) theta_H0
  read(iunit1,*) theta_B0
  read(iunit1,*) ntor
  close(iunit1)
!
  theta_H0 = theta_H0/twopi
  theta_B0 = theta_B0/twopi
!
  mdum_half = 2*mpolmax+1
  mdum_full = 4*mpolmax+2
!
  allocate(rho_pol(nrhopol),dummy1d(mdum_full))
  allocate(b_mn(nrhopol,-mpolmax:mpolmax))
!
  open(iunit1,file='RMZM_F.OUT')
  read(iunit1,*) boexp
  read(iunit1,*) boexp,boexp,boexp,boexp
  close(iunit1)
!
  open(iunit1,file='TORQUENTV.OUT')
  open(iunit2,file='DATATORQNTV.OUT')
  do irhopol=1,nrhopol
    read(iunit1,*) rho_pol(irhopol)
    read(iunit2,*) dummy1d
    b_mn(irhopol,:)=cmplx(dummy1d(1:mdum_half),dummy1d(mdum_half+1:mdum_full),8)*boexp
  enddo
  close(iunit1)
  close(iunit2)
!
! End read MARS data
!--------------------------------------------
! Read efit_to_boozer data
!
  open(iunit1,file='efit_to_boozer.inp',status='old',action='read')
  read(iunit1,*) nstep     !number of steps for field line integration
  read(iunit1,*) nlabel    !number of flux surfaces (equidistant in starting R) computed by original field line integration
  read(iunit1,*) ntheta    !grid size over poloidal angle for the discrete Fourier transform
  read(iunit1,*) idummy
  read(iunit1,*) nsurf     !number of flux surfaces in Boozer file
  read(iunit1,*) idummy
  read(iunit1,*) psimax    !psi at plasma boundary
! construction 1 (default): Boozer and Hamada angles from the flux-surface contours
! of the G-EQDSK (efit_sfl); 2: legacy, efit_to_boozer output and field-line integration
  call efit_to_boozer_read_options(iunit1, construction, sfl_cfg, write_gs_map)
  close(iunit1)
  if (construction == 1) then
    sfl_cfg%route = EFIT_SFL_COVARIANT
    call efit_sfl_load(sfl_B, sfl_cfg, sfl_status)
    call efit_sfl_raise(sfl_status)
    sfl_H = sfl_B
    sfl_H%cfg%route = EFIT_SFL_JACOBIAN_CONST
  else
!
  open(iunit1,file='box_size_axis.dat')
  read(iunit1,*) dummy
  read(iunit1,*) dummy
  read(iunit1,*) raxis,zaxis  !cylindrical coordinates of the magnetic axis
  close(iunit1)
!
  allocate(R_beg(0:nlabel),psi_pol(0:nlabel),psi_tor(0:nlabel))
  psi_pol(0) = 0.d0
  psi_tor(0) = 0.d0
  R_beg(0) = raxis
!
  open(iunit1,file='flux_functions.dat')
  read(iunit1,*)
  do ilabel=1,nlabel
    read(iunit1,*) R_beg(ilabel),dummy,dummy,psi_pol(ilabel),dummy,psi_tor(ilabel)
  enddo
  close(iunit1)
!
  psi_tor = psi_tor/psi_tor(nlabel)
  endif
!
! End read efit_to_boozer data
!--------------------------------------------
! Re-interpolate Hamada data to equidistant grid in s for Boozer file
!
  allocate(s(nsurf),R_beg_s(nsurf),b_mn_hs(nsurf,-mpolmax:mpolmax))
!
  hs=1.d0/dfloat(nsurf)
!
  do is=1,nsurf
    s(is)=hs*(dfloat(is)-0.5d0)
!
    if (construction == 1) then
      rho_pol_s = sqrt(efit_sfl_x_of_s(sfl_B, s(is)))
    else
    call binsrc(psi_tor(0:nlabel),0,nlabel,s(is),ibeg)
!
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,s(is),psi_tor(ibeg:iend),coef)
!
! poloidal radius and starting R value for s-grid of Bozzer file:
    rho_pol_s = sqrt(abs(sum(coef(0,:)*psi_pol(ibeg:iend))))
    R_beg_s(is) = sum(coef(0,:)*R_beg(ibeg:iend))
    endif
!
!
    call binsrc(rho_pol(1:nrhopol),1,nrhopol,rho_pol_s,ibeg)
!
    ibeg=max(1,ibeg-nplag/2)
    iend=ibeg+nplag-1
    if(iend.gt.nrhopol) then
      iend=nrhopol
      ibeg=iend+1-nplag
    endif
!
    call plag_coeff(nplag,nder,rho_pol_s,rho_pol(ibeg:iend),coef)
!
! Hamada Fourier coefficients for s-grid of Bozzer file:
    b_mn_hs(is,-mpolmax:mpolmax)=matmul(coef(0,:),b_mn(ibeg:iend,-mpolmax:mpolmax))
!
!    write (2001,*) s(is),dble(b_mn_hs(is,:))
!    write (2002,*) s(is),dimag(b_mn_hs(is,:))
  enddo
!
  fourier_series_factor = 1.d0  !1.d0 in case real part is taken for series with n > 0, 2.d0 otherwise
  b_mn_hs = b_mn_hs * fourier_series_factor
!
! End re-interpolate Hamada data to equidistant grid in s for Boozer file
!--------------------------------------------
! Write fake Boozer file
!
  iunit = 71
  dummy = 0.d0
!
  open(iunit,file='pert_fake.bc')
  write(iunit,*) 'CC Fake Boozer-coordinate data file for perturbation field'
  write(iunit,*) 'CC Fourier amplitudes over Boozer angles are set to those over Hamada angles'
  write(iunit,*) 'CC File format is by Andreas Martitsch for use in NEO2-QL'
  write(iunit,*) 'CC Original perturbation is computed by MARS in Hamada coordinates'
  write(iunit,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'
  write(iunit,'(4i6,e15.6,2f10.5)') 2*mpolmax, 0, nsurf, 1, dummy, dummy, dummy
!
  do is=1,nsurf
    write(iunit,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'
    write(iunit,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'
    write(iunit,'(6e17.8)') s(is), dummy, dummy, dummy, dummy, dummy
    write(iunit,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'
    do m=-mpolmax,mpolmax
      write(iunit,'(2i5,8e17.8)') m, ntor, dummy, dummy, dummy, dummy, dummy, dummy,   &
                     dble(b_mn_hs(is,m)), -dimag(b_mn_hs(is,m))
    enddo
  enddo
!
  close(iunit)
! End write fake Boozer file
!--------------------------------------------
! Convert the spectrum and write a proper Boozer file
!
  allocate(F_B(ntheta),F_H(ntheta),dF_B_dphi(ntheta))
  allocate(comfac_BH(ntheta),polfac_H(ntheta,-mpolmax:mpolmax),polfac_B(ntheta,-mpolmax:mpolmax))
  allocate(c_HB(-mpolmax:mpolmax,-mpolmax:mpolmax),b_mn_B(-mpolmax:mpolmax))
!
  open(iunit,file='pert.bc')
  write(iunit,*) 'CC Boozer-coordinate data file for perturbation field'
  write(iunit,*) 'CC Fourier amplitudes over Boozer angles got from conversion of Hamada ones'
  write(iunit,*) 'CC File format is by Andreas Martitsch for use in NEO2-QL'
  write(iunit,*) 'CC Original perturbation is computed by MARS in Hamada coordinates'
  write(iunit,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'
  write(iunit,'(4i6,e15.6,2f10.5)') 2*mpolmax, 0, nsurf, 1, dummy, dummy, dummy
!
  do is=1,nsurf
!
    if (construction == 1) then
! Along the field line from theta = 0, phi = q theta with theta the symmetry-flux angle
! and theta_B,H = theta + iota nu_B,H(theta), normalized as in the legacy branch.
      call efit_sfl_symflux(sfl_B, efit_sfl_x_of_s(sfl_B, s(is)), sfl_B%nt/4, sv_B, ct, sfl_status)
      if (sfl_status%code == 0) &
        call efit_sfl_symflux(sfl_H, efit_sfl_x_of_s(sfl_H, s(is)), sfl_H%nt/4, sv_H, ct, sfl_status)
      call efit_sfl_raise(sfl_status)
      aiota = sv_B%iota
      phimax = twopi/aiota
      do j=1,ntheta
        theta = twopi*dble(j)/dble(ntheta)
        call efit_sfl_symflux_eval(sv_B, theta, dummy4(1), dummy4(2), nu_B, dummy4(3), dnu_B)
        call efit_sfl_symflux_eval(sv_H, theta, dummy4(1), dummy4(2), nu_H, dummy4(3), dummy4(4))
        F_B(j) = (theta + aiota*nu_B)/twopi
        F_H(j) = (theta + aiota*nu_H)/twopi
        dF_B_dphi(j) = aiota*(1.d0 + aiota*dnu_B)/twopi
      enddo
    else
    call field_line_integration_for_coverter(R_beg_s(is),zaxis,ntheta,phimax,F_B,F_H,dF_B_dphi,ierr)
!
    aiota = twopi/phimax
    endif
    step_dft = phimax/dble(ntheta)
    F_B = F_B + theta_B0
    F_H = F_H + theta_H0
!
    comfac_BH = dF_B_dphi * exp(imun*dble(ntor)*phimax*(F_H-F_B))
    polfac_H(:,0) = (1.d0,0.d0)
    polfac_B(:,0) = (1.d0,0.d0)
    polfac_H(:,1) = exp(imun*twopi*F_H)
    polfac_H(:,-1) = conjg(polfac_H(:,1))
    polfac_B(:,-1) = exp(imun*twopi*F_B)
    polfac_B(:,1) = conjg(polfac_B(:,-1))
!
    do m=2,mpolmax
      polfac_H(:,m) = polfac_H(:,m-1)*polfac_H(:,1)
      polfac_H(:,-m) = conjg(polfac_H(:,m))
      polfac_B(:,m) = polfac_B(:,m-1)*polfac_B(:,1)
      polfac_B(:,-m) = conjg(polfac_B(:,m))
    enddo
!
! Conversion matrix Eq.(61)
    do m = -mpolmax,mpolmax
      do mm = -mpolmax,mpolmax
        c_HB(m,mm) = step_dft*sum(comfac_BH*polfac_B(:,m)*polfac_H(:,mm))
      enddo
    enddo
!
! Boozer coefficients Eq.(49)
    b_mn_B = matmul(c_HB,b_mn_hs(is,:))
!
    write(iunit,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'
    write(iunit,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'
    write(iunit,'(6e17.8)') s(is), aiota, dummy, dummy, dummy, dummy
    write(iunit,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'
    do m=-mpolmax,mpolmax
      write(iunit,'(2i5,8e17.8)') m, ntor, dummy, dummy, dummy, dummy, dummy, dummy,   &
                     dble(b_mn_b(m)), -dimag(b_mn_b(m))
    enddo
  enddo
!
  close(iunit)
! End convert the spectrum and write a proper Boozer file
!--------------------------------------------
end program hamada_to_boozer
