  program efit_to_boozer
!
  use field_sub
  use efit_to_boozer_mod
  use input_files, only : gfile
  use efit_sfl, only : efit_sfl_t, efit_sfl_config_t, efit_sfl_surface_t, efit_sfl_status_t, &
                       efit_sfl_transform, efit_sfl_diag_write, efit_sfl_gs_map, efit_sfl_raise
  use efit_sfl_magdata, only : efit_to_boozer_read_options, efit_sfl_magdata_init
!
  implicit none
  integer, parameter :: cdp = kind((1.0d0, 1.0d0))
!
  integer :: nstep,nsurfmax,i,is,it,nsurf,nt,mpol,inp_label,m,iunit
  integer :: iunit_lhs,iunit_rhs                                                            !<=25.07.2023
  double precision :: s,theta,hs,htheta,aiota,aJb,B_theta,B_phi,oneovernt,twoovernt,sqrtg00
  double precision :: psi,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta,sigma,     &
                      R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta
  double precision :: sigma_lhs,sigma_rhs                                                   !<=25.07.2023
  double precision :: phi,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: dB_phi_ds,dB_theta_ds,pprime
  double complex   :: four_ampl
!
  double precision, dimension(:), allocatable :: R_oft,Z_oft,al_oft,B_oft
  double precision, dimension(:), allocatable :: Rmn_c,Rmn_s,Zmn_c,Zmn_s,almn_c,almn_s,Bmn_c,Bmn_s
  double complex,   dimension(:), allocatable :: calE,calEm_times_Jb
!
! Construction of the angles: 1 = flux-surface contours of the G-EQDSK primaries
! (module efit_sfl, default), 2 = legacy field-line integration.
  integer, parameter :: CONSTRUCTION_CONTOUR = 1, CONSTRUCTION_LEGACY = 2
  integer :: construction,write_gs_map
  type(efit_sfl_t) :: sfl
  type(efit_sfl_status_t) :: sfl_status
  type(efit_sfl_config_t) :: sfl_cfg
  type(efit_sfl_surface_t), dimension(:), allocatable :: sfl_out
  double precision, dimension(:), allocatable :: s_out
  double precision, dimension(:,:), allocatable :: gs_res,gs_mask
  double precision :: flux_hdr,a_hdr,r_hdr,psitor_hdr
  character(len=512) :: override_line
  character(len=128) :: construction_note
!
  open(1,file='efit_to_boozer.inp',status='old',action='read')
  read (1,*) nstep    !number of integration steps
  read (1,*) nlabel   !grid size over radial variable
  read (1,*) ntheta   !grid size over poloidal angle
  read (1,*) nsurfmax !number of starting points between the
                      !magnetic axis and right box boundary
                      !when searching for the separatrix
  read (1,*) nsurf    !number of flux surfaces in Boozer file
  read (1,*) mpol     !number of poloidal modes in Boozer file
  read (1,*) psimax   !psi at plasma boundary
! Optional lines (defaults if absent): construction, route, allow_poor_equilibrium,
! override_reason, n_theta of the contours, write_gs_map.
  call efit_to_boozer_read_options(1, construction, sfl_cfg, write_gs_map)
  close(1)
!
  allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psi_pol(0:nlabel))
  allocate(psi_tor_vac(nlabel),psi_tor(0:nlabel),C_const(nlabel))
!
  allocate(R_spl(0:nspl,0:ntheta,nlabel),Z_spl(0:nspl,0:ntheta,nlabel),bmod_spl(0:nspl,0:ntheta,nlabel))
  allocate(sqgnorm_spl(0:nspl,0:ntheta,nlabel),Gfunc_spl(0:nspl,0:ntheta,nlabel))
!
! Primaries and consistency checks of the G-EQDSK (hard error unless overridden).
! The legacy construction reproduces the previous program without checks.
  if (construction == CONSTRUCTION_CONTOUR) then
    call efit_sfl_magdata_init(sfl, sfl_cfg, sigma, sfl_status)
    call efit_sfl_raise(sfl_status)
  else
  call flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta,      &
                        rmn,rmx,zmn,zmx,raxis,zaxis,sigma, &
                        rbeg,rsmall,qsaf,                  &
                        psi_pol(1:nlabel),psi_tor_vac,     &
                        psi_tor(1:nlabel),C_const,         &
                        R_spl(0,1:ntheta,:),               &
                        Z_spl(0,1:ntheta,:),               &
                        bmod_spl(0,1:ntheta,:),            &
                        sqgnorm_spl(0,1:ntheta,:),         &
                        Gfunc_spl(0,1:ntheta,:))
  endif
!
  call spline_magdata_in_symfluxcoord
!
  print *,'Splining done'
!
  override_line=''
  construction_note=''
  if (construction == CONSTRUCTION_CONTOUR) then
  open(1,file='efit_sfl_diag.dat')
  call efit_sfl_diag_write(sfl, 1)
  close(1)
  if (write_gs_map /= 0) then
    allocate(gs_res(sfl%nw,sfl%nh),gs_mask(sfl%nw,sfl%nh))
    call efit_sfl_gs_map(sfl, gs_res, gs_mask)
    open(1,file='efit_sfl_gsmap.dat')
    write (1,*) '# R [cm], Z [cm], (Delta* psi - gsrhs)/max|gsrhs|, inside psi_edge'
    do it=1,sfl%nh
      do i=1,sfl%nw
        write (1,'(4es24.15)') sfl%g%R_eqd(i), sfl%g%Z_eqd(it), gs_res(i,it), gs_mask(i,it)
      enddo
    enddo
    close(1)
  endif
  if (sfl%diag%override) write (override_line,'(a,i0,1x,es10.3,a,es10.3,a,a)') &
      '; sfl override: ', sfl%diag%status%code, sfl%diag%status%value, ' > ', &
      sfl%diag%status%tol, ': ', trim(sfl_cfg%override_reason)
!
  write (construction_note,'(a,i0,a)') &
      ' (contour SFL route ', sfl_cfg%route, ', pprime = file p''(psi) dpsi/ds)'
  endif
  if (construction == CONSTRUCTION_CONTOUR) then
    allocate(s_out(nsurf),sfl_out(nsurf))
    do is=1,nsurf
      s_out(is)=(dfloat(is)-0.5d0)/dfloat(nsurf)
    enddo
    call efit_sfl_transform(sfl, s_out, mpol, sfl_out)
    psitor_hdr=sfl%psitor_max
    a_hdr=sfl%rminor
    r_hdr=sfl%r0
    sigma=sfl%sig_theta
  else
    psitor_hdr=psitor_max
    a_hdr=rsmall(nlabel)
    r_hdr=raxis
  endif
  flux_hdr=sigma*psitor_hdr*1d-8*twopi
!
  nt=ntheta
!
  inp_label=1
  hs=1.d0/dfloat(nsurf)
  oneovernt=1.d0/dfloat(nt)
  twoovernt=2.d0*oneovernt
  htheta=twopi*oneovernt
!
  allocate(calE(nt),calEm_times_Jb(nt))
  allocate(R_oft(nt),Z_oft(nt),al_oft(nt),B_oft(nt))
  allocate(Rmn_c(0:mpol),Rmn_s(0:mpol),Zmn_c(0:mpol),Zmn_s(0:mpol))
  allocate(almn_c(0:mpol),almn_s(0:mpol),Bmn_c(0:mpol),Bmn_s(0:mpol))
!
!----------------------Begin comments Sergei 25.07.2023---------------------------------------
! Variable "sigma" indicates screw direction for the field line.
!
! If sigma = 1, this is a standard AUG case of the left-hand screw field line.
! For the (standard) left handed coordinates $(s,\vartheta,\varphi)$, iota > 0 in this
! case (contra-variant toroidal and poloidal field components have the same signs).
! The original formal of the output file "fromefit.bc" is the same with standard AUG file.
!
! If sigma = -1, field line is a right-hand srew, and the format of "fromefit.bc" file
! may not confirm with AUG output because it corresponds to the right handed
! coordinates $(s,\vartheta,\varphi)$ such that iota remains positive. Definitions
! of the currents, Jpol and Itor are kept independent of the coordinates system in this file,
! i.e. Jpol is total current through the disc with radius R flowing down (opposite to Z axis).
!
! To avoid confusion of coordinate systems, we write out two more output files,
! "fromefit_neo_lhs.bc" and "fromefit_neo_rhs.bc" with coordinate system indicated in the
! file name (left handed and right handed, respectively). In these files, we define
! poloidal current Jpol differently - it is the same with Jpol in "fromefit.bc" format
! in case of left handed coordinate system and has an opposite sign in case of right handed
! coordinate system so that current in toroidal field coils is always counted in "theta" direction.
! In case of standard AUG configurations (sigma = 1), file "fromefit_neo_lhs.bc" is the same with
! original format file "fromefit.bc".
! In the non-standard case (sigma=-1), all three files are different.
!
! All respective changes in the code are marked with !<=25.07.2023
!----------------------End comments Sergei 25.07.2023---------------------------------------
  sigma_lhs=sigma                                                                           !<=25.07.2023
  sigma_rhs=-sigma                                                                          !<=25.07.2023
  iunit=71
  iunit_lhs=72                                                                              !<=25.07.2023
  iunit_rhs=73                                                                              !<=25.07.2023
  open(iunit,file='fromefit.bc')
  open(iunit_lhs,file='fromefit_neo_lhs.bc')                                                !<=25.07.2023
  open(iunit_rhs,file='fromefit_neo_rhs.bc')                                                !<=25.07.2023
! NEO-2 (neo_sub.f90) and NEO-RT (boozer_read) skip exactly five header lines, so
! the construction and an override are recorded inside the existing CC lines.
  write(iunit,*) 'CC Boozer-coordinate data file generated from EFIT equilibrium file', &
                 trim(construction_note)
  write(iunit_lhs,*) 'CC Boozer-coordinate data file generated from EFIT equilibrium file', &
                     trim(construction_note)
  write(iunit_rhs,*) 'CC Boozer-coordinate data file generated from EFIT equilibrium file', &
                     trim(construction_note)
  write(iunit,*) 'CC Boozer file format: E. Strumberger'
  write(iunit_lhs,*) 'CC Boozer file format: left handed set (s,theta,phi)'                 !<=25.07.2023
  write(iunit_rhs,*) 'CC Boozer file format: right handed set (s,theta,phi)'                !<=25.07.2023
  write(iunit,*) 'CC Authors: S.V. Kasilov, C.G. Albert'
  write(iunit_lhs,*) 'CC Authors: S.V. Kasilov, C.G. Albert'                                !<=25.07.2023
  write(iunit_rhs,*) 'CC Authors: S.V. Kasilov, C.G. Albert'                                !<=25.07.2023
  write(iunit,*) 'CC Original EFIT equilibrium file: ',trim(gfile),trim(override_line)
  write(iunit_lhs,*) 'CC Original EFIT equilibrium file: ',trim(gfile),trim(override_line)
  write(iunit_rhs,*) 'CC Original EFIT equilibrium file: ',trim(gfile),trim(override_line)
  write(iunit,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'
  write(iunit_lhs,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'    !<=25.07.2023
  write(iunit_rhs,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'    !<=25.07.2023
  write(iunit,'(4i6,3es25.16)') mpol, 0, nsurf, 1, flux_hdr, a_hdr*1d-2, r_hdr*1d-2
  write(iunit_lhs,'(4i6,3es25.16)') mpol, 0, nsurf, 1, flux_hdr, a_hdr*1d-2, r_hdr*1d-2  !<=25.07.2023
  write(iunit_rhs,'(4i6,3es25.16)') mpol, 0, nsurf, 1, flux_hdr, a_hdr*1d-2, r_hdr*1d-2  !<=25.07.2023
!
  phi=0.d0
!
  do is=1,nsurf
    s=hs*(dfloat(is)-0.5d0)
    if (construction == CONSTRUCTION_CONTOUR) then
      aiota=sfl_out(is)%iota
      B_phi=sfl_out(is)%curr_g
      B_theta=sfl_out(is)%curr_i
      pprime=sfl_out(is)%pp_file
      sqrtg00=sfl_out(is)%sqrtg00
      Rmn_c=sfl_out(is)%rmnc
      Rmn_s=sfl_out(is)%rmns
      Zmn_c=sfl_out(is)%zmnc
      Zmn_s=sfl_out(is)%zmns
      almn_c=sfl_out(is)%vmnc
      almn_s=sfl_out(is)%vmns
      Bmn_c=sfl_out(is)%bmnc
      Bmn_s=sfl_out(is)%bmns
    else
    sqrtg00=0.d0
    do it=1,nt
      theta=htheta*dfloat(it)
!
      call magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta, &
                                       R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)
!
      call field_eq(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
      calE(it)=exp(cmplx(0.d0,theta+G/q, kind=cdp))
      aJb=R*(Br**2+Bp**2+Bz**2)/(C_norm*Bp)
      sqrtg00=sqrtg00+aJb/bmod**2
      ! Fourier quadrature must use the Jacobian of its interpolated angle.
      calEm_times_Jb(it)=cmplx(1.d0+dG_dtheta/q,0.d0, kind=cdp)
!
      R_oft(it)=R
      Z_oft(it)=Z
      al_oft(it)=-G/twopi
      B_oft(it)=bmod
    enddo
!
! iota, poloidal and toroidal covariant components:
    aiota=1.d0/q
    B_phi=Bp*R
    B_theta=q*(C_norm-B_phi)
    sqrtg00=(B_phi+B_theta/q)*sqrtg00*oneovernt
    dB_phi_ds=(dBpdR*dR_ds+dBpdZ*dZ_ds)*R+Bp*dR_ds
    dB_theta_ds=dq_ds*(C_norm-B_phi)+q*(dC_norm_ds-dB_phi_ds)
    pprime=-(dB_phi_ds+aiota*dB_theta_ds)/(2.d0*twopi*sqrtg00)
!
! Fourier amplitudes of R, Z, lambda and B:
    Rmn_c(0)=real(sum(R_oft*calEm_times_Jb))*oneovernt
    Rmn_s(0)=0.d0
    Zmn_c(0)=real(sum(Z_oft*calEm_times_Jb))*oneovernt
    Zmn_s(0)=0.d0
    almn_c(0)=real(sum(al_oft*calEm_times_Jb))*oneovernt
    almn_s(0)=0.d0
    Bmn_c(0)=real(sum(B_oft*calEm_times_Jb))*oneovernt
    Bmn_s(0)=0.d0
!
    do m=1,mpol
      calEm_times_Jb=calEm_times_Jb*calE
!
      four_ampl=sum(R_oft*calEm_times_Jb)*twoovernt
      Rmn_c(m)=real(four_ampl)
      Rmn_s(m)=dimag(four_ampl)
      four_ampl=sum(Z_oft*calEm_times_Jb)*twoovernt
      Zmn_c(m)=real(four_ampl)
      Zmn_s(m)=dimag(four_ampl)
      four_ampl=sum(al_oft*calEm_times_Jb)*twoovernt
      almn_c(m)=real(four_ampl)
      almn_s(m)=dimag(four_ampl)
      four_ampl=sum(B_oft*calEm_times_Jb)*twoovernt
      Bmn_c(m)=real(four_ampl)
      Bmn_s(m)=dimag(four_ampl)
    enddo
    endif
!
! test Fourier expansion
!do it=1,nt
!  theta=htheta*dfloat(it)
!  bmod=0.d0
!  G=0.d0
!  do m=0,mpol
!    bmod=bmod+almn_c(m)*cos(m*theta)+almn_s(m)*sin(m*theta)
!    G=G+almn_c(m)*cos(m*(theta-aiota*al_oft(it)))+almn_s(m)*sin(m*(theta-aiota*al_oft(it)))
!  enddo
!  write(1001,*) theta,al_oft(it),G,bmod   !columns 2 and 3 should fit, column 4 - series over symmetry angle with Boozer amplitudes
!enddo
!close(1001)
!pause
! end test Fourier expansion
!
    write(iunit,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'
    write(iunit_lhs,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'  !<=25.07.2023
    write(iunit_rhs,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'  !<=25.07.2023
    write(iunit,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'
    write(iunit_lhs,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'  !<=25.07.2023
    write(iunit_rhs,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'  !<=25.07.2023
    ! Radial geometry derivatives require full precision in all surface records.
    write(iunit,'(6es25.16)') s, aiota, -B_phi*5.d0, -sigma*B_theta*5.d0, pprime*1.d-1,              &
                                                      -sigma*sqrtg00*psitor_hdr*1d-6*twopi**2
    write(iunit_lhs,'(6es25.16)') s, aiota*sigma_lhs, -B_phi*5.d0, -sigma*B_theta*5.d0, pprime*1.d-1,            &               !<=25.07.2023
                                                      -abs(sqrtg00*psitor_hdr)*1d-6*twopi**2                                   !<=25.07.2023
    write(iunit_rhs,'(6es25.16)') s, aiota*sigma_rhs, B_phi*5.d0, -sigma*B_theta*5.d0, pprime*1.d-1,             &               !<=25.07.2023
                                                      abs(sqrtg00*psitor_hdr)*1d-6*twopi**2                                    !<=25.07.2023
    write(iunit,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'
    write(iunit_lhs,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &                       !<=25.07.2023
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'                                    !<=25.07.2023
    write(iunit_rhs,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &                       !<=25.07.2023
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'                                    !<=25.07.2023
    do m=0,mpol
      write(iunit,'(2i5,8es25.16)') m,0,Rmn_c(m)*1d-2, Rmn_s(m)*1d-2, Zmn_c(m)*1d-2, Zmn_s(m)*1d-2,  &
                     almn_c(m), almn_s(m), Bmn_c(m)*1d-4, Bmn_s(m)*1d-4
      write(iunit_lhs,'(2i5,8es25.16)') m,0,Rmn_c(m)*1d-2, Rmn_s(m)*1d-2*sigma_lhs, Zmn_c(m)*1d-2, Zmn_s(m)*1d-2*sigma_lhs,  &   !<=25.07.2023
                     almn_c(m), almn_s(m)*sigma_lhs, Bmn_c(m)*1d-4, Bmn_s(m)*1d-4*sigma_lhs                                    !<=25.07.2023
      write(iunit_rhs,'(2i5,8es25.16)') m,0,Rmn_c(m)*1d-2, Rmn_s(m)*1d-2*sigma_rhs, Zmn_c(m)*1d-2, Zmn_s(m)*1d-2*sigma_rhs,  &   !<=25.07.2023
                     almn_c(m), almn_s(m)*sigma_rhs, Bmn_c(m)*1d-4, Bmn_s(m)*1d-4*sigma_rhs                                    !<=25.07.2023
    enddo

  enddo
!
  close(iunit)
  close(iunit_lhs)      !<=25.07.2023
  close(iunit_rhs)      !<=25.07.2023
!
!
!
!
!
  open(1,form='formatted',file='box_size_axis.dat')
  write (1,*) rmn,rmx, '<= rmn, rmx (cm)'
  write (1,*) zmn,zmx, '<= zmn, zmx (cm)'
  write (1,*) raxis,zaxis, '<= raxis, zaxis (cm)'
  close(1)
!
  open(1,form='formatted',file='flux_functions.dat')
  write (1,*) '# R_beg, r,  q, psi_pol, psi_tor_vac, psi_tor'
  do i=1,nlabel
    write (1,*) rbeg(i),rsmall(i),qsaf(i),psi_pol(i),psi_tor_vac(i),psi_tor(i)*psitor_max
  enddo
  close(1)
!
  deallocate(rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor,C_const)
!
  end program efit_to_boozer
