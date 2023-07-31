!> Module for handling boozer files - reading, writing, data.
!>
!> Defines subroutines to read/write boozer files and a type for
!> handling the data. Via the type reading/writing can also be done.
!>
!> A boozer file is a textfile that contains data to reconstruct flux
!> surfaces of a magnetic fusion device. The actual data consists of
!> values that are constant on a flux surface, and values that vary on a
!> flux surface. The latter is given in terms of fourier coefficients.
!> The general layout starts with some comment lines. After these there
!> are two lines. The first one is a 'table header' for the values in
!> the second line. These values are 'sizes' for the file and general
!> device parameters.
!> After these initial lines follow the blocks for the flux surfaces.
!> Within the block, first the values constant on the flux surface are
!> given. This is done in three lines, with the first two again act as
!> 'table header', with the names and the units (if any) of the values
!> in the third row. The fourier coefficients of the values that vary on
!> the flux surface, have one line of 'table header', with the names of
!> the variables and in square brackets the unit (if any for the
!> quantity). The first two columns are expected to hold the poloidal
!> and toroidal mode numbers.
module boozer

  implicit none

  type :: BoozerFile
  end type BoozerFile

contains

  subroutine read_boozer_head()
  end subroutine read_boozer_head


  subroutine read_boozer_file()
  end subroutine read_boozer_file


  !> \brief Write header of a boozer file (comments and global variables).
  !>
  !> Write the header of a boozer file to a given file unit.
  !> The header in this context includes comments at the beginning, and
  !> the rows. The first with the 'name' of the global variables and the
  !> second with the values itself.
  !>
  !> Note that this subroutine assumes it gets quantities in SI units.
  !> This is important in two ways. First, the units that are written in
  !> the header are in SI units (meter, Tesla). Second, the quantities
  !> given are written to file as is, i.e. with no conversion.
  !>
  !> input:
  !> ------
  !> iunit: integer, file unit to which to write.
  !> gfile: character array, in case the file is due to conversion from
  !>   another file, should be set to the name of the file.
  !> mpol, ntor, nsurf, nper
  !> flux: float, magnetic flux (unit: Tesla square meter)
  !> rminor: float, minor radius (unit: meter)
  !> rmajor: float, major radius (unit: meter)
  !>
  !> output:
  !> -------
  !> no formal output
  !>
  !> side effects:
  !> -------------
  !> writes data to the file
  subroutine write_boozer_head(iunit, mpol, ntor, nsurf, nper, &
                             & flux, rminor, rmajor, gfile)
    use libneo_kinds, only : real_kind

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: mpol, ntor, nsurf, nper

    real(kind=real_kind), intent(in) :: flux
    real(kind=real_kind), intent(in) :: rminor, rmajor

    character(len=100), intent(in) :: gfile

    write(iunit,*) 'CC Boozer-coordinate data file generated from EFIT equilibrium file'
    write(iunit,*) 'CC Boozer file format: E. Strumberger'
    write(iunit,*) 'CC Authors: S.V. Kasilov, C.G. Albert'
    write(iunit,*) 'CC Original EFIT equilibrium file: ', trim(gfile)
    write(iunit,*) 'm0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'
    write(iunit,'(4i6,e15.6,2f10.5)') mpol, ntor, nsurf, nper, &
        & flux, rminor, rmajor
  end subroutine write_boozer_head


  !> \brief Write flux surface header of a boozer file given the currents.
  !>
  !> Write the header of a boozer flux surface block file to a given
  !> file unit.
  !> The header in this context includes two rows. The first with the
  !> 'name' of the global variables and the second with the values
  !> itself.
  !>
  !> Note that this subroutine assumes it gets quantities in SI units.
  !> This is important in two ways. First, the units that are written in
  !> the header are in SI units. Second, the quantities given are
  !> written to file as is, i.e. with no conversion.
  !>
  !> input:
  !> ------
  !> iunit: integer, file unit to which to write.
  !> s: float, flux surface label (boozer coordinate).
  !> iota: float, rotational transform of the flux surface.
  !> Jpol_over_nper: float, (unit: ampere)
  !> Itor: float, (unit: ampere)
  !> pprime: float, (unit: pascal)
  !> sqert_g: float, (unit: (dV/ds)/nper)
  !>
  !> output:
  !> -------
  !> no formal output
  !>
  !> side effects:
  !> -------------
  !> writes data to the file
  subroutine write_boozer_block_head_current(iunit, s, iota, Jpol_over_nper, Itor, pprime, sqrt_g)
    use libneo_kinds, only : real_kind

    implicit none

    integer, intent(in) :: iunit

    real(kind=real_kind), intent(in) :: s, iota, &
        & Jpol_over_nper, Itor, pprime, sqrt_g

    write(iunit,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'
    write(iunit,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'
    write(iunit,'(6e17.8)') s, iota, Jpol_over_nper, Itor, &
        & pprime, sqrt_g
  end subroutine write_boozer_block_head_current


  !> \brief Write flux surface header of a boozer file given magnetic field components.
  !>
  !> Write the header of a boozer flux surface block file to a given
  !> file unit.
  !> The header in this context includes two rows. The first with the
  !> 'name' of the global variables and the second with the values
  !> itself.
  !>
  !> Note that this subroutine assumes it gets quantities in SI units.
  !> This is important in two ways. First, the units that are written in
  !> the header are in SI units. Second, the quantities given are
  !> written to file as is, i.e. with no conversion.
  !>
  !> input:
  !> ------
  !> iunit: integer, file unit to which to write.
  !> s: float, flux surface label (boozer coordinate).
  !> iota: float, rotational transform of the flux surface.
  !> bsubvB: float, magnetic field component
  !> bsubuB: float, magnetic field component
  !> pprime: float, (unit: pascal)
  !> sqert_g: float, (unit: (dV/ds)/nper)
  !>
  !> output:
  !> -------
  !> no formal output
  !>
  !> side effects:
  !> -------------
  !> writes data to the file
  subroutine write_boozer_block_head(iunit, s, iota, bsubvB, bsubuB, pprime, vp, enfp)
    use libneo_kinds, only : real_kind
    use math_constants, only : MU_0, PI

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: enfp
    real(kind=real_kind), intent(in) :: s, iota, bsubvB, bsubuB, &
        & pprime, vp

    call write_boozer_block_head_current(iunit, s, iota, &
        & -2.0*PI/MU_0*bsubvB/enfp, -2.0*PI/MU_0*bsubuB, pprime, &
        & -4.0*PI**2*vp/enfp)
  end subroutine write_boozer_block_head


  !> \brief Write flux surface data of a boozer file.
  !>
  !> First write column 'names' (with unit in brackets if any), and then
  !> the fourier coefficients describing the flux surface.
  !>
  !> Note that this subroutine assumes it gets quantities in SI units.
  !> This is important in two ways. First, the units that are written in
  !> the header are in SI units (meter and Tesla). Second, the
  !> quantities given are written to file as is, i.e. with no
  !> conversion.
  !>
  !> input:
  !> ------
  !> iunit: integer, file unit to which to write.
  !> total_number_modes: integer, determines the size of the arrays.
  !> m: array of integers, contains the poloidal mode numbers of the
  !>   fourier coefficients.
  !> n: array of integers, as m, but toroidal mode numbers.
  !> Rmn_c: array of floats, cosine coefficients for radial position of
  !>   the flux surface (unit: meters).
  !> Rmn_s: array of floats, as Rmnc_s but sine coefficients.
  !> Zmn_c, Zmn_s: arrays of floats, same as Rmn_c, Rmn_s respectively,
  !>   but for vertical position of the flux surface (unit: meters).
  !> almn_c, almn_s: arrays of floats, same as Rmn_c, Rmn_s respectively,
  !>   but for normalized rotation velocity (no unit).
  !> Bmn_c, Bmn_s: arrays of floats, same as Rmn_c, Rmn_s respectively,
  !>   but for magnitude of magnetic field (unit: Tesla).
  !>
  !> output:
  !> -------
  !> no formal output
  !>
  !> side effects:
  !> -------------
  !> writes data to the file
  subroutine write_boozer_block_data(iunit, total_number_modes, &
      & m, n, Rmn_c, Rmn_s, Zmn_c, Zmn_s, almn_c, almn_s, Bmn_c, Bmn_s)
    use libneo_kinds, only : real_kind

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: total_number_modes
    integer, dimension(1:total_number_modes), intent(in) :: m, n
    real(kind=real_kind), dimension(1:total_number_modes), intent(in) :: &
        & Rmn_c, Rmn_s, Zmn_c, Zmn_s, almn_c, almn_s, Bmn_c, Bmn_s

    integer k

    write(iunit,*) '    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]'   &
                   //'         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'
    do k=1,total_number_modes
      write(iunit,'(2i5,8e17.8)') m(k), n(k), &
          & Rmn_c(k), Rmn_s(k), &
          & Zmn_c(k), Zmn_s(k), &
          & almn_c(k), almn_s(k), &
          & Bmn_c(k), Bmn_s(k)
    end do

  end subroutine write_boozer_block_data


  subroutine write_boozer_block_data_complex(iunit, total_number_modes, m, n, rmn, zmn, vmn, bmn, enfp)
    use libneo_kinds, only : complex_kind

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: total_number_modes, enfp
    integer, dimension(1:total_number_modes), intent(in) :: m, n
    complex(kind=complex_kind), dimension(1:total_number_modes), intent(in) :: &
        & rmn, zmn, vmn, bmn

    call write_boozer_block_data(iunit, total_number_modes, m, n/enfp, &
        & real(rmn), -aimag(rmn), real(zmn), -aimag(zmn), &
        & real(vmn), -aimag(vmn), real(bmn), -aimag(bmn))
  end subroutine write_boozer_block_data_complex


  !> \brief Convert an efit file to a boozer file.
  !>
  !> input:
  !> ------
  !> namelist_unit: integer, optional, file unit number from which to
  !>   read the namelist ('efit_to_boozer_nml') of this subroutine. If
  !>   not given, a file 'efit_to_boozer.inp' is expected to exist and
  !>   contain the namelist.
  subroutine efit_to_boozer(namelist_unit)

    use libneo_kinds, only : complex_kind, real_kind
!~     use efit_to_boozer_mod
!~     use input_files, only : gfile

    implicit none

    integer, intent(in), optional :: namelist_unit

    character(len=1024) :: gfile

    integer :: namelist_unit_loc
    integer :: nstep,nsurfmax,i,is,it,nsurf,nt,mpol,inp_label,m,iunit
    integer, dimension(:), allocatable :: m_c, n_c
    real(kind=real_kind) :: s,theta,hs,htheta,aiota,aJb,B_theta,B_phi,oneovernt,twoovernt,sqrtg00
    real(kind=real_kind) :: psi,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta,sigma,     &
                        R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta
    real(kind=real_kind) :: phi,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(kind=real_kind) :: dB_phi_ds,dB_theta_ds,pprime
    complex(kind=complex_kind)   :: four_ampl

    real(kind=real_kind), dimension(:), allocatable :: R_oft,Z_oft,al_oft,B_oft
    real(kind=real_kind), dimension(:), allocatable :: Rmn_c,Rmn_s,Zmn_c,Zmn_s,almn_c,almn_s,Bmn_c,Bmn_s
    complex(kind=complex_kind),   dimension(:), allocatable :: calE,calEm_times_Jb

!~   module efit_to_boozer_mod
    integer, parameter :: nspl  = 3 !spline order in poloidal angle interpolation
    integer, parameter :: nplag = 4 !stencil for Largange polynomial interpolation (polynomial order + 1)
    integer, parameter :: nder  = 1 !number of derivatives from Largange polynomial interpolation
    logical :: load=.true.
    integer :: nlabel,ntheta
    real(kind=real_kind) :: rmn,rmx,zmn,zmx,raxis,zaxis,h_theta,twopi,psipol_max,psitor_max
    real(kind=real_kind), dimension(nplag)        :: R_lag,Z_lag,sqrtg_lag,bmod_lag,G_lag
    real(kind=real_kind), dimension(nplag)        :: dbmod_dt_lag,dR_dt_lag,dZ_dt_lag,dG_dt_lag
    real(kind=real_kind), dimension(0:nder,nplag) :: coef
    real(kind=real_kind), dimension(:),     allocatable :: rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor,C_const
    real(kind=real_kind), dimension(:,:,:), allocatable :: R_spl,Z_spl,bmod_spl,sqgnorm_spl,Gfunc_spl
!~   end module efit_to_boozer_mod

    namelist /efit_to_boozer_nml/ nstep, nlabel, ntheta, nsurfmax, nsurf, mpol

    if (present(namelist_unit)) then
      namelist_unit_loc = namelist_unit
    else
      open(newunit=namelist_unit_loc, file='efit_to_boozer.inp', status='old')
    end if

    read(namelist_unit_loc, nml=efit_to_boozer_nml)
    close(namelist_unit_loc)

    allocate(rbeg(nlabel),rsmall(nlabel),qsaf(nlabel),psi_pol(0:nlabel))
    allocate(psi_tor_vac(nlabel),psi_tor(0:nlabel),C_const(nlabel))

    allocate(R_spl(0:nspl,0:ntheta,nlabel),Z_spl(0:nspl,0:ntheta,nlabel),bmod_spl(0:nspl,0:ntheta,nlabel))
    allocate(sqgnorm_spl(0:nspl,0:ntheta,nlabel),Gfunc_spl(0:nspl,0:ntheta,nlabel))

!~     call flint_for_Boozer(nstep,nsurfmax,nlabel,ntheta,      & !~ found in field_line_int_for_boozer.f90
!~                           rmn,rmx,zmn,zmx,raxis,zaxis,sigma, &
!~                           rbeg,rsmall,qsaf,                  &
!~                           psi_pol(1:nlabel),psi_tor_vac,     &
!~                           psi_tor(1:nlabel),C_const,         &
!~                           R_spl(0,1:ntheta,:),               &
!~                           Z_spl(0,1:ntheta,:),               &
!~                           bmod_spl(0,1:ntheta,:),            &
!~                           sqgnorm_spl(0,1:ntheta,:),         &
!~                           Gfunc_spl(0,1:ntheta,:))

!~     call spline_magdata_in_symfluxcoord !~ found in spline_and_interpolate_magdata.f90

    print *,'Splining done'

    nt=ntheta

    inp_label=1
    hs=1.d0/real(nsurf, kind=real_kind)
    oneovernt=1.d0/real(nt, kind=real_kind)
    twoovernt=2.d0*oneovernt
    htheta=twopi*oneovernt

    allocate(calE(nt),calEm_times_Jb(nt))
    allocate(R_oft(nt),Z_oft(nt),al_oft(nt),B_oft(nt))
    allocate(m_c(0:mpol))
    m_c = (/ (m,m=0,mpol) /)
    n_c = (/ (0,m=0,mpol) /)
    allocate(Rmn_c(0:mpol),Rmn_s(0:mpol),Zmn_c(0:mpol),Zmn_s(0:mpol))
    allocate(almn_c(0:mpol),almn_s(0:mpol),Bmn_c(0:mpol),Bmn_s(0:mpol))

    open(newunit=iunit,file='fromefit.bc')
    call write_boozer_head(iunit, mpol, 0, nsurf, 1, &
        & sigma*psitor_max*1d-8*twopi, rsmall(nlabel)*1d-2, raxis*1d-2, &
        & trim(gfile))

    phi=0.d0

    do is=1,nsurf
      s=hs*(real(is, kind=real_kind)-0.5d0)
      do it=1,nt
        theta=htheta*real(it, kind=real_kind)

!~         call magdata_in_symfluxcoord_ext(inp_label,s,psi,theta,q,dq_ds,C_norm,dC_norm_ds,sqrtg,bmod,dbmod_dtheta, & !~ found in spline_and_interpolate_magdata.f90
!~                                        R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,G,dG_ds,dG_dtheta)

!~         call field_eq(R,phi,Z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ) !~ found in field_divB0.f90

        calE(it)=exp(cmplx(0.d0,theta+G/q, kind=complex_kind))
        aJb=R*(Br**2+Bp**2+Bz**2)/(C_norm*Bp)
        calEm_times_Jb(it)=cmplx(aJb,0.d0, kind=complex_kind)

        R_oft(it)=R
        Z_oft(it)=Z
        al_oft(it)=-G/twopi
        B_oft(it)=bmod
      end do

      ! iota, poloidal and toroidal covariant components:
      aiota=1.d0/q
      B_phi=Bp*R
      B_theta=q*(C_norm-B_phi)
      sqrtg00=(B_phi+B_theta/q)*real(sum(calEm_times_Jb/B_oft**2))*oneovernt
      dB_phi_ds=(dBpdR*dR_ds+dBpdZ*dZ_ds)*R+Bp*dR_ds
      dB_theta_ds=dq_ds*(C_norm-B_phi)+q*(dC_norm_ds-dB_phi_ds)
      pprime=-(dB_phi_ds+aiota*dB_theta_ds)/(2.d0*twopi*sqrtg00)

      ! Fourier amplitudes of R, Z, lambda and B:
      Rmn_c(0)=real(sum(R_oft*calEm_times_Jb))*oneovernt
      Rmn_s(0)=0.d0
      Zmn_c(0)=real(sum(Z_oft*calEm_times_Jb))*oneovernt
      Zmn_s(0)=0.d0
      almn_c(0)=real(sum(al_oft*calEm_times_Jb))*oneovernt
      almn_s(0)=0.d0
      Bmn_c(0)=real(sum(B_oft*calEm_times_Jb))*oneovernt
      Bmn_s(0)=0.d0

      do m=1,mpol
        calEm_times_Jb=calEm_times_Jb*calE

        four_ampl=sum(R_oft*calEm_times_Jb)*twoovernt
        Rmn_c(m)=real(four_ampl)
        Rmn_s(m)=aimag(four_ampl)
        four_ampl=sum(Z_oft*calEm_times_Jb)*twoovernt
        Zmn_c(m)=real(four_ampl)
        Zmn_s(m)=aimag(four_ampl)
        four_ampl=sum(al_oft*calEm_times_Jb)*twoovernt
        almn_c(m)=real(four_ampl)
        almn_s(m)=aimag(four_ampl)
        four_ampl=sum(B_oft*calEm_times_Jb)*twoovernt
        Bmn_c(m)=real(four_ampl)
        Bmn_s(m)=aimag(four_ampl)
      end do

      call write_boozer_block_head_current(iunit, s, aiota, -B_phi*5.d0, &
          & -sigma*B_theta*5.d0, pprime*1.d-1, -sigma*sqrtg00*psitor_max*1d-6*twopi**2)
      call write_boozer_block_data(iunit, mpol+1, m_c, n_c, &
          & Rmn_c*1d-2, Rmn_s*1d-2, Zmn_c*1d-2, Zmn_s*1d-2, &
          & almn_c, almn_s, Bmn_c*1d-4, Bmn_s*1d-4)
    end do

    close(iunit)

    open(1,form='formatted',file='box_size_axis.dat')
    write (1,*) rmn,rmx, '<= rmn, rmx (cm)'
    write (1,*) zmn,zmx, '<= zmn, zmx (cm)'
    write (1,*) raxis,zaxis, '<= raxis, zaxis (cm)'
    close(1)

    open(1,form='formatted',file='flux_functions.dat')
    write (1,*) '# R_beg, r,  q, psi_pol, psi_tor_vac, psi_tor'
    do i=1,nlabel
      write (1,*) rbeg(i),rsmall(i),qsaf(i),psi_pol(i),psi_tor_vac(i),psi_tor(i)*psitor_max
    end do
    close(1)

    deallocate(rbeg,rsmall,qsaf,psi_pol,psi_tor_vac,psi_tor,C_const)

  end subroutine efit_to_boozer


  !> \brief Get data required for boozer header from vmec file.
  subroutine getHeadDataVmecNc(infilename, nfp, psi_tor_a, a, R0, m0b, n0b)
    use libneo_kinds, only : real_kind
    use nctools_module, only : nc_get, nc_inq_dim, nc_open

    implicit none

    character(len=100), intent(in) :: infilename
    real(kind=real_kind), intent(out) :: psi_tor_a, a, R0
    integer, intent(out) :: nfp, m0b, n0b

    integer :: fileunit
    integer :: empol, entor
    integer :: size_m, size_n
    integer, dimension(:), allocatable :: m,n
    real(kind=real_kind), dimension(:), allocatable :: phipf


    call nc_open(infilename, fileunit)

    call nc_inq_dim(fileunit, 'xm', size_m)
    call nc_inq_dim(fileunit, 'xm', size_n)

    allocate(m(size_m), n(size_n))

    call nc_get(fileunit, 'nfp', nfp)
    call nc_get(fileunit, 'Aminor_p', a)
    call nc_get(fileunit, 'Rmajor_p', R0)
    call nc_get(fileunit, 'phipf', phipf)
    call nc_get(fileunit, 'xm', m)
    call nc_get(fileunit, 'xn', n)
    psi_tor_a = phipf(1)
    empol = maxval(abs(m))
    entor = maxval(abs(n))
    m0b = 2*empol
    n0b = 2*entor

    if (allocated(m)) deallocate(m)
    if (allocated(n)) deallocate(n)

    close(unit=fileunit)

  end subroutine getHeadDataVmecNc

end module boozer
