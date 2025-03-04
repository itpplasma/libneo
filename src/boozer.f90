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
  use math_constants, only : twopi

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
    use libneo_kinds, only : dp

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: mpol, ntor, nsurf, nper

    real(dp), intent(in) :: flux
    real(dp), intent(in) :: rminor, rmajor

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
    use libneo_kinds, only : dp

    implicit none

    integer, intent(in) :: iunit

    real(dp), intent(in) :: s, iota, &
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
    use libneo_kinds, only : dp
    use math_constants, only : MU_0, PI

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: enfp
    real(dp), intent(in) :: s, iota, bsubvB, bsubuB, &
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
    use libneo_kinds, only : dp

    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: total_number_modes
    integer, dimension(1:total_number_modes), intent(in) :: m, n
    real(dp), dimension(1:total_number_modes), intent(in) :: &
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
        & rmn%re, -rmn%im, zmn%re, -zmn%im, &
        & vmn%re, -vmn%im, bmn%re, -bmn%im)
  end subroutine write_boozer_block_data_complex

  !> \brief Get data required for boozer header from vmec file.
  subroutine getHeadDataVmecNc(infilename, nfp, psi_tor_a, a, R0, m0b, n0b)
    use libneo_kinds, only : dp
    use nctools_module, only : nc_get, nc_inq_dim, nc_open

    implicit none

    character(len=100), intent(in) :: infilename
    real(dp), intent(out) :: psi_tor_a, a, R0
    integer, intent(out) :: nfp, m0b, n0b

    integer :: fileunit
    integer :: empol, entor
    integer :: size_m, size_n
    integer, dimension(:), allocatable :: m,n
    real(dp), dimension(:), allocatable :: phipf


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
