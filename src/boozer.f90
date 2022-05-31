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


  !> \brief Write flux surface header of a boozer file.
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
  subroutine write_boozer_block_head(iunit, s, iota, Jpol_over_nper, Itor, pprime, sqrt_g)
    use libneo_kinds, only : real_kind
    use math_constants, only : TWOPI

    implicit none

    integer, intent(in) :: iunit

    real(kind=real_kind), intent(in) :: s, iota, &
        & Jpol_over_nper, Itor, pprime, sqrt_g

    write(iunit,*) '        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)'
    write(iunit,*) '                                             [A]           [A]             [Pa]         (dV/ds)/nper'
    write(iunit,'(6e17.8)') s, iota, Jpol_over_nper, Itor, &
        & pprime, sqrt_g
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

end module boozer
