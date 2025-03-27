!> \brief Module for handling input and output.
!>
!> This module is intended for handling the input and output of code.
!> This may include interfaces to handle input/output independent of the
!> requested/selected file format.
!> It includes routines to read in (and maybe write) specific file
!> formats, like efit and boozer. Note that the input is just read, not
!> processed.
module io
  use libneo_kinds, only : dp

  implicit none

  private

  public get_free_unit, efit_data_type, boozer_data_type

  !> Storing last used file unit number. Used for getting a free one.
  integer, save :: iunit=100

  ! Formats used for reading the data.
  character(len=*), parameter :: format_efit_header = '(6a8,3i4)'
  character(len=*), parameter :: format_five_rows_doubles = '(5(e16.9))'
  character(len=*), parameter :: format_boozer_output_data = &
    & '(2(i4, 1X), 8(e16.9, 1X))'
  character(len=*), parameter :: format_boozer_head_param = &
    & '(4(1X, i5), 3(1X, e16.9))'
  character(len=*), parameter :: format_boozer_flux_head = &
    & '(6(1X, e16.9))'

  !> \brief Class representing efit data (-file).
  !>
  !> Supports reading the data from a file, writing it to a file and
  !> getting some of the scalars.
  type efit_data_type
    private

    integer :: nwEQD, nhEQD
    integer :: n_bndyxy, nlimEQD

    real(dp) :: psiSep, bt0, rzero

    real(dp), dimension(:), allocatable :: fpol, pres
    real(dp), dimension(:), allocatable :: ffprim
    real(dp), dimension(:), allocatable :: pprime
    real(dp), dimension(:), allocatable :: qpsi
    real(dp), dimension(:,:), allocatable :: psiRZ
    real(dp), dimension(:), allocatable :: LCFS, limEQD
    ! These two are for storing the grid coordinates. Regular -> only one dimension.
    real(dp), dimension(:), allocatable :: rad, zet

    real(dp) :: xdim,zdim,r1,zmid,rmaxis,zmaxis
    real(dp) :: plas_cur, psiAxis

    !> This array of strings collects the part at the beginning of the
    !> first line. It is not clear if this contains some relevant
    !> information or not, so we store it to be save (it is just 60
    !> bytes per file and it is not expected that there will be more
    !> than one or two around).
    character(len=10) :: dummy(6)
  contains
    procedure :: read_data => read_data_of_efit_file
    procedure :: read_dimension => read_dimension_of_efit_file
    procedure :: write_data => write_data_of_efit_file

    procedure :: get_nwEQD => get_nwEQD_
    procedure :: get_nhEQD => get_nhEQD_
    procedure :: get_psiSep => get_psiSep_
    procedure :: get_bt0 => get_bt0_
    procedure :: get_rzero => get_rzero_
    procedure :: get_rad => get_rad_
    procedure :: get_zet => get_zet_
    procedure :: get_psiRZ => get_psiRZ_

    final :: finalize_efit_class_object
  end type efit_data_type

  !> \brief Class representing boozer data (-file).
  !>
  !> Supports reading the data from a file, writing it to a file and
  !> getting some of the scalars.
  type boozer_data_type
    private

    integer :: m0b, n0b, nsurf, nper
    real(dp) :: flux ![Tm^2]
    real(dp) :: a ![m]
    real(dp) :: R ![m]

    real(dp), dimension(:), allocatable :: s, iota
    real(dp), dimension(:), allocatable :: Jpol_nper ! [A]
    real(dp), dimension(:), allocatable :: Itor ![A]
    real(dp), dimension(:), allocatable :: pprime ![Pa]
    real(dp), dimension(:), allocatable :: sqrt_g00 !(dV/ds)/nper

    integer(dp), dimension(:,:), allocatable :: m, n
    real(dp), dimension(:,:), allocatable :: rmnc, rmns, zmnc, zmns ! [m]
    real(dp), dimension(:,:), allocatable :: vmnc,vmns  ! [ ]
    real(dp), dimension(:,:), allocatable :: bmnc, bmns ! [T]
  contains
    procedure :: read_data => read_data_of_boozer_file
    procedure :: write_data => write_data_of_boozer_file

    procedure :: get_m0b => get_m0b_
    procedure :: get_n0b => get_n0b_
    procedure :: get_nsurf => get_nsurf_
    procedure :: get_nper => get_nper_
    procedure :: get_flux => get_flux_
    procedure :: get_a => get_a_
    procedure :: get_R => get_R_
  end type boozer_data_type

contains

  !> \brief Get a free unit number.
  !>
  !> This function will return a unit number, not attached to a file.
  !>
  !> At the moment subsequent calls will only result in values equal to
  !> or larger than previous values, but do not rely on this behaviour,
  !> as it may change.
  !>
  !> \return A unit number (integer) which is not attached to a file.
  function get_free_unit()
    integer :: get_free_unit

    logical :: unit_numer_is_connected_to_file

    unit_numer_is_connected_to_file = .true.
    do
      inquire(unit=iunit, opened=unit_numer_is_connected_to_file)
      if(.not. unit_numer_is_connected_to_file) exit
      iunit = iunit + 1
    end do

    get_free_unit = iunit
  end function get_free_unit

  !> \brief Read first dimensions from efit file, i.e. those in the first line.
  !>
  !> This subroutine will read the two dimensions written in an efit
  !> file with given filename and return them.
  !>
  !> \param filename: input, name of the file from which to read the
  !>   dimensions.
  !> \param nwEQD, nhEQD: output, the dimensions read from the file.
  subroutine read_dimension_of_efit_file(this, filename)

    implicit none

    class(efit_data_type), intent(inout) :: this

    character(len=*), intent(in) :: filename

    integer :: idum, i
    character(len=10) :: dummy(6)


    open(unit=iunit,file=trim(filename),status='old',action='read')
    read(iunit,2000)(dummy(i),i=1,6), idum, this%nwEQD, this%nhEQD
    close(iunit)

    return

    2000  format(6a8,3i4)
  end subroutine read_dimension_of_efit_file

  !> \brief Read content from a efit file with given filename.
  !>
  !> Note that the two arrays LCFS and limEQD will only be allocated, if
  !> the corresponding sizes are larger than zero.
  !>
  !> \param filename: input, name of the file from which to read the
  !>   data.
  subroutine read_data_of_efit_file(this, filename)

    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(inout) :: this

    character(len=*), intent(in) :: filename

    ! loop variables
    integer :: i,j

    ! File unit.
    integer :: gunit

    real (dp) :: xdum
    integer :: idum

    gunit = get_free_unit()

    open(unit=gunit,file=trim(filename),status='old',action='read')

    ! Read in first line with sizes.
    read(gunit,fmt=format_efit_header) (this%dummy(i),i=1,6), idum, this%nwEQD, this%nhEQD

    ! Read first four lines with values:
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & this%xdim, this%zdim, this%rzero, this%r1, this%zmid
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & this%rmaxis, this%zmaxis, this%psiAxis, this%psiSep, this%bt0
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & this%plas_cur, this%psiAxis, xdum, this%rmaxis, xdum
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & this%zmaxis, xdum, this%psiSep, xdum, xdum

    ! Allocate first set of arrays. Done after reading first four lines,
    ! in case there was a problem, arrays are not already allocted.
    allocate(this%fpol(this%nwEQD), this%pres(this%nwEQD), this%ffprim(this%nwEQD))
    allocate(this%pprime(this%nwEQD), this%qpsi(this%nwEQD))
    allocate(this%psiRZ(this%nwEQD,this%nhEQD))

    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & (this%fpol(i),i=1,this%nwEQD)
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & (this%pres(i),i=1,this%nwEQD)
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & (this%ffprim(i),i=1,this%nwEQD)
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & (this%pprime(i),i=1,this%nwEQD)
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & ((this%psiRZ(i,j),i=1,this%nwEQD),j=1,this%nhEQD)
    read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
      & (this%qpsi(i),i=1,this%nwEQD)

    ! Boundary Data, first read size, allocate the arrays, last read the data.
    read(gunit,*,end=55,err=250) this%n_bndyxy, this%nlimEQD

    if (this%n_bndyxy > 0) then
      allocate(this%LCFS(2*this%n_bndyxy))
      read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
        & (this%LCFS(i),i=1,2*this%n_bndyxy)
    end if

    if (this%nlimEQD > 0) then
      allocate(this%limEQD(2*this%nlimEQD))
      read(gunit,fmt=format_five_rows_doubles,end=55,err=250) &
        & (this%limEQD(i),i=1,2*this%nlimEQD)
    end if

    close(gunit)

    allocate(this%rad(this%nwEQD))
    call set_array_equidistant(this%nwEQD, this%xdim, this%r1, this%rad)

    allocate(this%zet(this%nhEQD))
    call set_array_equidistant(this%nhEQD, this%zdim, this%zmid - this%zdim/2.0, this%zet)

    ! Without the return, the print statements for the error labels
    ! would lead to compiler errors.
    return

    ! Error labels
    55    print *, 'Error in read_data_of_efit_file: Early EOF in ',trim(filename); STOP
    250   print *, 'Error in read_data_of_efit_file: Error reading ',trim(filename); STOP

  end subroutine read_data_of_efit_file

  !> \brief Write data of efit class to file with given filename.
  !>
  !> \param filename: input, data is written to file with this name.
  subroutine write_data_of_efit_file(this, filename)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(inout) :: this

    character(len=*), intent(in) :: filename

    ! loop variables
    integer :: i,j

    ! File unit.
    integer :: gunit

    real (dp), parameter :: xdum = 0.0
    integer, parameter :: idum = 0

    gunit = get_free_unit()

    open(unit=gunit,file=trim(filename),status='replace',action='write')

    write (gunit,fmt=format_efit_header) &
      & (this%dummy(i),i=1,6), idum, this%nwEQD, this%nhEQD
    write (gunit,fmt=format_five_rows_doubles) &
      & this%xdim, this%zdim, this%rzero, this%r1, this%zmid
    write (gunit,fmt=format_five_rows_doubles) &
      & this%rmaxis, this%zmaxis, this%psiAxis, this%psiSep, this%bt0
    write (gunit,fmt=format_five_rows_doubles) &
      & this%plas_cur, this%psiAxis, xdum, this%rmaxis, xdum
    write (gunit,fmt=format_five_rows_doubles) &
      & this%zmaxis, xdum, this%psiSep, xdum, xdum
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%fpol(i),i=1,this%nwEQD)
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%pres(i),i=1,this%nwEQD)
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%ffprim(i),i=1,this%nwEQD)
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%pprime(i),i=1,this%nwEQD)
    write (gunit,fmt=format_five_rows_doubles) &
      & ((this%psirz(i,j),i=1,this%nwEQD),j=1,this%nhEQD)
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%qpsi(i),i=1,this%nwEQD)
    write (gunit,fmt='(2i5)') this%n_bndyxy, this%nlimEQD
    write (gunit,fmt=format_five_rows_doubles) &
      & (this%LCFS(i),i=1,2*this%n_bndyxy)
    write (gunit,fmt=format_five_rows_doubles) &
      (this%limEQD(i),i=1,2*this%nlimEQD)

    close(gunit)

  end subroutine write_data_of_efit_file

  function get_nwEQD_(this)
    implicit none

    class(efit_data_type), intent(in) :: this

    integer :: get_nwEQD_

    get_nwEQD_ = this%nwEQD
  end function get_nwEQD_

  function get_nhEQD_(this)
    implicit none

    class(efit_data_type), intent(in) :: this

    integer :: get_nhEQD_

    get_nhEQD_ = this%nhEQD
  end function get_nhEQD_

  function get_psiSep_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_psiSep_

    get_psiSep_ = this%psiSep
  end function get_psiSep_

  function get_bt0_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_bt0_

    get_bt0_ = this%bt0
  end function get_bt0_

  function get_rzero_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_rzero_

    get_rzero_ = this%rzero
  end function get_rzero_

  function get_rad_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_rad_(this%nwEQD)

    get_rad_ = this%rad
  end function get_rad_

  function get_zet_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_zet_(this%nhEQD)

    get_zet_ = this%zet
  end function get_zet_

  function get_psiRZ_(this)
    use libneo_kinds, only : dp

    implicit none

    class(efit_data_type), intent(in) :: this

    real(dp) :: get_psiRZ_(this%nwEQD,this%nhEQD)

    get_psiRZ_ = this%psiRZ
  end function get_psiRZ_

  !> The array will contain values from origin up to origin+width (end
  !> points included).
  !> \note It is called origin as offset is a keyword.
  subroutine set_array_equidistant(number_of_points, width, origin, array)
    use libneo_kinds, only : dp

    implicit none

    integer, intent(in) :: number_of_points
    real(dp), intent(in) :: width, origin
    real(dp), intent(out) :: array(number_of_points)

    integer :: j

    do j=1, number_of_points
      array(j) = origin + (j-1)*(width/(number_of_points-1))
    end do

  end subroutine set_array_equidistant

  !> \brief Destructor for efit class.
  subroutine finalize_efit_class_object(this)
    type (efit_data_type) :: this

    if (allocated(this%fpol)) deallocate(this%fpol)
    if (allocated(this%pres)) deallocate(this%pres)
    if (allocated(this%ffprim)) deallocate(this%ffprim)
    if (allocated(this%pprime)) deallocate(this%pprime)
    if (allocated(this%qpsi)) deallocate(this%qpsi)
    if (allocated(this%rad)) deallocate(this%rad)
    if (allocated(this%psiRZ)) deallocate(this%psiRZ)
    if (allocated(this%LCFS)) deallocate(this%LCFS)
    if (allocated(this%limEQD)) deallocate(this%limEQD)
  end subroutine finalize_efit_class_object

  !> \brief Read content from boozer file with given filename.
  !>
  !> \param filename: input, name of the file from whciht to read the
  !>   data.
  !> \param inp_swi_: optional input, switchin version/type of the
  !>   boozer file. So far only values 8 and 9 (default) are supported.
  subroutine read_data_of_boozer_file(this, filename, inp_swi_)

    implicit none

    class (boozer_data_type) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: inp_swi_

    integer :: inp_swi
    integer :: i, j
    integer :: i_alloc
    integer :: r_un
    integer :: m_max_pert, n_max_pert, mnmax_pert
    character(len=150) :: char_dummy

    if (present(inp_swi_)) then
      inp_swi = inp_swi_
    else
      inp_swi = 9
    end if

    r_un = get_free_unit()
    open(unit=r_un,file=trim(filename),status='old',action='read')

    ! Reads lines starting with 'CC' and the first line after these,
    ! which should be the header for the global variables.
    do
      read (r_un,*) char_dummy
      if (char_dummy(1:2) /= 'CC') then
        exit
      end if
    end do
    ! Read global variables array sizes/flux/radius
    read (r_un,*) this%m0b, this%n0b, this%nsurf, this%nper, this%flux, &
      & this%a, this%R

    m_max_pert = this%m0b + 1
    n_max_pert = this%n0b + 1
    mnmax_pert = m_max_pert*n_max_pert

    ! allocate storage arrays
    allocate(this%m(this%nsurf, mnmax_pert), this%n(this%nsurf, mnmax_pert), stat = i_alloc)
    if (i_alloc /= 0) then
      stop "Allocation for the arrays containing the mode numbers " // &
            & "of the perturbation field failed!"
    end if

    allocate(this%s(this%nsurf), this%iota(this%nsurf), &
        & this%Jpol_nper(this%nsurf), this%Itor(this%nsurf), &
        & this%pprime(this%nsurf), this%sqrt_g00(this%nsurf), &
        & stat = i_alloc)
    if (i_alloc /= 0) then
      stop "Allocation for the real arrays containing " // &
           & "the radial dependent fields failed!"
    end if

    allocate(this%rmnc(this%nsurf, mnmax_pert), this%zmnc(this%nsurf, mnmax_pert),&
         this%vmnc(this%nsurf, mnmax_pert), this%bmnc(this%nsurf, mnmax_pert), stat = i_alloc)
    if (i_alloc /= 0) then
      stop "Allocation for the Fourier arrays for the perturbation field failed!"
    end if

    allocate(this%rmns(this%nsurf, mnmax_pert), this%zmns(this%nsurf, mnmax_pert),&
         this%vmns(this%nsurf, mnmax_pert), this%bmns(this%nsurf, mnmax_pert), stat = i_alloc)
    if (i_alloc /= 0) then
      stop "Allocation for the Fourier arrays for the perturbation field failed!"
    end if
    ! Might not be set in the following loop, thus initialize them here.
    this%rmns = 0
    this%zmns = 0
    this%vmns = 0
    this%bmns = 0

    ! read input arrays
    do i =1, this%nsurf
      read(r_un,*) char_dummy ! Header: variable names
      read(r_un,*) char_dummy ! Header: units
      read(r_un,*) this%s(i), this%iota(i), this%Jpol_nper(i), &
        & this%Itor(i), this%pprime(i), this%sqrt_g00(i) ! Header: values
      read(r_un,*) char_dummy ! variable names for block.

      do j=1,mnmax_pert

        ! If clause inside the loop is not efficient, but reduces code
        ! doubling, also, this part is not considered time-sensitive.
        if (inp_swi .EQ. 8) then ! NEW IPP TOKAMAK
          read(r_un,*) this%m(i,j), this%n(i,j), &
            & this%rmnc(i,j), this%zmnc(i,j), this%vmnc(i,j),&
            & this%bmnc(i,j)
        elseif (inp_swi .EQ. 9) then ! ASDEX-U (E. Strumberger)
          read(r_un,*) this%m(i,j), this%n(i,j), &
            & this%rmnc(i,j), this%rmns(i,j), this%zmnc(i,j), this%zmns(i,j),&
            & this%vmnc(i,j), this%vmns(i,j), this%bmnc(i,j), this%bmns(i,j)
        else
          print *,'FATAL: There is yet no other input type for the perturbed field defined'
          print *,'  possible values for inp_swi: 8, 9'
          stop
        end if
      end do
    end do

    close(r_un)

  end subroutine read_data_of_boozer_file

  !> \brief Write data of boozer class to file with the given filename.
  !>
  !> \param filename: input, data is written to file with this name.
  subroutine write_data_of_boozer_file(this, filename, inp_swi_)

    implicit none

    class (boozer_data_type) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: inp_swi_

    integer :: inp_swi
    integer :: i, j
    integer :: r_un
    integer :: m_max_pert, n_max_pert, mnmax_pert

    if (present(inp_swi_)) then
      inp_swi = inp_swi_
    else
      inp_swi = 9
    end if

    r_un = get_free_unit()
    open(unit=r_un,file=trim(filename),status='replace',action='write')




    write (r_un,'(a)') 'CC Boozer-coordinate data file'
    write (r_un,'(a)') 'CC Version: 01'
    write (r_un,'(a)') 'CC Author:  neo-2'
    write (r_un,'(a)') 'CC shot:'
    write (r_un,'(a)') ' m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]'
    write (r_un,format_boozer_head_param) this%m0b, this%n0b, this%nsurf, this%nper, this%flux, &
      & this%a, this%R

    m_max_pert = this%m0b + 1
    n_max_pert = this%n0b + 1
    mnmax_pert = m_max_pert*n_max_pert

    ! read input arrays
    do i =1, this%nsurf
      write(r_un,'(a)') '        s               iota           Jpol/nper' // &
                    & '          Itor            pprime         sqrt g(0,0)'
      write(r_un,'(a)') '                                          [A]' // &
                    & '           [A]             [Pa]         (dV/ds)/nper'
      write(r_un,format_boozer_flux_head) this%s(i), this%iota(i), this%Jpol_nper(i), &
        & this%Itor(i), this%pprime(i), this%sqrt_g00(i) ! Header: values
      write(r_un,'(a)') '    m    n      rmnc [m]         rmns [m]         ' // &
                    & 'zmnc [m]         zmns [m]         vmnc [ ]         ' // &
                    & 'vmns [ ]         bmnc [T]         bmns [T]'

      do j=1,mnmax_pert

        ! If clause inside the loop is not efficient, but reduces code
        ! doubling, also, this part is not considered time-sensitive.
        if (inp_swi .EQ. 8) then ! NEW IPP TOKAMAK
          write(r_un,format_boozer_output_data) this%m(i,j), this%n(i,j), &
            & this%rmnc(i,j), this%zmnc(i,j), this%vmnc(i,j),&
            & this%bmnc(i,j)
        elseif (inp_swi .EQ. 9) then ! ASDEX-U (E. Strumberger)
          write(r_un,format_boozer_output_data) this%m(i,j), this%n(i,j), &
            & this%rmnc(i,j), this%rmns(i,j), this%zmnc(i,j), this%zmns(i,j),&
            & this%vmnc(i,j), this%vmns(i,j), this%bmnc(i,j), this%bmns(i,j)
        else
          print *,'FATAL: There is yet no other input type for the perturbed field defined'
          print *,'  possible values for inp_swi: 8, 9'
          stop
        end if
      end do
    end do

    close(r_un)

  end subroutine write_data_of_boozer_file

  function get_m0b_(this)
    implicit none

    class(boozer_data_type), intent(in) :: this

    integer :: get_m0b_

    get_m0b_ = this%m0b
  end function get_m0b_

  function get_n0b_(this)
    implicit none

    class(boozer_data_type), intent(in) :: this

    integer :: get_n0b_

    get_n0b_ = this%n0b
  end function get_n0b_

  function get_nsurf_(this)
    implicit none

    class(boozer_data_type), intent(in) :: this

    integer :: get_nsurf_

    get_nsurf_ = this%nsurf
  end function get_nsurf_

  function get_nper_(this)
    implicit none

    class(boozer_data_type), intent(in) :: this

    integer :: get_nper_

    get_nper_ = this%nper
  end function get_nper_

  function get_flux_(this)
    use libneo_kinds, only : dp

    implicit none

    class(boozer_data_type), intent(in) :: this

    real(dp) :: get_flux_

    get_flux_ = this%flux
  end function get_flux_

  function get_a_(this)
    use libneo_kinds, only : dp

    implicit none

    class(boozer_data_type), intent(in) :: this

    real(dp) :: get_a_

    get_a_ = this%a
  end function get_a_

  function get_R_(this)
    use libneo_kinds, only : dp

    implicit none

    class(boozer_data_type), intent(in) :: this

    real(dp) :: get_R_

    get_R_ = this%R
  end function get_R_

end module io
