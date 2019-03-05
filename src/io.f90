!> \brief Module for handling input and output.
!>
!> This module is intended for handling the input and output of code.
!> This may include interfaces to handle input/output independent of the
!> requested/selected file format.
!> It includes routines to read in (and maybe write) specific file
!> formats, like efit and boozer. Note that the input is just read, not
!> processed.
module io
  use libneo_kinds, only : real_kind

  implicit none

  private

  public get_free_unit

  !> Storing last used file unit number. Used for getting a free one.
  integer, save :: iunit=100

  !> \brief Class representing efit data (-file).
  type efit_data_type
    private

    integer :: nwEQD, nhEQD
    integer :: n_bndyxy, nlimEQD

    real(kind=real_kind) :: psiSep, bt0, rzero

    real(kind=real_kind), dimension(:), allocatable :: fpol, pres
    real(kind=real_kind), dimension(:), allocatable :: ffprim
    real(kind=real_kind), dimension(:), allocatable :: pprime
    real(kind=real_kind), dimension(:), allocatable :: qpsi
    real(kind=real_kind), dimension(:,:), allocatable :: psiRZ
    real(kind=real_kind), dimension(:), allocatable :: LCFS, limEQD
    ! These two are for storing the grid coordinates. Regular -> only one dimension.
    real(kind=real_kind), dimension(:), allocatable :: rad, zet

    real(kind=real_kind) :: xdim,zdim,r1,zmid,rmaxis,zmaxis
    real(kind=real_kind) :: plas_cur, psiAxis

    !> This array of strings collects the part at the beginning of the
    !> first line. It is not clear if this contains some relevant
    !> information or not, so we store it to be save (it is just 60
    !> bytes per file and it is not expected that there will be more
    !> than one or two around).
    character(len=10) :: dummy(6)
  contains
    procedure :: read_data => read_data_of_efit_file
    procedure :: read_dimension => read_dimension_of_efit_file

    final :: finalize_efit_class_object
  end type efit_data_type

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
    55    print *, 'READ_EQDIM1: Early EOF in',trim(filename); STOP
    250   print *, 'READ_EQDIM1: Error reading ',trim(filename); STOP
  end subroutine read_dimension_of_efit_file

  !> \brief Read content from a efit file with given filename.
  !>
  !>
  !> \param filename: input, name of the file from which to read the
  !>   data.
  !> \param nwEQD, nhEQD: output, the dimensions read from the file.
  subroutine read_data_of_efit_file(this, filename)

    use libneo_kinds, only : real_kind

    implicit none

    class(efit_data_type), intent(inout) :: this

    character(len=*), intent(in) :: filename

    ! loop variables
    integer :: i,j

    ! File unit.
    integer :: gunit

    real (kind=real_kind) :: xdum
    integer :: idum

    gunit = get_free_unit()

    open(unit=gunit,file=trim(filename),status='old',action='read')

    ! Read in first line with sizes.
    read(gunit,fmt=2000) (this%dummy(i),i=1,6), idum, this%nwEQD, this%nhEQD

    ! Read first four lines with values:
    read(gunit,fmt=2010,end=55,err=250) this%xdim, this%zdim, this%rzero, this%r1, this%zmid
    read(gunit,fmt=2010,end=55,err=250) this%rmaxis, this%zmaxis, this%psiAxis, this%psiSep, this%bt0
    read(gunit,fmt=2010,end=55,err=250) this%plas_cur, this%psiAxis, xdum, this%rmaxis, xdum
    read(gunit,fmt=2010,end=55,err=250) this%zmaxis, xdum, this%psiSep, xdum, xdum

    ! Allocate first set of arrays. Done after reading first four lines,
    ! in case there was a problem, arrays are not already allocted.
    allocate(this%fpol(this%nwEQD), this%pres(this%nwEQD), this%ffprim(this%nwEQD))
    allocate(this%pprime(this%nwEQD), this%qpsi(this%nwEQD))
    allocate(this%psiRZ(this%nwEQD,this%nhEQD))

    read(gunit,fmt=2010,end=55,err=250) (this%fpol(i),i=1,this%nwEQD)
    read(gunit,fmt=2010,end=55,err=250) (this%pres(i),i=1,this%nwEQD)
    read(gunit,fmt=2010,end=55,err=250) (this%ffprim(i),i=1,this%nwEQD)
    read(gunit,fmt=2010,end=55,err=250) (this%pprime(i),i=1,this%nwEQD)
    read(gunit,fmt=2010,end=55,err=250) ((this%psiRZ(i,j),i=1,this%nwEQD),j=1,this%nhEQD)
    read(gunit,fmt=2010,end=55,err=250) (this%qpsi(i),i=1,this%nwEQD)

    ! Boundary Data, first read size, allocate the arrays, last read the data.
    read(gunit,*,end=55,err=250) this%n_bndyxy, this%nlimEQD
    allocate(this%LCFS(2*this%n_bndyxy))
    allocate(this%limEQD(2*this%nlimEQD))
    read(gunit,2010,end=55,err=250)(this%LCFS(i),i=1,2*this%n_bndyxy)
    read(gunit,2010,end=55,err=250)(this%limEQD(i),i=1,2*this%nlimEQD)

    close(gunit)

    call set_array_equidistant(this%nwEQD, this%xdim, this%r1, this%rad)

    call set_array_equidistant(this%nhEQD, this%zdim, this%zmid - this%zdim/2.0, this%zet)

    ! Without the return, the print statements for the error labels
    ! would lead to compiler errors.
    return

    ! Formats used for reading the data.
    2000  format(6a8,3i4)
    2010  format(5(e16.9))

    ! Error labels
    55    print *, 'Error in read_data_of_efit_file: Early EOF in ',trim(filename); STOP
    250   print *, 'Error in read_data_of_efit_file: Error reading ',trim(filename); STOP

  end subroutine read_data_of_efit_file

  !> The array will contain values from origin up to origin+width (end
  !> points included).
  !> \note It is called origin as offset is a keyword.
  subroutine set_array_equidistant(number_of_points, width, origin, array)
    use libneo_kinds, only : real_kind

    implicit none

    integer, intent(in) :: number_of_points
    real(kind=real_kind), intent(in) :: width, origin
    real(kind=real_kind), intent(out) :: array(number_of_points)

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

end module io
