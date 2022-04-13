module hdf5_tools
  !**********************************************************
  ! Compilation of useful HDF-5 wrapper functions
  ! Date:   30.11.2015
  !**********************************************************

  !**********************************************************
  ! Include hdf5 module, hdf5 lite interface and
  ! ISO_C_BINDING for long-integer support
  !**********************************************************
  use hdf5
  use h5lt
  implicit none

  !**********************************************************
  ! Definition of complex type
  !**********************************************************
  integer, parameter :: dpp = kind(1.0d0)
  integer, parameter :: dcp = dpp
  type complex_t
     real(kind=dpp) re
     real(kind=dpp) im
  end type complex_t

  !**********************************************************
  ! Records the error code of the HDF-5 functions
  !**********************************************************
  integer :: h5error

  !**********************************************************
  ! Overwrite existing datasets with h5_add?
  !**********************************************************
  logical :: h5overwrite = .false.

  !**********************************************************
  ! Wrapper functions to add content
  !**********************************************************
  interface h5_add
     module procedure h5_add_int
     module procedure h5_add_int_1
     module procedure h5_add_int_1_nobounds
     module procedure h5_add_int_2
     module procedure h5_add_int_2_nobounds
     module procedure h5_add_int_3
     !module procedure h5_add_int8_2
     module procedure h5_add_double_0
     module procedure h5_add_double_1
     module procedure h5_add_double_1_nobounds
     module procedure h5_add_double_2
     module procedure h5_add_double_3
     module procedure h5_add_double_4
     module procedure h5_add_double_5
     module procedure h5_add_complex_1
     module procedure h5_add_complex_2
     module procedure h5_add_string
     module procedure h5_add_logical
  end interface h5_add

  !**********************************************************
  ! Wrapper functions to read content
  !**********************************************************
  interface h5_get
     module procedure h5_get_int
     module procedure h5_get_int_1
     module procedure h5_get_int_2
     module procedure h5_get_int_3
     !module procedure h5_get_int8_2
     module procedure h5_get_double_0
     module procedure h5_get_double_1
     module procedure h5_get_double_2
     module procedure h5_get_double_3
     module procedure h5_get_double_4
     module procedure h5_get_double_5
     module procedure h5_get_double_4_hyperslab
     module procedure h5_get_complex_1
     module procedure h5_get_complex_2
     module procedure h5_get_string
     module procedure h5_get_logical
  end interface h5_get

  !**********************************************************
  ! Fortran supports arrays with arbitrary bounds.
  ! Therefore, all h5_add functions save these boundaries.
  !**********************************************************
  interface h5_get_bounds
     module procedure h5_get_bounds_1
     module procedure h5_get_bounds_2
  end interface h5_get_bounds

  !**********************************************************
  ! Define unlimited dimension
  !**********************************************************
  interface h5_define_unlimited
     module procedure h5_define_unlimited_array
     module procedure h5_define_unlimited_matrix
  end interface h5_define_unlimited

  !**********************************************************
  ! Wrapper functions to append content to
  ! unlimited dimensions
  !**********************************************************
  interface h5_append
     module procedure h5_append_int_0
     module procedure h5_append_double_0
     module procedure h5_append_double_4
  end interface h5_append

contains

  !**********************************************************
  ! Initialize HDF-5 Fortran interface
  !**********************************************************
  subroutine h5_init()
    !write (*,*) "Initializing HDF5 interface"
    call h5open_f(h5error)
    call h5eset_auto_f(1, h5error)
  end subroutine h5_init

  !**********************************************************
  ! Deinitialize HDF-5 Fortran interface
  !**********************************************************
  subroutine h5_deinit()
    !write (*,*) "Deinitializing HDF5 interface"
    call h5close_f(h5error)
  end subroutine h5_deinit

  !**********************************************************
  ! Check error code
  !**********************************************************
  subroutine h5_check()
    if (h5error < 0) then
       write (*,*) "HDF5 Error"
       call H5Eprint_f(h5error)
       error stop
    end if
  end subroutine h5_check

  subroutine h5_disable_exceptions()
    call h5eset_auto_f(0, h5error)
  end subroutine h5_disable_exceptions

 subroutine h5_enable_exceptions()
    call h5eset_auto_f(1, h5error)
  end subroutine h5_enable_exceptions

  !**********************************************************
  ! Create file
  !**********************************************************
  subroutine h5_create(filename, h5id, opt_fileformat_version)
    character(len=*)           :: filename
    integer(HID_T)             :: h5id
    integer, optional          :: opt_fileformat_version

    integer                    :: fileformat_version

    if (.not. present(opt_fileformat_version)) then
       fileformat_version = 1
    else
       fileformat_version = opt_fileformat_version
    end if

    !write (*,*) "Creating HDF5 file: ", trim(filename)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, h5id, h5error)
    call h5_add(h5id, 'version', fileformat_version)
    call h5_check()

  end subroutine h5_create

  !**********************************************************
  ! Close file
  !**********************************************************
  subroutine h5_close(h5id)
    integer(HID_T)             :: h5id

    call h5fclose_f(h5id, h5error)
  end subroutine h5_close

  !**********************************************************
  ! Open file if exists
  !**********************************************************
  subroutine h5_open(filename, h5id)
    character(len=*)      :: filename
    integer(HID_T)        :: h5id
    logical               :: f_exists

    !write (*,*) "Opening HDF5 file: ", trim(filename)

    inquire (file=filename, exist=f_exists)
    if (f_exists) then
       call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, h5id, h5error)
    else
       write (*,*) "HDF5 file does not exist:", filename
       error stop
    end if
  end subroutine h5_open

  !**********************************************************
  ! Open file to read-write
  !**********************************************************
  subroutine h5_open_rw(filename, h5id, opt_fileformat_version)
    character(len=*)      :: filename
    integer(HID_T)        :: h5id
    logical               :: f_exists
    integer, optional     :: opt_fileformat_version
    integer               :: fileformat_version

    !write (*,*) "Opening HDF5 file: ", filename
    inquire (file=filename, exist=f_exists)
    if (f_exists) then
       call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, h5id, h5error)
    else
       fileformat_version = 1
       if (present(opt_fileformat_version)) fileformat_version = opt_fileformat_version
       call h5_create(filename, h5id, fileformat_version)
       !write (*,*) "File does not exist!"
       !error stop
    end if
  end subroutine h5_open_rw

  !**********************************************************
  ! Define group
  !**********************************************************
  subroutine h5_define_group(h5id, grpname, h5grpid)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer(HID_T)       :: h5grpid

    call h5gcreate_f(h5id, trim(adjustl(grpname)), h5grpid, h5error)
    call h5_check()
  end subroutine h5_define_group

  !**********************************************************
  ! Open group
  !**********************************************************
  subroutine h5_open_group(h5id, grpname, h5id_grp)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer(HID_T)       :: h5id_grp

    !write (*,*) "Opening group ", grpname, "."
    call h5gopen_f(h5id, trim(adjustl(grpname)), h5id_grp, h5error)
    call h5_check()
  end subroutine h5_open_group

  !**********************************************************
  ! Close group. This is important since
  ! there is a limited number of open groups at one time
  !**********************************************************
  subroutine h5_close_group(h5id_grp)
    integer(HID_T) :: h5id_grp

    call h5gclose_f(h5id_grp, h5error)
  end subroutine h5_close_group

  !**********************************************************
  ! Get number of group members. Useful to iterate
  ! over all groups
  !**********************************************************
  subroutine h5_get_nmembers(h5id, grpname, nmembers)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname
    integer              :: nmembers

    call h5gn_members_f(h5id, grpname, nmembers, h5error)
  end subroutine h5_get_nmembers

  !**********************************************************
  ! Get information about an object. Is it a group,
  ! a variable, an attribute, ...
  !**********************************************************
  subroutine h5_get_objinfo(h5id, grpname, idx, objname, type_obj)
    integer(HID_T)       :: h5id
    character(len=*)     :: grpname, objname
    integer              :: idx
    integer              :: type_obj

    call h5gget_obj_info_idx_f(h5id, grpname , idx, objname, type_obj, &
         & h5error)
  end subroutine h5_get_objinfo

  function h5_isvalid(h5id) result(valid)
    integer(HID_T)       :: h5id
    logical              :: valid

    call h5iis_valid_f(h5id, valid, h5error)
    call h5_check()
  end function h5_isvalid

  subroutine h5_obj_exists(h5id, name_obj, exists)
    integer(HID_T)       :: h5id
    character(len=*)     :: name_obj
    logical              :: exists
    !integer              :: storage_type, nlinks, max_corder

    !call h5_disable_exceptions()
    !call h5gget_info_by_name_f(h5id, trim(name), &
    !     storage_type, nlinks, max_corder, h5error)
    !exists = h5error .eq. 0
    !call h5_enable_exceptions()
    exists = h5_exists(h5id, name_obj)
  end subroutine h5_obj_exists

  function h5_exists(h5id, name_obj) result(exists)
    integer(HID_T)       :: h5id
    character(len=*)     :: name_obj
    logical              :: exists

    call h5lexists_f(h5id, name_obj, exists, h5error)
    call h5_check()
  end function h5_exists

  subroutine h5_write_opened_obj_count(h5id)
    integer(HID_T)       :: h5id
    integer(SIZE_T)      :: obj_count

    call h5fget_obj_count_f(h5id, H5F_OBJ_ALL_F, obj_count, h5error)
    write (*,*) "Opened HDF5 objects: ", obj_count

  end subroutine h5_write_opened_obj_count

  subroutine h5_delete(h5id, name_obj)
    integer(HID_T)       :: h5id
    character(len=*)     :: name_obj
    logical              :: exists

    !write (*,*) "Trying to delete ", name
    call h5lexists_f(h5id, name_obj, exists, h5error)
    if (exists) then
       call h5ldelete_f(h5id, name_obj, h5error)
       call h5_check()
    end if
    !write (*,*) exists
  end subroutine h5_delete

  !**********************************************************
  ! Helper routine to recursively create parent groups
  ! for given dataset if they do not exist yet.
  !**********************************************************
  recursive subroutine h5_create_parent_groups(h5id, dataset)
    integer(HID_T), intent(in) :: h5id
    character(len = *), intent(in) :: dataset
    integer :: sep
    integer(HID_T) :: grp_id

    sep = index(dataset, '/')
    if (sep > 1) then
       if (h5_exists(h5id, dataset(:sep-1))) then
          call h5_open_group(h5id, dataset(:sep-1), grp_id)
       else
          call h5_define_group(h5id, dataset(:sep-1), grp_id)
       end if
       if (len_trim(dataset) > sep) call h5_create_parent_groups(grp_id, dataset(sep+1:))
       call h5_close_group(grp_id)
    end if
  end subroutine h5_create_parent_groups

  !**********************************************************
  ! Define matrix with unlimited dimensions. Used for
  ! appending data with an unknown number of elements.
  !**********************************************************
  subroutine h5_define_unlimited_matrix(h5id, dataset, datatype, dims, dsetid)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(:)  :: dims
    integer(HID_T)         :: datatype
    integer(HID_T)         :: dsetid

    integer                :: rank
    integer(SIZE_T), dimension(:), allocatable :: maxdims, startdims
    integer(HID_T)         :: dspaceid
    integer(HID_T)         :: crp_list
    integer                :: k

    rank = size(dims,1)
    allocate(maxdims(rank), startdims(rank))
    maxdims = dims
    startdims = dims

    do k = lbound(maxdims,1), ubound(maxdims,1)
       if (maxdims(k) == -1) maxdims(k) = H5S_UNLIMITED_F
       if (maxdims(k) == -1) startdims(k)  = 1
    end do
    !write (*,*) "Defining chunk: ", startdims
    call h5screate_simple_f(rank, startdims, dspaceid, h5error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, h5error)
    call h5pset_chunk_f(crp_list, rank, startdims, h5error)
    call h5dcreate_f(h5id, dataset, datatype, dspaceid, dsetid, h5error, crp_list )
    call h5sclose_f(dspaceid, h5error)
    call h5pclose_f(crp_list, h5error)

    deallocate(maxdims)
    deallocate(startdims)

  end subroutine h5_define_unlimited_matrix

  !**********************************************************
  ! Same as h5_define_unlimited_matrix
  !**********************************************************
  subroutine h5_define_unlimited_array(h5id, dataset, datatype, dsetid)
    integer(HID_T)       :: h5id
    character(len=*)     :: dataset
    integer(HID_T)       :: datatype
    integer(HID_T)       :: dsetid

    integer(HID_T)        :: dspaceid
    integer(HID_T)        :: crp_list
    integer               :: rank = 1
    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: maxdims

    maxdims = (/H5S_UNLIMITED_F/)

    call h5screate_simple_f(rank, dims, dspaceid, h5error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, h5error)
    call h5pset_chunk_f(crp_list, rank, dims, h5error)
    call h5dcreate_f(h5id, dataset, datatype, dspaceid, &
         dsetid, h5error, crp_list )
    call h5sclose_f(dspaceid, h5error)
    call h5pclose_f(crp_list, h5error)

  end subroutine h5_define_unlimited_array

  !**********************************************************
  ! Get bounds of array
  !**********************************************************
  subroutine h5_get_bounds_1(h5id, dataset, lb1, ub1)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(1)  :: lbounds, ubounds
    integer, intent(inout) :: lb1, ub1
    logical                :: attr_exists

    lbounds(1) = 0
    ubounds(1) = 0

    call h5aexists_by_name_f(h5id, dataset, 'lbounds', attr_exists, h5error)
    if (attr_exists) call h5ltget_attribute_int_f(h5id, dataset,'lbounds', lbounds(1), h5error)
    call h5aexists_by_name_f(h5id, dataset, 'ubounds', attr_exists, h5error)
    if (attr_exists) call h5ltget_attribute_int_f(h5id, dataset,'ubounds', ubounds(1), h5error)

    lb1 = lbounds(1)
    ub1 = ubounds(1)

  end subroutine h5_get_bounds_1

  !**********************************************************
  ! Get bounds of matrix
  !**********************************************************
  subroutine h5_get_bounds_2(h5id, dataset, lb1, lb2, ub1, ub2)
    integer(HID_T)         :: h5id
    character(len=*)       :: dataset
    integer, dimension(1:2):: lbounds, ubounds
    integer, intent(out)   :: lb1, lb2, ub1, ub2
    logical                :: attr_exists

    lbounds(1) = 0
    lbounds(2) = 0
    ubounds(1) = 0
    ubounds(2) = 0

    call h5aexists_by_name_f(h5id, dataset, 'lbounds', attr_exists, h5error)
    if (attr_exists) call h5ltget_attribute_int_f(h5id, dataset,'lbounds', lbounds, h5error)
    call h5aexists_by_name_f(h5id, dataset, 'ubounds', attr_exists, h5error)
    if (attr_exists) call h5ltget_attribute_int_f(h5id, dataset,'ubounds', ubounds, h5error)

    lb1 = lbounds(1)
    lb2 = lbounds(2)
    ub1 = ubounds(1)
    ub2 = ubounds(2)

    !write (*,*) "get_bounds: ", dataset, lbounds, ubounds

  end subroutine h5_get_bounds_2

  !**********************************************************
  ! Append integer scalar to unlimited dimension
  !**********************************************************
  subroutine h5_append_int_0(dsetid, value, offset)
    integer(HID_T)       :: dsetid
    integer              :: value
    integer              :: offset

    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: size
    integer(HID_T)                :: memspace
    integer                       :: rank = 1
    integer(HID_T)                :: dspaceid

    integer(HSIZE_T), dimension(1) :: offsetd
    integer(HSIZE_T), dimension(1) :: countd

    size    = (/offset/)
    offsetd = (/offset-1/)
    countd  = (/1/)
    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, countd, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)

  end subroutine h5_append_int_0

  !**********************************************************
  ! Append double scalar to unlimited dimension
  !**********************************************************
 subroutine h5_append_double_0(dsetid, value, offset)
    integer(HID_T)       :: dsetid
    double precision     :: value
    integer              :: offset

    integer(SIZE_T), dimension(1) :: dims = (/1/)
    integer(SIZE_T), dimension(1) :: size
    integer(HID_T)                :: memspace
    integer                       :: rank = 1
    integer(HID_T)                :: dspaceid

    integer(HSIZE_T), dimension(1) :: offsetd
    integer(HSIZE_T), dimension(1) :: countd

    size    = (/offset/)
    offsetd = (/offset-1/)
    countd  = (/1/)

    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, countd, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)

  end subroutine h5_append_double_0

  !**********************************************************
  ! Append double matrix to unlimited dimension
  !**********************************************************
  subroutine h5_append_double_4(dsetid, value, offset)
    integer(HID_T)                     :: dsetid
    double precision, dimension(:,:,:) :: value
    integer                            :: offset

    integer(SIZE_T), dimension(4)  :: dims
    integer(SIZE_T), dimension(4)  :: size
    integer(HID_T)                 :: memspace
    integer                        :: rank = 4
    integer(HID_T)                 :: dspaceid

    integer(HSIZE_T), dimension(4) :: offsetd

    size    = (/shape(value), offset/)
    dims    = (/shape(value), 1/)
    offsetd = (/0, 0, 0, offset-1/)
    !offsetd = offset
    !countd  = shape(value)

    !write (*,*) "Size:", size
    !write (*,*) "Dims:", dims
    !write (*,*) "Offset:", offsetd
    call h5dset_extent_f(dsetid, size, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, dims, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)

  end subroutine h5_append_double_4

  !**********************************************************
  ! Add logical
  !**********************************************************
  subroutine h5_add_logical(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    logical                        :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit
    integer(HSIZE_T)               :: dims(1) = (/0/)
    integer                        :: internalvalue

    if (h5overwrite) call h5_delete(h5id, dataset)
    internalvalue = 0
    if (value) internalvalue = 1
    call h5ltmake_dataset_int_f(h5id, dataset, 0,dims, (/internalvalue/), h5error)
    call h5_set_optional_attributes(h5id, dataset, comment, unit)
    call h5_check()

  end subroutine h5_add_logical

  !**********************************************************
  ! Add integer scalar
  !**********************************************************
  subroutine h5_add_int(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    integer                        :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit
    integer(HSIZE_T)               :: dims(1) = (/0/)

    if (h5overwrite) call h5_delete(h5id, dataset)
    call h5ltmake_dataset_int_f(h5id, dataset, 0,dims, (/value/), h5error)
    call h5_set_optional_attributes(h5id, dataset, comment, unit)
    call h5_check()

  end subroutine h5_add_int

  !**********************************************************
  ! Add integer array
  !**********************************************************
  subroutine h5_add_int_1(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:)             :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes(h5id, dataset, comment, unit)

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_int_1

  subroutine h5_add_int_1_nobounds(h5id, dataset, value, comment, unit, default)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, allocatable, dimension(:)    :: value
    !integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer, optional        :: default
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    if (h5overwrite) call h5_delete(h5id, dataset)
    if ( allocated(value) ) then
       allocate(dims(rank))
       dims = ubound(value) - lbound(value) + 1
       size = rank
       call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbound(value), size, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubound(value), size, h5error)

       call h5_set_optional_attributes(h5id, dataset, comment, unit)

       deallocate(dims)
    else
       ! value is not allocated
       if (present(default) ) then
          call h5_add(h5id, dataset, default)
       else
          call h5_add(h5id, dataset, 0)
       end if
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', 'value not allocated', h5error)
    end if
    call h5_check()
  end subroutine h5_add_int_1_nobounds

  !**********************************************************
  ! Add integer matrix
  !**********************************************************
  subroutine h5_add_int_2(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:,:)           :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 2

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes(h5id, dataset, comment, unit)

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_int_2

  subroutine h5_add_int_2_nobounds(h5id, dataset, value, comment, unit, default)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, allocatable, dimension(:,:)    :: value
    !integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer, optional        :: default
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 2

    if (h5overwrite) call h5_delete(h5id, dataset)
    if ( allocated(value) ) then
       allocate(dims(rank))
       dims = ubound(value) - lbound(value) + 1
       size = rank
       call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbound(value), size, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubound(value), size, h5error)

       call h5_set_optional_attributes(h5id, dataset, comment, unit)

       deallocate(dims)
    else
       ! value is not allocated
       if (present(default) ) then
          call h5_add(h5id, dataset, default)
       else
          call h5_add(h5id, dataset, 0)
       end if
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', 'value not allocated', h5error)
    end if
    call h5_check()
  end subroutine h5_add_int_2_nobounds

  !**********************************************************
  ! Add integer array
  !**********************************************************
  subroutine h5_add_int_3(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:,:,:)         :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 3

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_int_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_int_3

  !**********************************************************
  ! Add long integer matrix. This function makes use of the
  ! HDF-5 Fortran 2003 interface, since the default HDF-5
  ! functions to not support integer(kind=8).
  ! This is documentated at
  ! https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=758694
  !**********************************************************
!!$  subroutine h5_add_int8_2(h5id, dataset, value, lbounds, ubounds, comment, unit)
!!$    integer(HID_T)                           :: h5id
!!$    character(len=*)                         :: dataset
!!$    integer(kind=8), dimension(:,:),target   :: value
!!$    integer, dimension(:)                    :: lbounds, ubounds
!!$    character(len=*), optional               :: comment
!!$    character(len=*), optional               :: unit
!!$    integer(HSIZE_T), dimension(:), allocatable    :: dims
!!$    integer(SIZE_T)                          :: size
!!$    integer                                  :: rank = 2
!!$    integer(HID_T)                           :: dspace_id, dset_id
!!$    integer(kind=8), dimension(:,:), pointer :: test
!!$    type(C_PTR)                              :: f_ptr
!!$    integer(HID_T)                           :: h5_kind_type_i
!!$
!!$    allocate(dims(rank))
!!$    dims = ubounds - lbounds + 1
!!$    size = rank
!!$    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)
!!$
!!$    call h5screate_simple_f(rank, dims, dspace_id, h5error)
!!$    call h5dcreate_f(h5id, dataset, h5_kind_type_i, dspace_id, dset_id, h5error)
!!$
!!$    test => value
!!$    f_ptr = c_loc(test(1,1))
!!$    call h5dwrite_f(dset_id, h5_kind_type_i, f_ptr, h5error)
!!$
!!$    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
!!$    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)
!!$
!!$    call h5_set_optional_attributes(h5id, dataset, comment, unit)
!!$
!!$    call h5dclose_f(dset_id, h5error)
!!$    call h5sclose_f(dspace_id, h5error)
!!$
!!$    deallocate(dims)
!!$
!!$    call h5_check()
!!$  end subroutine h5_add_int8_2

  !**********************************************************
  ! Get long integer matrix
  !**********************************************************
!!$  subroutine h5_get_int8_2(h5id, dataset, value)
!!$    integer(HID_T)                           :: h5id
!!$    character(len=*)                         :: dataset
!!$    integer(kind=8), dimension(:,:), target  :: value
!!$    integer                                  :: lb1, lb2, ub1, ub2
!!$    integer(HSIZE_T), dimension(2)           :: dims
!!$    integer(HID_T)                           :: dspace_id, dset_id
!!$    integer(kind=8), dimension(:,:), pointer :: test
!!$    integer(HID_T)                           :: h5_kind_type_i
!!$    type(C_PTR)                              :: f_ptr
!!$
!!$    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)
!!$
!!$    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
!!$    dims = (/ub1-lb1+1, ub2-lb2+1/)
!!$
!!$    call h5dopen_f(h5id, dataset, dset_id, h5error)
!!$    call h5dget_space_f(dset_id, dspace_id, h5error)
!!$    test => value
!!$    f_ptr = c_loc(test(1,1))
!!$
!!$    call h5dread_f(dset_id, h5_kind_type_i, f_ptr, h5error)
!!$
!!$    call h5dclose_f(dset_id, h5error)
!!$    call h5sclose_f(dspace_id, h5error)
!!$
!!$    call h5_check()
!!$
!!$  end subroutine h5_get_int8_2

  !**********************************************************
  ! Get integer scalar
  !**********************************************************
  subroutine h5_get_int(h5id, dataset, value)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    integer, intent(out)           :: value
    integer, dimension(1)          :: buf
    integer(HSIZE_T)               :: dims(1) = (/0/)

    call h5ltread_dataset_int_f_1(h5id, dataset, buf, dims, h5error)
    value = buf(1)

    call h5_check()

  end subroutine h5_get_int

  !**********************************************************
  ! Get integer array
  !**********************************************************
  subroutine h5_get_int_1(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    integer, dimension(:)             :: value
    integer                           :: lb1, ub1
    integer(HSIZE_T), dimension(1)    :: dims

    call h5_get_bounds(h5id, dataset, lb1, ub1)
    dims = (/ub1 - lb1 + 1/)
    !write (*,*) dims, 'UB', ub1, lb1
    call h5ltread_dataset_int_f(h5id, dataset, value, dims, h5error)

    call h5_check()

  end subroutine h5_get_int_1

  !**********************************************************
  ! Add double scalar
  !**********************************************************
  subroutine h5_add_double_0(h5id, dataset, value, comment, unit, accuracy)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    double precision               :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit
    double precision, optional     :: accuracy

    integer(HSIZE_T)               :: dims(1) = (/0/)

    if (h5overwrite) call h5_delete(h5id, dataset)
    call h5ltmake_dataset_double_f(h5id, dataset, 0, dims, (/value/), h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    call h5_check()

  end subroutine h5_add_double_0

  !**********************************************************
  ! Add double array
  !**********************************************************
  subroutine h5_add_double_1(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:)    :: value
    integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    double precision, optional        :: accuracy

    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)
    call h5_check()

  end subroutine h5_add_double_1

  subroutine h5_add_double_1_nobounds(h5id, dataset, value, comment, unit, default, accuracy)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, allocatable, dimension(:)    :: value
    !integer, dimension(:)             :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit
    double precision, optional        :: default
    double precision, optional        :: accuracy

    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                   :: size
    integer                           :: rank = 1

    if (h5overwrite) call h5_delete(h5id, dataset)
    if ( allocated(value) ) then
       allocate(dims(rank))
       dims = ubound(value) - lbound(value) + 1
       size = rank
       call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbound(value), size, h5error)
       call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubound(value), size, h5error)

       call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

       deallocate(dims)
    else
       ! value is not allocated
       if (present(default) ) then
          call h5_add(h5id, dataset, default)
       else
          call h5_add(h5id, dataset, 0.0d0)
       end if
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', 'value not allocated', h5error)
    end if
    call h5_check()
  end subroutine h5_add_double_1_nobounds

  !**********************************************************
  ! Get integer matrix
  !**********************************************************
  subroutine h5_get_int_2(h5id, dataset, value)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    integer, dimension(:,:)               :: value
    integer(HSIZE_T), dimension(2)        :: dims

    dims = shape(value)
    call h5ltread_dataset_int_f(h5id, dataset, value, dims, h5error)

    !PRINT *,dims(1),dims(2),dims(3)

    call h5_check()
  end subroutine h5_get_int_2

  !**********************************************************
  ! Get integer array
  !**********************************************************
  subroutine h5_get_int_3(h5id, dataset, value)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    integer, dimension(:,:,:)             :: value
    integer(HSIZE_T), dimension(3)        :: dims

    dims = shape(value)
    call h5ltread_dataset_int_f(h5id, dataset, value, dims, h5error)
    call h5_check()
  end subroutine h5_get_int_3

  !**********************************************************
  ! Get double scalar
  !**********************************************************
  subroutine h5_get_double_0(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, intent(out)     :: value
    double precision, dimension(1)    :: buf
    integer(HSIZE_T), dimension(1)    :: dims = (/0/)

    call h5ltread_dataset_double_f(h5id, dataset, buf, dims, h5error)
    value = buf(1)

    call h5_check()

  end subroutine h5_get_double_0

  !**********************************************************
  ! Get double array
  !**********************************************************
  subroutine h5_get_double_1(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:)    :: value
    integer                           :: lb1, ub1
    integer(HSIZE_T), dimension(1)    :: dims

    call h5_get_bounds(h5id, dataset, lb1, ub1)
    !if (.not. unlimited) then
    !   dims = (/ub1 - lb1 + 1/)
    !else
       dims = shape(value)
       !write (*,*) "Unlimited dimension ", dims
    !end if
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)
    call h5_check()
  end subroutine h5_get_double_1

  !**********************************************************
  ! Get double matrix
  !**********************************************************
   subroutine h5_get_double_2(h5id, dataset, value)
    integer(HID_T)                    :: h5id
    character(len=*)                  :: dataset
    double precision, dimension(:,:)  :: value
    integer                           :: lb1, lb2, ub1, ub2
    integer(HSIZE_T), dimension(2)    :: dims

    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
    dims = (/ub1-lb1+1, ub2-lb2+1/)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    call h5_check()
  end subroutine h5_get_double_2

  !**********************************************************
  ! Get double 3-dim-matrix
  !**********************************************************
  subroutine h5_get_double_3(h5id, dataset, value)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    double precision, dimension(:,:,:)    :: value
    integer(HSIZE_T), dimension(3)        :: dims

    dims = shape(value)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    !PRINT *,dims(1),dims(2),dims(3)

    call h5_check()
  end subroutine h5_get_double_3

  !**********************************************************
  ! Get double 4-dim-matrix
  !**********************************************************
  subroutine h5_get_double_4(h5id, dataset, value)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    double precision, dimension(:,:,:,:)  :: value
    integer(HSIZE_T), dimension(4)        :: dims

    dims = shape(value)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    call h5_check()
  end subroutine h5_get_double_4

  !**********************************************************
  ! Get double 5-dim-matrix
  !**********************************************************
  subroutine h5_get_double_5(h5id, dataset, value)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    double precision, dimension(:,:,:,:,:):: value
    integer(HSIZE_T), dimension(5)        :: dims

    dims = shape(value)
    call h5ltread_dataset_double_f(h5id, dataset, value, dims, h5error)

    call h5_check()
  end subroutine h5_get_double_5

  !**********************************************************
  ! Get double 4-dim-matrix
  !**********************************************************
  subroutine h5_get_double_4_hyperslab(h5id, dataset, value, offset_par, count_par)
    integer(HID_T)                        :: h5id
    character(len=*)                      :: dataset
    double precision, dimension(:,:,:,:)  :: value
    integer,dimension(4)                  :: offset_par, count_par
    integer(HSIZE_T), dimension(4)        :: offset, count, startout
    integer(HSIZE_T), dimension(4)        :: dims
    integer(HID_T)                        :: dset_id, dataspace, memspace
    integer                               :: memrank = 4

    offset = offset_par
    count = count_par
    dims = shape(value)

    call h5dopen_f(h5id, dataset, dset_id, h5error)
    call h5dget_space_f(dset_id, dataspace, h5error)

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
         offset, count, h5error)

    call h5screate_simple_f(memrank, dims, memspace, h5error)

    startout = (/0,0,0,0/)
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
         startout, count, h5error)

    call H5dread_f(dset_id, H5T_NATIVE_DOUBLE, value, dims,h5error, &
         memspace, dataspace)
    call h5_check()

    call h5sclose_f(dataspace, h5error)
    call h5sclose_f(memspace, h5error)
    call h5dclose_f(dset_id, h5error)

    !write (*,*) "Read", value
    !write (*,*) "Shape", shape(value)
    !write (*,*) "Dimns", offset, count, startout, dims
    call h5_check()
  end subroutine h5_get_double_4_hyperslab

  !**********************************************************
  ! Get 1-dim complex double matrix
  !**********************************************************
  subroutine h5_get_complex_1(h5id, dataset, val)
    integer(HID_T), intent(in)                     :: h5id
    character(len=*), intent(in)                   :: dataset
    complex(kind=dcp), dimension(:), intent(inout) :: val
#if defined(__INTEL_COMPILER)
    real, allocatable, dimension(:)                :: temp
#endif
    integer                                        :: lb1, ub1
    integer, parameter                             :: rank = 1
    integer(HSIZE_T), dimension(rank)              :: dims
    integer(SIZE_T)                                :: re_size, im_size, t_size
    integer(SIZE_T)                                :: offset
    integer(HID_T)                                 :: type_id
    integer(HID_T)                                 :: dset_id, dt_re_id, dt_im_id

    !**********************************************************
    ! Get sizes
    !**********************************************************
    call h5tget_size_f(H5T_NATIVE_DOUBLE, re_size, h5error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE, im_size, h5error)
    t_size = re_size + im_size

    !**********************************************************
    ! Create compound type
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, t_size, type_id, h5error)
    offset = 0
    call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, h5error)
    offset = offset + re_size
    call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Create sub datasets
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, re_size, dt_re_id, h5error)
    offset = 0
    call h5tinsert_f(dt_re_id, "real", offset, H5T_NATIVE_DOUBLE, h5error)

    call h5tcreate_f(H5T_COMPOUND_F, im_size, dt_im_id, h5error)
    offset = 0
    call h5tinsert_f(dt_im_id, "imag", offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Get dimension of value to be retrieved
    !**********************************************************
    call h5_get_bounds(h5id, dataset, lb1, ub1)
    dims = (/ ub1 - lb1 + 1 /)

    !**********************************************************
    ! Read data
    !**********************************************************
    call h5dopen_f(h5id, dataset, dset_id, h5error)

#if defined(__INTEL_COMPILER)
    if (allocated(temp)) deallocate(temp)
    allocate(temp(lbound(val,1):ubound(val,1)))
    call h5dread_f(dset_id, dt_re_id, temp, dims, h5error)
    val%re = temp
    call h5dread_f(dset_id, dt_im_id, temp, dims, h5error)
    val%im = temp
    deallocate(temp)
#else
    call h5dread_f(dset_id, dt_re_id, val%re, dims, h5error)
    call h5dread_f(dset_id, dt_im_id, val%im, dims, h5error)
#endif

    call h5dclose_f(dset_id, h5error)
    call h5tclose_f(type_id, h5error)
    call h5tclose_f(dt_re_id, h5error)
    call h5tclose_f(dt_im_id, h5error)

    call h5_check()
  end subroutine h5_get_complex_1

  !**********************************************************
  ! Get 2-dim complex double matrix
  !**********************************************************
  subroutine h5_get_complex_2(h5id, dataset, val)
    integer(HID_T), intent(in)                       :: h5id
    character(len=*), intent(in)                     :: dataset
    complex(kind=dcp), dimension(:,:), intent(inout) :: val
#if defined(__INTEL_COMPILER)
    real, allocatable, dimension(:,:) :: temp
#endif
    integer                                          :: lb1, ub1, lb2, ub2
    integer, parameter                               :: rank = 2
    integer(HSIZE_T), dimension(rank)                :: dims
    integer(SIZE_T)                                  :: re_size, im_size, t_size
    integer(SIZE_T)                                  :: offset
    integer(HID_T)                                   :: type_id
    integer(HID_T)                                   :: dset_id, dt_re_id, dt_im_id

    !**********************************************************
    ! Get sizes
    !**********************************************************
    call h5tget_size_f(H5T_NATIVE_DOUBLE, re_size, h5error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE, im_size, h5error)
    t_size = re_size + im_size

    !**********************************************************
    ! Create compound type
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, t_size, type_id, h5error)
    offset = 0
    call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, h5error)
    offset = offset + re_size
    call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Create sub datasets
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, re_size, dt_re_id, h5error)
    offset = 0
    call h5tinsert_f(dt_re_id, "real", offset, H5T_NATIVE_DOUBLE, h5error)

    call h5tcreate_f(H5T_COMPOUND_F, im_size, dt_im_id, h5error)
    offset = 0
    call h5tinsert_f(dt_im_id, "imag", offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Get dimension of value to be retrieved
    !**********************************************************
    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
    dims = (/ ub1 - lb1 + 1, ub2 - lb2 + 1 /)

    !**********************************************************
    ! Read data
    !**********************************************************
    call h5dopen_f(h5id, dataset, dset_id, h5error)

#if defined(__INTEL_COMPILER)
    if (allocated(temp)) deallocate(temp)
    allocate(temp(lbound(val, 1):ubound(val, 1), lbound(val, 2):ubound(val, 2)))
    call h5dread_f(dset_id, dt_re_id, temp, dims, h5error)
    val%re = temp
    call h5dread_f(dset_id, dt_im_id, temp, dims, h5error)
    val%im = temp
    deallocate(temp)
#else
    call h5dread_f(dset_id, dt_re_id, val%re, dims, h5error)
    call h5dread_f(dset_id, dt_im_id, val%im, dims, h5error)
#endif

    call h5dclose_f(dset_id, h5error)
    call h5tclose_f(type_id, h5error)
    call h5tclose_f(dt_re_id, h5error)
    call h5tclose_f(dt_im_id, h5error)

    call h5_check()
  end subroutine h5_get_complex_2

  !**********************************************************
  ! Get logical scalar
  !**********************************************************
  subroutine h5_get_logical(h5id, dataset, value)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    logical, intent(out)           :: value
    integer, dimension(1)          :: buf
    integer(HSIZE_T)               :: dims(1) = (/0/)

    call h5ltread_dataset_int_f_1(h5id, dataset, buf, dims, h5error)

    value = .false.
    if (buf(1) .eq. 1) then
       value = .true.
    end if
    call h5_check()

  end subroutine h5_get_logical

  !**********************************************************
  ! Add double matrix
  !**********************************************************
  subroutine h5_add_double_2(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    double precision, dimension(:,:)            :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 2

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5_check()
  end subroutine h5_add_double_2

  !**********************************************************
  ! Add 3-dim double matrix
  !**********************************************************
  subroutine h5_add_double_3(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    double precision, dimension(:,:,:)          :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 3

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5_check()
  end subroutine h5_add_double_3

  !**********************************************************
  ! Add 4-dim double matrix
  !**********************************************************
  subroutine h5_add_double_4(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    double precision, dimension(:,:,:,:)        :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 4

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5_check()
  end subroutine h5_add_double_4

  !**********************************************************
  ! Add 5-dim double matrix
  !**********************************************************
  subroutine h5_add_double_5(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    double precision, dimension(:,:,:,:,:)      :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 5

    if (h5overwrite) call h5_delete(h5id, dataset)
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    call h5ltmake_dataset_double_f(h5id, dataset, rank, dims, value, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5_check()
  end subroutine h5_add_double_5

  !**********************************************************
  ! Add 1-dim complex double matrix
  !**********************************************************
  subroutine h5_add_complex_1(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    complex(kind=dcp), dimension(:)             :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 1
    integer(SIZE_T)                             :: re_size, im_size, t_size
    integer(SIZE_T)                             :: offset
    integer(HID_T)                              :: type_id
    integer(HID_T)                              :: dspace_id, dset_id, dt_re_id, dt_im_id

    if (h5overwrite) call h5_delete(h5id, dataset)
    !**********************************************************
    ! Get sizes
    !**********************************************************
    call h5tget_size_f(H5T_NATIVE_DOUBLE, re_size, h5error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE, im_size, h5error)
    t_size = re_size + im_size

    !**********************************************************
    ! Create compound type
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, t_size, type_id, h5error)
    offset = 0
    call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, h5error)
    offset = offset + re_size
    call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Get dimension of value to be stored
    !**********************************************************
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank

    !**********************************************************
    ! Create dataset
    !**********************************************************
    call h5screate_simple_f(rank, dims, dspace_id, h5error)
    call h5dcreate_f(h5id, dataset, type_id, dspace_id, dset_id, h5error)

    !**********************************************************
    ! Create sub datasets
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, re_size, dt_re_id, h5error)
    offset = 0
    call h5tinsert_f(dt_re_id, "real", offset, H5T_NATIVE_DOUBLE, h5error)

    call h5tcreate_f(H5T_COMPOUND_F, im_size, dt_im_id, h5error)
    offset = 0
    call h5tinsert_f(dt_im_id, "imag", offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Write data
    !**********************************************************
    call h5dwrite_f(dset_id, dt_re_id, real(value),  dims, h5error)
    call h5dwrite_f(dset_id, dt_im_id, aimag(value), dims, h5error)

    !**********************************************************
    ! Set additional attributes
    !**********************************************************
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)
    call h5tclose_f(type_id, h5error)
    call h5tclose_f(dt_re_id, h5error)
    call h5tclose_f(dt_im_id, h5error)

    call h5_check()
  end subroutine h5_add_complex_1

  !**********************************************************
  ! Add 2-dim complex double matrix
  !**********************************************************
  subroutine h5_add_complex_2(h5id, dataset, value, lbounds, ubounds, comment, unit, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    complex(kind=dcp), dimension(:,:)           :: value
    integer, dimension(:)                       :: lbounds, ubounds
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit
    double precision, optional                  :: accuracy

    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(SIZE_T)                             :: size
    integer                                     :: rank = 2
    integer(SIZE_T)                             :: re_size, im_size, t_size
    integer(SIZE_T)                             :: offset
    integer(HID_T)                             :: type_id
    integer(HID_T)                              :: dspace_id, dset_id, dt_re_id, dt_im_id

    if (h5overwrite) call h5_delete(h5id, dataset)
    !**********************************************************
    ! Get sizes
    !**********************************************************
    call h5tget_size_f(H5T_NATIVE_DOUBLE, re_size, h5error)
    call h5tget_size_f(H5T_NATIVE_DOUBLE, im_size, h5error)
    t_size = re_size + im_size

    !**********************************************************
    ! Create compound type
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, t_size, type_id, h5error)
    offset = 0
    call h5tinsert_f(type_id, 'real', offset, H5T_NATIVE_DOUBLE, h5error)
    offset = offset + re_size
    call h5tinsert_f(type_id, 'imag', offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Get dimension of value to be stored
    !**********************************************************
    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank

    !**********************************************************
    ! Create dataset
    !**********************************************************
    call h5screate_simple_f(rank, dims, dspace_id, h5error)
    call h5dcreate_f(h5id, dataset, type_id, dspace_id, dset_id, h5error)

    !**********************************************************
    ! Create sub datasets
    !**********************************************************
    call h5tcreate_f(H5T_COMPOUND_F, re_size, dt_re_id, h5error)
    offset = 0
    call h5tinsert_f(dt_re_id, "real", offset, H5T_NATIVE_DOUBLE, h5error)

    call h5tcreate_f(H5T_COMPOUND_F, im_size, dt_im_id, h5error)
    offset = 0
    call h5tinsert_f(dt_im_id, "imag", offset, H5T_NATIVE_DOUBLE, h5error)

    !**********************************************************
    ! Write data
    !**********************************************************
    call h5dwrite_f(dset_id, dt_re_id, real(value),  dims, h5error)
    call h5dwrite_f(dset_id, dt_im_id, aimag(value), dims, h5error)

    !**********************************************************
    ! Set additional attributes
    !**********************************************************
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    call h5_set_optional_attributes_float(h5id, dataset, comment, unit, accuracy)

    deallocate(dims)

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)
    call h5tclose_f(type_id, h5error)
    call h5tclose_f(dt_re_id, h5error)
    call h5tclose_f(dt_im_id, h5error)

    call h5_check()
  end subroutine h5_add_complex_2

  !**********************************************************
  ! Add string
  !**********************************************************
  subroutine h5_add_string(h5id, dataset, value, comment, unit)
    integer(HID_T)                 :: h5id
    character(len=*)               :: dataset
    character(len=*)               :: value
    character(len=*), optional     :: comment
    character(len=*), optional     :: unit

    if (h5overwrite) call h5_delete(h5id, dataset)
    call h5ltmake_dataset_string_f(h5id, dataset, value, h5error)

    call h5_set_optional_attributes(h5id, dataset, comment, unit)

    call h5_check()
  end subroutine h5_add_string

  !**********************************************************
  ! Get string
  !**********************************************************
  subroutine h5_get_string(h5id, dataset, value)
    integer(HID_T), intent(in)     :: h5id
    character(len=*), intent(in)   :: dataset
    character(len=*), intent(out)  :: value
    integer                        :: k

    call h5ltread_dataset_string_f(h5id, dataset, value, h5error)
    k = index(value, char(0))
    if (k > 0) value(k:) = ' '  ! pad with blanks starting from null byte
    call h5_check()
  end subroutine h5_get_string

  !**********************************************************
  ! Set optional attributes for floating point type fields.
  !
  ! Set the optional attributes that are common to most fields, e.g.
  ! comment, plus some that are specific for floating point type (e.g.
  ! double and complex) fields, e.g. accuracy.
  !**********************************************************
  subroutine h5_set_optional_attributes_float(h5id, dataset, comment, unit_, accuracy)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit_
    double precision, optional    :: accuracy

    integer(SIZE_T), parameter :: size_ = 1

    call h5_set_optional_attributes(h5id, dataset, comment, unit_)

    if (present(accuracy)) then
      call h5ltset_attribute_double_f(h5id, dataset, 'accuracy', (/ accuracy /), size_, h5error)
    end if

  end subroutine h5_set_optional_attributes_float

  !**********************************************************
  ! Set optional attributes common to (most) fields.
  !
  ! Set optional attributes that are common to (most) fields, e.g.
  ! comment.
  !**********************************************************
  subroutine h5_set_optional_attributes(h5id, dataset, comment, unit_)
    integer(HID_T)                              :: h5id
    character(len=*)                            :: dataset
    character(len=*), optional                  :: comment
    character(len=*), optional                  :: unit_

    if (present(comment)) then
      call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit_)) then
      call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit_, h5error)
    end if

  end subroutine h5_set_optional_attributes

end module hdf5_tools
