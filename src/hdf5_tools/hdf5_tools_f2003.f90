module hdf5_tools_f2003
  !**********************************************************
  ! Compilation of useful HDF-5 wrapper functions
  ! These functions require HDF5 libraries compiled with
  ! Fortran 2003 support.
  ! Date:   06.08.2015
  !**********************************************************

  !**********************************************************
  ! Include hdf5 module, hdf5 lite interface and
  ! ISO_C_BINDING for long-integer support
  !**********************************************************
  use hdf5_tools
  use ISO_C_BINDING

  implicit none
  
  !**********************************************************
  ! Wrapper functions to read content
  !**********************************************************
  interface h5_get
     module procedure h5_get_int8_1
     module procedure h5_get_int8_2
  end interface h5_get

  !**********************************************************
  ! Wrapper functions to add content
  !**********************************************************
  interface h5_add
     module procedure h5_add_int8_1
     module procedure h5_add_int8_2
  end interface h5_add

contains

  subroutine h5_copy(h5id_src, srcname, h5id_dest, destname)
    integer(HID_T)   :: h5id_src, h5id_dest
    character(len=*) :: srcname, destname

    call h5ocopy_f(h5id_src, srcname, h5id_dest, destname, h5error)
  end subroutine h5_copy

  !**********************************************************
  ! Add long integer matrix. This function makes use of the
  ! HDF-5 Fortran 2003 interface, since the default HDF-5
  ! functions to not support integer(kind=8).
  ! This is documentated at
  ! https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=758694
  !**********************************************************
  subroutine h5_add_int8_2(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                           :: h5id
    character(len=*)                         :: dataset
    integer(kind=8), dimension(:,:),target   :: value
    integer, dimension(:)                    :: lbounds, ubounds
    character(len=*), optional               :: comment
    character(len=*), optional               :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                          :: size
    integer                                  :: rank = 2
    integer(HID_T)                           :: dspace_id, dset_id
    integer(kind=8), dimension(:,:), pointer :: test
    type(C_PTR)                              :: f_ptr
    integer(HID_T)                           :: h5_kind_type_i

    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)

    call h5screate_simple_f(rank, dims, dspace_id, h5error)
    call h5dcreate_f(h5id, dataset, h5_kind_type_i, dspace_id, dset_id, h5error)

    test => value
    f_ptr = c_loc(test(1,1))
    call h5dwrite_f(dset_id, h5_kind_type_i, f_ptr, h5error)
  
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)   

    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)
    
    deallocate(dims)

    call h5_check()
  end subroutine h5_add_int8_2

  subroutine h5_add_int8_1(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T)                           :: h5id
    character(len=*)                         :: dataset
    integer(kind=8), dimension(:),target     :: value
    integer, dimension(:)                    :: lbounds, ubounds
    character(len=*), optional               :: comment
    character(len=*), optional               :: unit
    integer(HSIZE_T), dimension(:), allocatable    :: dims
    integer(SIZE_T)                          :: size
    integer                                  :: rank = 1
    integer(HID_T)                           :: dspace_id, dset_id
    integer(kind=8), dimension(:), pointer   :: test
    type(C_PTR)                              :: f_ptr
    integer(HID_T)                           :: h5_kind_type_i

    allocate(dims(rank))
    dims = ubounds - lbounds + 1
    size = rank
    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)

    call h5screate_simple_f(rank, dims, dspace_id, h5error)
    call h5dcreate_f(h5id, dataset, h5_kind_type_i, dspace_id, dset_id, h5error)

    test => value
    f_ptr = c_loc(test(1))
    call h5dwrite_f(dset_id, h5_kind_type_i, f_ptr, h5error)
  
    call h5ltset_attribute_int_f(h5id, dataset, 'lbounds', lbounds, size, h5error)
    call h5ltset_attribute_int_f(h5id, dataset, 'ubounds', ubounds, size, h5error)

    if (present(comment)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'comment', comment, h5error)
    end if
    if (present(unit)) then
       call h5ltset_attribute_string_f(h5id, dataset, 'unit', unit, h5error)
    end if

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)

    deallocate(dims)

    call h5_check()
  end subroutine h5_add_int8_1

  !**********************************************************
  ! Get long integer matrix
  !**********************************************************
  subroutine h5_get_int8_2(h5id, dataset, value)
    integer(HID_T)                           :: h5id
    character(len=*)                         :: dataset
    integer(kind=8), dimension(:,:), target  :: value
    integer                                  :: lb1, lb2, ub1, ub2
    integer(HSIZE_T), dimension(2)           :: dims
    integer(HID_T)                           :: dspace_id, dset_id
    integer(kind=8), dimension(:,:), pointer :: test
    integer(HID_T)                           :: h5_kind_type_i 
    type(C_PTR)                              :: f_ptr

    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)
    
    call h5_get_bounds(h5id, dataset, lb1, lb2, ub1, ub2)
    dims = (/ub1-lb1+1, ub2-lb2+1/)

    call h5dopen_f(h5id, dataset, dset_id, h5error)
    call h5dget_space_f(dset_id, dspace_id, h5error)
    test => value
    f_ptr = c_loc(test(1,1))

    call h5dread_f(dset_id, h5_kind_type_i, f_ptr, h5error)

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)

    call h5_check()

  end subroutine h5_get_int8_2

  subroutine h5_get_int8_1(h5id, dataset, value)
    integer(HID_T)                           :: h5id
    character(len=*)                         :: dataset
    integer(kind=8), dimension(:), target    :: value
    integer                                  :: lb1, ub1
    integer(HSIZE_T), dimension(1)           :: dims
    integer(HID_T)                           :: dspace_id, dset_id
    integer(kind=8), dimension(:), pointer   :: test
    integer(HID_T)                           :: h5_kind_type_i
    type(C_PTR)                              :: f_ptr

    h5_kind_type_i = h5kind_to_type(8,H5_INTEGER_KIND)
    
    call h5_get_bounds(h5id, dataset, lb1, ub1)
    dims = (/ub1-lb1+1/)

    call h5dopen_f(h5id, dataset, dset_id, h5error)
    call h5dget_space_f(dset_id, dspace_id, h5error)
    test => value
    f_ptr = c_loc(test(1))

    call h5dread_f(dset_id, h5_kind_type_i, f_ptr, h5error)

    call h5dclose_f(dset_id, h5error)
    call h5sclose_f(dspace_id, h5error)

    call h5_check()

  end subroutine h5_get_int8_1

end module hdf5_tools_f2003
