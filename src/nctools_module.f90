module nctools_module

  use netcdf
  implicit none

  !> \brief Interface for subroutines to define variables in nc-file.
  !>
  !> This interface is for defining a variable with given name. The
  !> different variants are for different types/numbers of dimensions.
  !> With the exception of multidim, the argumentlist is for all the
  !> same. First is the id of the netcdf file in which to define the
  !> variable. Second is the name of the variable to define. Third is
  !> the field. It is only used to determine the size, if necessary. The
  !> type of this varies between the variants. Forth parameter returns
  !> the id of the variable. Fifth and sixth parameters are optional
  !> strings, which contain a comment and unit for the variable,
  !> respectively.
  !> The multidim variant differs from this by having instead of the
  !> data parameter two parameters. The first is a netcdf kind
  !> parameter, the second is an integer array determining the
  !> dimensions of the field.
  interface nc_define
    module procedure nc_define_int_0
    module procedure nc_define_int_1
    module procedure nc_define_int_2
    module procedure nc_define_double_0
    module procedure nc_define_double_1
    module procedure nc_define_double_2
    module procedure nc_define_long_2
    module procedure nc_define_multidim
    module procedure nc_define_char_1
  end interface nc_define


  !> \brief Interface to define variables with an unlimited dimension.
  !>
  !> Define an array with given type and dimensions and one unlimited
  !> dimension.
  interface nc_define_unlimited
    module procedure nc_define_multidim_unlimited
  end interface nc_define_unlimited


  !> \brief Interface for subroutine to get the size of an array.
  !>
  !> An interface to inquire the size of the dimensions of an 1 or 2
  !> dimensional array.
  !> The subroutines share the same basic format of the interface, whith
  !> three parameters. The first is an integer for the id of the netcdf
  !> file. The second is a character array, i.e. string, with the name
  !> of the field. The third is an integer (1d case) or an array with
  !> two elements (2d case), which will hold the size of the array after
  !> the call.
  interface nc_inq_dim
     module procedure nc_inq_dim_1
     module procedure nc_inq_dim_2
  end interface nc_inq_dim


  !> \brief Interface for subroutines to get id and sizes of an array.
  !>
  !> Contrary to nc_inq_dim these subroutines return the id of the field
  !> and the dimensions, but the latter as multiple variables, not as
  !> array.
  !> The common part of the argument list is first the id of the
  !> nc-file, then the name of the field for which one wants to get the
  !> bounds. Third variable returns the id of the field. Then there are
  !> two additional variables for each dimension. The first of these
  !> will return the lower and upper bound (in this order) of the first
  !> dimension of the array. The next pair for the second dimension,
  !> etc.
  interface nc_inq_id_bounds
    module procedure nc_inq_id_bounds_1
    module procedure nc_inq_id_bounds_2
  end interface nc_inq_id_bounds


  !> \brief Interface for subroutines to add a variable.
  !>
  !> This interface combines several subroutines to add a variable with
  !> with given value as group with given name. Optionally also a
  !> comment and a unit can be given.
  !> All subroutines share the type of Interface. First parameter is
  !> nc-file ID, second is name of the group, third the value(s) to
  !> store. As fourth and fifth argument a comment, and the unit of the
  !> value can be given. The third parameter is what varies between the
  !> individual routines in type and/or dimension.
  interface nc_add
    module procedure nc_add_int_0
    module procedure nc_add_double_0
    module procedure nc_add_double_1
    module procedure nc_add_double_2
    module procedure nc_add_char_1
  end interface nc_add


  !> \brief Interface for subroutines to get a variable.
  !>
  !> This interface combines several subroutines to get a variable with
  !> given name from a netcdf file, given as id.
  !> All subroutines share the type of interface. There are three
  !> parameters. The first is an integer, for the netcdf id of the file,
  !> from which to read. Second parameter is a 1D character array, i.e.
  !> a string, giving the name of the field to read. The third parameter
  !> is the variable in which to store the result, it varies based on
  !> the actual subroutine called. It varies in type (integer and
  !> double) and in the number of dimensions (0, i.e. scalar, 1d,...).
  interface nc_get
     module procedure nc_get_int_0
     module procedure nc_get_int_1
     module procedure nc_get_double_0
     module procedure nc_get_double_1
     module procedure nc_get_double_2
     module procedure nc_get_double_3
  end interface nc_get

contains

  !> \brief Define a group and get the id.
  !>
  !> input:
  !> ------
  !> ncid: integer, id of the file where to create the group.
  !> grpname: string, name of the group to create.
  !>
  !> output:
  !> -------
  !> ncid_grp: integer, id of the created group.
  !>
  !> sideeffects:
  !> ------------
  !> Changes the file.
  subroutine nc_defineGroup(ncid, grpname, ncid_grp)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: grpname
    integer, intent(out) :: ncid_grp

    call nf90_check(nf90_def_grp(ncid, trim(grpname), ncid_grp))

  end subroutine nc_defineGroup


  subroutine nc_define_int_0(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_INT, varid = varid))

    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_define_int_0


  subroutine nc_define_int_1(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, dimension(:), allocatable, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid

    if (allocated(var)) then
      call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
      call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

      call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
      if (present(comment)) then
        call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
      end if
      if (present(unit)) then
        call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
      end if
    end if
  end subroutine nc_define_int_1


  subroutine nc_define_int_2(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, dimension(:,:), allocatable, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid(2)

    if (allocated(var)) then
      call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
      call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
      call nf90_check(nf90_def_var(ncid, name, NF90_INT, dimid, varid))

      call nf90_check(nf90_put_att(ncid, varid, 'lbound1', lbound(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
      if (present(comment)) then
        call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
      end if
      if (present(unit)) then
        call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
      end if
    end if
  end subroutine nc_define_int_2


  subroutine nc_define_double_0(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, varid = varid))

    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_define_double_0


  subroutine nc_define_double_1(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:), allocatable, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid

    if (allocated(var)) then
      call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
      call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

      call nf90_check(nf90_put_att(ncid, varid, 'lbound', LBOUND(var)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound', UBOUND(var)))
      if (present(comment)) then
        call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
      end if
      if (present(unit)) then
        call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
      end if
    end if
  end subroutine nc_define_double_1


  subroutine nc_define_double_1_nonalloc(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:), intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid

    call nf90_check(nf90_def_dim(ncid, name // "_dim", size(var,1), dimid))
    call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

    call nf90_check(nf90_put_att(ncid, varid, 'lbound', lbound(var)))
    call nf90_check(nf90_put_att(ncid, varid, 'ubound', ubound(var)))
    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if

  end subroutine nc_define_double_1_nonalloc


  subroutine nc_define_double_2(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:,:), allocatable, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid(2)

    if (allocated(var)) then
      call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
      call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
      call nf90_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid, varid))

      call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
      if (present(comment)) then
        call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
      end if
      if (present(unit)) then
        call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
      end if
    end if
  end subroutine nc_define_double_2


  subroutine nc_define_long_2(ncid, name, var, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer(kind=8), dimension(:,:), allocatable, intent(in) :: var
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid(2)

    if (allocated(var)) then
      call nf90_check(nf90_def_dim(ncid, name // "_dim1", size(var,1), dimid(1)))
      call nf90_check(nf90_def_dim(ncid, name // "_dim2", size(var,2), dimid(2)))
      call nf90_check(nf90_def_var(ncid, name, NF90_INT64, dimid, varid))

      call nf90_check(nf90_put_att(ncid, varid, 'lbound1', LBOUND(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound1', UBOUND(var,1)))
      call nf90_check(nf90_put_att(ncid, varid, 'lbound2', LBOUND(var,2)))
      call nf90_check(nf90_put_att(ncid, varid, 'ubound2', UBOUND(var,2)))
      if (present(comment)) then
        call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
      end if
      if (present(unit)) then
        call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
      end if
    end if
  end subroutine nc_define_long_2


  subroutine nc_define_char_1(ncid, name, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    call nf90_check(nf90_def_var(ncid, name, NF90_STRING, varid = varid))

    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_define_char_1


  subroutine nc_define_multidim(ncid, name, field_type, dims, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: field_type
    integer, dimension(:), intent(in) :: dims
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer, dimension(:), allocatable :: dimid
    integer :: i
    character(len=32) :: dimname
    allocate(dimid(size(dims)))

    do i=1,size(dims)
      write (dimname, '(A,A,I1)') name, '_dim', i
      call nf90_check(nf90_def_dim(ncid, dimname, dims(i), dimid(i)))
    end do
    call nf90_check(nf90_def_var(ncid, name, field_type, dimid, varid))

    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if

    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if

    deallocate(dimid)

  end subroutine nc_define_multidim


  subroutine nc_define_multidim_unlimited(ncid, name, field_type, varid, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: field_type
    integer, intent(out) :: varid
    character(len=*), intent(in), optional :: comment, unit

    integer :: dimid

    call nf90_check(nf90_def_dim(ncid, name // '_dim', NF90_UNLIMITED, dimid))
    call nf90_check(nf90_def_var(ncid, name, field_type, dimid, varid))
    if (present(comment)) then
      call nf90_check(nf90_put_att(ncid, varid, "comment", comment))
    end if
    if (present(unit)) then
      call nf90_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
  end subroutine nc_define_multidim_unlimited


  subroutine nc_add_int_0(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_int_0(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_int_0


  subroutine nc_add_double_0(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_double_0(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_double_0


  subroutine nc_add_double_1(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:), allocatable, intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_double_1(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_double_1


  subroutine nc_add_double_1_nonalloc(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:), intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_double_1_nonalloc(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_double_1_nonalloc


  subroutine nc_add_double_2(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    double precision, dimension(:,:), allocatable, intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_double_2(ncid, name, var, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_double_2


  subroutine nc_add_char_1(ncid, name, var, comment, unit)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: var
    character(len=*), intent(in), optional :: comment, unit

    integer :: ierr, varid

    ierr = nf90_redef(ncid)
    call nc_define_char_1(ncid, name, varid, comment, unit)
    ierr = nf90_enddef(ncid)
    call nf90_check(nf90_put_var(ncid, varid, var))

  end subroutine nc_add_char_1


  subroutine nc_get_int_0(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    integer, intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_int_0

  subroutine nc_get_int_1(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    integer, dimension(:), intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_int_1

  subroutine nc_get_double_0(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_0

  subroutine nc_get_double_1(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:), intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_1

  subroutine nc_get_double_2(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:), intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_2

  subroutine nc_get_double_3(ncid, name, var)
    integer :: ncid
    character(len=*) :: name
    double precision, dimension(:,:,:), intent(out) :: var
    integer :: varid

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_var(ncid, varid, var))
  end subroutine nc_get_double_3


  !> \brief Check if a group exists already in the netcdf file.
  !>
  !> input:
  !> ------
  !> ncid: integer, id of the file to check for the group.
  !> name: string, name of the group to check for.
  !>
  !> output:
  !> -------
  !> grp_ncid: integer, if the group is found, this contains the id of
  !>   the group. If the group is not found, the value is undefined.
  !> found: logical, optional. If present will be set to true if the
  !>   group was found, and false if not. If not present, then usual
  !>   error check will be applied.
  subroutine nc_inquire_group(ncid, name, grp_ncid, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: grp_ncid
    logical, intent(out), optional :: found

    integer :: ierr

    ierr = nf90_inq_ncid(ncid, trim(name), grp_ncid)
    if (present(found)) then
      found = (ierr == NF90_NOERR)
    else
      call nf90_check(ierr)
    end if

  end subroutine nc_inquire_group


  subroutine nc_inq_dim_1(ncid, name, len)
    integer :: ncid, varid
    character(len=*)      :: name
    integer, dimension(1) :: dimids
    integer :: len

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_inquire_variable(ncid, varid, dimids = dimids))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(1), len = len))
  end subroutine nc_inq_dim_1

  subroutine nc_inq_dim_2(ncid, name, len)
    integer :: ncid, varid
    character(len=*)      :: name
    integer, dimension(2) :: dimids
    integer, dimension(2) :: len

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_inquire_variable(ncid, varid, dimids = dimids))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(1), len = len(1)))
    call nf90_check(nf90_inquire_dimension(ncid, dimids(2), len = len(2)))
  end subroutine nc_inq_dim_2


  subroutine nc_inq_id_bounds_1(ncid, name, varid, lb, ub)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid, lb, ub

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_att(ncid, varid, "lbound", lb))
    call nf90_check(nf90_get_att(ncid, varid, "ubound", ub))
  end subroutine nc_inq_id_bounds_1


  subroutine nc_inq_id_bounds_2(ncid, name, varid, lb1, ub1, lb2, ub2)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: varid, lb1, ub1, lb2, ub2

    call nf90_check(nf90_inq_varid(ncid, name, varid))
    call nf90_check(nf90_get_att(ncid, varid, "lbound1", lb1))
    call nf90_check(nf90_get_att(ncid, varid, "ubound1", ub1))
    call nf90_check(nf90_get_att(ncid, varid, "lbound2", lb2))
    call nf90_check(nf90_get_att(ncid, varid, "ubound2", ub2))
  end subroutine nc_inq_id_bounds_2


  !> \brief Find/create group with given name.
  !>
  !> Check if a netcdf file contains a given group, and create the group
  !> if not. In both cases return the group id.
  !>
  !> input:
  !> ------
  !> ncid: integer, id of the netcdf file to check.
  !> name: string, name of the group for which to check.
  !>
  !> output:
  !> -------
  !> grpid: integer, will contain the id of the requested group.
  !> found: logical, will be true if the group with 'name' did
  !>   already exist, false if not.
  !>
  !> sideeffects:
  !> ------------
  !> Creates the group if it does not already exist.
  subroutine nc_findGroup(ncid, name, grpid, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out) :: grpid
    logical, intent(out) :: found
    character(len=256) :: filename
    integer :: ierr

    found = .true.
    ierr = nf90_inq_ncid(ncid, trim(name), grpid);
    if (ierr /= NF90_NOERR) then
       write(filename,'(100A)') trim(adjustl(name)), '.nc'
       call nf90_check(nf90_open(filename, NF90_NOWRITE, grpid))
       found = .false.
    end if

  end subroutine nc_findGroup


  !> \brief Subroutine for checking status of operation.
  !>
  !> Check the status of an netcdf operation. In case of an error, the
  !> error message is printed via netcdf and the program stoped.
  !> An optional argument gives the possibility to ignore errors.
  !>
  !> input:
  !> ------
  !> status: integer, return/status value of an netcdf subroutine/
  !>   operation, that should be checked.
  !> optException: logical, optional, when false treat an error not as
  !>   an exception (in the programming sense), i.e. continue without
  !>   error message. [.true.]
  subroutine nf90_check(status, optException)
    integer, intent(in) :: status
    logical, intent(in), optional :: optException
    logical :: exception

    exception = .true.
    if (present(optException)) exception = optException

    if(status /= nf90_noerr) then
      if (exception) then
        print *, trim(nf90_strerror(status))
        stop
      end if
    end if
  end subroutine nf90_check

  !> \brief Wrapper for closing a netcdf file.
  !>
  !> input:
  !> ------
  !> ncid: netcdf id, got from e.g. nc_open, for the file to close.
  !>
  !> output:
  !> -------
  !> none
  !>
  !> sideeffects:
  !> none intended, but there might be some due to error checking.
  subroutine nc_close(ncid)
    integer, intent(in) :: ncid

    call nf90_check(nf90_close(ncid))
  end subroutine nc_close

  !> \brief Wrapper for opening a netcdf file in read mode.
  !>
  !> input:
  !> ------
  !> filename: character array, with the name(+path) of the file to open.
  !> optException: logical, optional, if false then errors are not
  !>   treated as exceptions.
  !>
  !> output:
  !> -------
  !> ncid: integer, netcdf id of the opened file, that can be used to
  !>   inquire information or get data.
  !>
  !> sideeffects:
  !> none intended, but there might be some due to error checking.
  subroutine nc_open(filename, ncid, optException)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ncid
    logical, intent(in), optional :: optException

    logical :: exception

    exception = .true.
    if (present(optException)) exception = optException

    call nf90_check(nf90_open(filename, NF90_NOWRITE, ncid), exception)
  end subroutine nc_open


  !> \brief Create an empty netcdf-4 file.
  !>
  !> input:
  !> ------
  !> filename: string, name of the file to create.
  !> fileformat_version: string, optional, version information for the
  !>   file you are creating (not netcdf version). Added as an attribute
  !>   of the created file. ['1.0']
  !>
  !> output:
  !> -------
  !> ncid: integer, netcdf id of the opened file, that can be used to
  !>   inquire information or get data.
  !>
  !> sideeffects:
  !> ------------
  !> Creates file on the disk.
  subroutine nc_create(filename, ncid, fileformat_version)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ncid
    character(len=*), intent(in out), optional :: fileformat_version

    if (.not. present(fileformat_version)) then
       fileformat_version = '1.0'
    end if

    write (*,*) "Creating NetCDF-4 File: ", filename
    call nf90_check(nf90_create(filename, NF90_NETCDF4, ncid))
    call nf90_check(nf90_put_att(ncid, NF90_GLOBAL, 'Version', fileformat_version))
  end subroutine nc_create


  !> \brief Wrapper for nf90_enddef.
  !>
  !> With check of error code.
  !>
  !> input:
  !> ------
  !> ncid: integer, id of the file for which to do the operation.
  !>
  !> sideeffects:
  !> ------------
  !> Changes state of file.
  subroutine nc_enddef(ncid)
    integer, intent(in) :: ncid

    call nf90_check(nf90_enddef(ncid))
  end subroutine nc_enddef

end module nctools_module
