module nctools_module

  use netcdf
  implicit none

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
