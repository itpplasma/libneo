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
  
  subroutine nf90_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop
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
    integer :: ncid
    
    call nf90_check(nf90_close(ncid))
  end subroutine nc_close

  !> \brief Wrapper for opening a netcdf file in read mode.
  !>
  !> input:
  !> ------
  !> filename: character array, with the name(+path) of the file to open.
  !>
  !> output:
  !> -------
  !> ncid: integer, netcdf id of the opened file, that can be used to
  !>   inquire information or get data.
  !>
  !> sideeffects:
  !> none intended, but there might be some due to error checking.
  subroutine nc_open(filename, ncid)
    character(len=*) :: filename
    integer :: ncid

    call nf90_check(nf90_open(filename, NF90_NOWRITE, ncid))
  end subroutine nc_open

end module nctools_module
