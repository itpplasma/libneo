program test_hdf5_tools
  use hdf5_tools

  use hdf5_tools, only : h5_add, h5_get, h5_init, h5_deinit
  use hdf5_tools, only : h5_create, h5_open, h5_close

  implicit none

  character(len=*), parameter :: filename = "test_file.h5"
  character(len=*), parameter :: groupname = "testgroup"
  real, parameter :: tolerance = 1.0e-12

  integer(HID_T) :: h5id, h5grpid

  character(len=12) :: dataset_int = "data_integer"
  character(len=11) :: dataset_dou = "data_double"
  character(len=12) :: dataset_com = "data_complex"
  character(len=12) :: dataset_log = "data_logical"
  integer :: int_read, int_write = 4
  double precision, dimension(2,2) :: double_read, double_write = reshape((/2.0, 3.0, 5.0, -1.0 /), (/2, 2/))
  complex(kind=dcp), dimension(3) :: complex_read, complex_write
  logical :: logical_read, logical_write = .false.
  logical :: error_found = .false.

  complex_write(1)%re = 2.0
  complex_write(1)%im = 3.0
  complex_write(2)%re = 5.0
  complex_write(2)%im =-1.0
  complex_write(3)%re = 0.4
  complex_write(3)%im = 1.3e-3

  call h5_init()

  ! Create file and write data.
  call h5_create(filename, h5id)

  call h5_add(h5id, dataset_int, int_write, unit='5')

  call h5_define_group(h5id, groupname, h5grpid)

  call h5_add(h5grpid, dataset_dou, double_write, lbounds=(/1,1/), ubounds=(/2,2/), comment='Heinz', accuracy=1.0d-100)

  call h5_open_group(h5id, groupname, h5grpid)

  call h5_add(h5grpid, dataset_com, complex_write, lbounds=(/1/), ubounds=(/3/), comment='Karl', unit='kg / pc')

  call h5_close_group(h5grpid)

  call h5_add(h5id, dataset_log, logical_write)

  call h5_close(h5id)

  ! Open file again and read data.
  ! Note: different order of reads compared to writes is intentional.
  call h5_open(filename, h5id)

  call h5_get(h5id, dataset_log, logical_read)

  call h5_open_group(h5id, groupname, h5grpid)

  call h5_get(h5grpid, dataset_com, complex_read)

  call h5_get(h5grpid, dataset_dou, double_read)

  call h5_close_group(h5grpid)

  call h5_get(h5id, dataset_int, int_read)

  call h5_close(h5id)

  call h5_deinit()

  ! Make sure that data read is same as those written.
  if (int_write /= int_read) then
    write(*,*) 'STOP integer does not match.'
    error_found = .true.
  end if

  if (sum(abs(double_read - double_write)) > tolerance) then
    write(*,*) 'STOP double does not match.'
    error_found = .true.
  end if

  if (sum(abs(complex_read - complex_write)) > tolerance) then
    write(*,*) 'STOP complex does not match.'
    error_found = .true.
  end if

  if (logical_write .neqv. logical_read) then
    write(*,*) 'STOP logical does not match.'
    error_found = .true.
  end if

  if (error_found) then
    stop 'STOP error occured'
  end if

end program test_hdf5_tools
