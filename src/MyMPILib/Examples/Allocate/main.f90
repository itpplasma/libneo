program allocatetest
  use mpiprovider_module
  implicit none
  
  integer :: MPI_COMM_WORKERS
  integer :: myrank, ierr
  double precision, dimension(:), allocatable :: testmat
  integer :: multiplier
  
  ! Initialize MPI (mpro is the singleton name of the MPIProvider)
  !call mpro%init()

  myrank = mpro%getRank()
  if (myrank .eq. 2) then
     write (*,*) myrank
     multiplier = 1d6
  else
     multiplier = 1d4
  end if
  allocate(testmat(1000*multiplier*(1)))
  testmat = 1.1d0
  call sleep(120)
  !write (*,*) testmat(1)
  !deallocate(testmat)
  
  ! Close all connections by the use of the MPIProvider
  !call mpro%deinit(.false.)

end program allocatetest
