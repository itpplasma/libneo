!> Demonstration program for numerical integration
program simpleintegrate
  use mpiprovider_module
  implicit none
  
  integer :: MPI_COMM_WORKERS
  double complex, dimension(1:5, 1:4, 0:3) :: A
  double precision, dimension(1:5) :: B
  integer :: myrank, ierr
  
  ! Initialize MPI (mpro is the singleton name of the MPIProvider)
  call mpro%init()

  !if (mpro%getNumProcs() .eq. 1) then

     myrank = mpro%getRank()
     B = 0d0
     B(myrank) = myrank*10
     !A = 0
     !A(:,:,myrank) = 10*(myrank+1)
     !A(:,:,myrank) = 10*myrank

     !call MPI_ALLGATHER(A(1:5, myrank), 5, MPI_DOUBLE_COMPLEX, A, 5, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
     !call mpro%allgather(A(:,:,myrank), A)
     call mpro%allgather(B(myrank), B)
     write (*,*) myrank, B
     stop
     !write (*,*) myrank, A(2,2,:)

     ! Try a spawn
     !call MPI_Comm_spawn('./gathertest', MPI_ARGV_NULL, 3, &
     !MPI_INFO_NULL, 0, MPI_COMM_WORLD, MPI_COMM_WORKERS, & 
     !MPI_ERRCODES_IGNORE, ierr)

     !call mpro%barrier()
  !else

     !write (*,*) "This seems to be a spwaned worker"
     !write (*,*) "I am ", mpro%getRank(), mpro%getNumProcs()

     !call mpro%barrier()
  !end if

  ! Close all connections by the use of the MPIProvider
  call mpro%deinit(.false.)

end program simpleintegrate
