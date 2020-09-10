!> Module for Singleton class mpiprovider
module mpiprovider_module

  !Include OpenMPI Libs
  use mpi

  use commandline_parser_module
  use configFile_parser_module
  use wuList_module
  use mpelog_module
  use mylog_module
  use packBuffer_module

  implicit none

  ! Some constants
  integer, parameter :: wuTag     = 1
  integer, parameter :: comTag    = 2
  include "../Internal/version.f90"

  !> Storage structure for workunits. Belongs to the scheduler, but in the singleton each workunit has access to it
  type :: globalSchedulerStorage
    ! Main lists for workunit processing
    type(workunitlist), allocatable :: waitingWorkunits   !< Workunits which have not been sent to a client yet
    type(workunitlist), allocatable :: pendingWorkunits   !< Workunits which are already on clients
    type(workunitlist), allocatable :: processedWorkunits !< Workunits which are completely done

    !> UID for next packable object
    integer    :: nextUID = 1

  end type globalSchedulerStorage

  !> Singleton class for communicating with MPI
  type, private :: mpiprovider

    !> Some MPI status variables
    integer :: ierr, rank, numprocs
    integer :: mpi_comm_clients

    !> Current procname of process
    character(len=MPI_MAX_PROCESSOR_NAME) :: procname

    !> Activate/Deactivate MPE-profiling
    logical :: mpe_log = .true.

    !> The buffer for packing and unpacking
    type(packBuffer), public :: packBuffer

    !> Storage for workunit lists
    type(globalSchedulerStorage) :: storage

    !> Define the initWorkunit
    class(workunit), pointer :: initWU => null()

    !> Performance analysis
    double precision :: wuTime = 0
    integer :: wuCount = 0
    double precision :: initTime = 0
    double precision :: deinitTime = 0
    double precision :: schedInitTime = 0
    double precision :: schedDeinitTime = 0
    integer :: balanceCount = 0
    double precision :: drTime = 0
    double precision :: drPackTime = 0
    integer :: drCount = 0

    !> Benchmarks
    integer :: meanDurWu      !< How long should a workunit take (mean value)
    integer :: meanDurInitWu  !< How long should the creating of an initWU-value take (mean value)
    real  :: clientSpeed      !< How fast should the client be (for simulating faster and slower clients)
  contains
    !> Constructor of mpiprovider
    procedure :: init => init_mpiprovider

    !> Allocate send and receive buffer
    procedure :: allocateBuffers => allocateBuffers_mpiprovider

    !> Destructor of mpiprovider
    procedure :: deinit => deinit_mpiprovider

    !> Some getter functions for MPI status variables
    procedure :: getProcName => getProcName_mpiprovider
    procedure :: getRank => getRank_mpiprovider
    procedure :: getNumProcs => getNumProcs_mpiprovider
    procedure :: isParallel => isParallel_mpiprovider
    procedure :: isMaster => isMaster_mpiprovider

    !> Wrapper commands
    procedure :: send => send_mpiprovider
    procedure :: recv => recv_mpiprovider
    procedure :: isend => isend_mpiprovider

    procedure :: allgather_complex_1 => allgather_complex_1_mpiprovider
    procedure :: allgather_complex_2 => allgather_complex_2_mpiprovider
    procedure :: allgather_complex_3 => allgather_complex_3_mpiprovider
    procedure :: allgather_complex_4 => allgather_complex_4_mpiprovider
    procedure :: allgather_double_1 => allgather_double_1_mpiprovider
    procedure :: allgather_double_2 => allgather_double_2_mpiprovider
    procedure :: allgather_double_3 => allgather_double_3_mpiprovider
    procedure :: allgather_double_4 => allgather_double_4_mpiprovider
    procedure :: allgather_double_5 => allgather_double_5_mpiprovider


    procedure :: allgather_inplace_complex_1 => allgather_inplace_complex_1_mpiprovider
    procedure :: allgather_inplace_complex_2 => allgather_inplace_complex_2_mpiprovider
    procedure :: allgather_inplace_complex_3 => allgather_inplace_complex_3_mpiprovider
    procedure :: allgather_inplace_complex_4 => allgather_inplace_complex_4_mpiprovider
    procedure :: allgather_inplace_complex_5 => allgather_inplace_complex_5_mpiprovider
    procedure :: allgather_inplace_double_1 => allgather_inplace_double_1_mpiprovider
    procedure :: allgather_inplace_double_2 => allgather_inplace_double_2_mpiprovider
    procedure :: allgather_inplace_double_3 => allgather_inplace_double_3_mpiprovider
    procedure :: allgather_inplace_double_4 => allgather_inplace_double_4_mpiprovider
    procedure :: allgather_inplace_double_5 => allgather_inplace_double_5_mpiprovider

    generic, public :: allgather => allgather_complex_1, allgather_complex_2, &
         allgather_complex_3, allgather_complex_4, &
         allgather_double_1, allgather_double_2, &
         allgather_double_3, allgather_double_4, &
         allgather_double_5

   generic, public :: allgather_inplace => &
        allgather_inplace_complex_1, allgather_inplace_complex_2, &
        allgather_inplace_complex_3, allgather_inplace_complex_4, &
        allgather_inplace_complex_5, &
        allgather_inplace_double_1, allgather_inplace_double_2, &
        allgather_inplace_double_3, allgather_inplace_double_4, &
        allgather_inplace_double_5

    procedure :: barrier => barrier_mpiprovider

    !> For benchmark analysis
    procedure :: randomSleep => randomSleep_mpiprovider
    procedure, private :: initRandomSeed => initRandomSeed_mpiprovider

  end type mpiprovider

  !> Singleton
  type(mpiprovider) :: mpro

contains

  !> Constructor of mpiprovider
  subroutine init_mpiprovider(this)
    class(mpiprovider) :: this
    character(len=MPI_MAX_PROCESSOR_NAME) :: procname
    integer :: ierr, myid, namelen, numprocs
    integer :: mpi_world_group, mpi_client_group
    integer, dimension(1) :: group_excl_ranks

    ! Initialize MPI
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "An error occurred in MPI_INIT", ierr
      stop
    end if

    ! Get number of processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "An error occurred in MPI_COMM_SIZE", ierr
      stop
    end if

    ! Get rank
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "An error occurred in MPI_COMM_RANK", ierr
      stop
    end if

    ! Get processor name
    call MPI_GET_PROCESSOR_NAME(procname, namelen, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "An error occurred in MPI_GET_PROCESSOR_NAME", ierr
      stop
    end if

    ! Create communicator for clients
    if (numprocs > 1) then
        group_excl_ranks(1) = 0
        call MPI_COMM_GROUP(MPI_COMM_WORLD, mpi_world_group, ierr)
        call MPI_GROUP_EXCL(mpi_world_group, 1, group_excl_ranks, mpi_client_group, ierr)
        call MPI_COMM_CREATE(MPI_COMM_WORLD, mpi_client_group, this%mpi_comm_clients, ierr)
        call MPI_GROUP_FREE(mpi_client_group, ierr)
        call MPI_GROUP_FREE(mpi_world_group, ierr)
    end if

    ! Set local variables
    this%ierr = ierr
    this%procname = procname
    this%rank = myid
    this%numprocs = numprocs

    this%initTime = MPI_WTime()

    ! Initialize packBuffer
    this%packBuffer = packBuffer() !New constructor possible since gfortran-4.7

    ! Initialize random number generator
    call this%initRandomSeed()

    ! Display state
    write (*,"(A, I3, A, A10, A)") "MPI-Provider on process number ", myid, " on node ", trim(procname), " initialized!"

  ! Print version information
  if (myid == 0) then
     write (*,*) ''
     write (*,*) "--------- MyMPILib Git Revision --------"
     write (*,*) MyMPILib_Version
     write (*,*) "----------------------------------------"
     write (*,*) ''
     if (len_trim(MyMPILib_Version_Additional) /= 0) then
          write (*,*) "################################### MyMPILib Git Additional Information ##################################"
          write (*,*) MyMPILib_Version_Additional
          write (*,*) "##########################################################################################################"
          write (*,*) ''
     end if
  end if
  end subroutine init_mpiprovider

  !> Allocate the buffers, set normBufferSize before calling this function
  subroutine allocateBuffers_mpiprovider(this, buffersize)
    class(mpiprovider) :: this
    integer :: buffersize

    ! Are buffers required or sequential mode?
    if (this%numprocs > 1) then
      if (this%getRank() == 0) write (*,*) "Allocating ", buffersize, "MBytes for the pack-buffers."
      call this%packBuffer%allocateBuffers(buffersize)
    end if

  end subroutine allocateBuffers_mpiprovider

  !> Destructor of mpiprovider
  subroutine deinit_mpiprovider(this, showStats_opt)
    class(mpiprovider) :: this
    logical, optional  :: showStats_opt
    logical            :: showStats
    integer :: ierr

    call this%packBuffer%clear()
    if (this%numprocs > 1) then
      if (allocated(this%packBuffer%sendBuffer)) deallocate(this%packBuffer%sendBuffer)
      if (allocated(this%packBuffer%recvBuffer)) deallocate(this%packBuffer%recvBuffer)
    end if

    this%deinitTime = MPI_WTime()

    showStats = .true.
    if (present(showStats_opt)) showStats = showStats_opt

    if (mpro%isMaster() .and. showStats) then
      write (*,*) "Runtime analysis"
      write (*,*) "----------------"
      write (*,"(A, F16.2, A)") "Complete runtime:       ", this%deinitTime - this%initTime, " s"
      if (this%numprocs > 1) then
         write (*,"(A, F16.2, A)") "Scheduler runtime:      ", this%scheddeinitTime - this%schedinitTime, " s"
         write (*,"(A, F16.2, A)") "Time before scheduling: ", this%schedInitTime - this%initTime, " s"
         write (*,"(A, F16.2, A)") "Time after scheduling:  ", this%deinitTime - this%scheddeinitTime, " s"
      end if
    end if

    call mpi_finalize(ierr)
    this%ierr = ierr

  end subroutine deinit_mpiprovider

  !> Wrapper function for isend
  function isend_mpiprovider(this, dest, tag, buffer) result(res)
    class(mpiprovider) :: this
    integer, intent(in) :: dest, tag
    character(len=*) :: buffer
    integer :: ierr
    integer :: res

    call mlog%logEvent(mpe_e_sendA)
    call MPI_ISEND(buffer, len(buffer), MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, res, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ISEND", ierr
      stop
    end if
    call mlog%logEvent(mpe_e_sendB, dest, tag)
  end function isend_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_double_1_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double precision   :: sendbuf
    double precision, dimension(:) :: recvbuf
    integer :: sendcount, recvcount
    integer :: ierr

    sendcount = 1
    recvcount = sendcount

    !write (*,*) "This is allgather_double_1_mpiprovider"
    !write (*,*) "Sendbuf: ", sendbuf
    !write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
    !write (*,*) "Recvbuf after gather: ", recvbuf

  end subroutine allgather_double_1_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_double_2_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:)   :: sendbuf
    double precision, dimension(:,:) :: recvbuf
    integer, dimension(1) :: lbound_send, ubound_send
    integer, dimension(2) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = ubound_send(1)-lbound_send(1)+1
    recvcount = sendcount

   ! write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ALLGATHER", ierr
      stop
    end if
  end subroutine allgather_double_2_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_double_3_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:)   :: sendbuf
    double precision, dimension(:,:,:) :: recvbuf
    integer, dimension(2) :: lbound_send, ubound_send
    integer, dimension(3) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = (ubound_send(2)-lbound_send(2)+1) * (ubound_send(1)-lbound_send(1)+1)
    recvcount = sendcount

   ! write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ALLGATHER", ierr
      stop
    end if
  end subroutine allgather_double_3_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_double_4_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:,:)   :: sendbuf
    double precision, dimension(:,:,:,:) :: recvbuf
    integer, dimension(3) :: lbound_send, ubound_send
    integer, dimension(4) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = (ubound_send(3)-lbound_send(3)+1) * &
                (ubound_send(2)-lbound_send(2)+1) * &
                (ubound_send(1)-lbound_send(1)+1)
    recvcount = sendcount

    !write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ALLGATHER", ierr
      stop
    end if
  end subroutine allgather_double_4_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_double_5_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:,:,:)   :: sendbuf
    double precision, dimension(:,:,:,:,:) :: recvbuf
    integer, dimension(4) :: lbound_send, ubound_send
    integer, dimension(5) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = &
         (ubound_send(4)-lbound_send(4)+1) * &
         (ubound_send(3)-lbound_send(3)+1) * &
         (ubound_send(2)-lbound_send(2)+1) * &
         (ubound_send(1)-lbound_send(1)+1)

    recvcount = sendcount

    !write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_double_5_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_complex_1_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double complex   :: sendbuf
    double complex, dimension(:) :: recvbuf
    integer :: sendcount, recvcount
    integer :: ierr

    sendcount = 1
    recvcount = sendcount

    !write (*,*) "This is allgather_double_1_mpiprovider"
    !write (*,*) "Sendbuf: ", sendbuf
    !write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
    !write (*,*) "Recvbuf after gather: ", recvbuf

  end subroutine allgather_complex_1_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_complex_2_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:)   :: sendbuf
    double complex, dimension(:,:) :: recvbuf
    integer, dimension(1) :: lbound_send, ubound_send
    integer, dimension(2) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = ubound_send(1)-lbound_send(1)+1
    recvcount = sendcount

    ! write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_complex_2_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_complex_3_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:)   :: sendbuf
    double complex, dimension(:,:,:) :: recvbuf
    integer, dimension(2) :: lbound_send, ubound_send
    integer, dimension(3) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = (ubound_send(2)-lbound_send(2)+1) * (ubound_send(1)-lbound_send(1)+1)
    recvcount = sendcount

   ! write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ALLGATHER", ierr
      stop
    end if
  end subroutine allgather_complex_3_mpiprovider

  !> Wrapper function for allgather
  subroutine allgather_complex_4_mpiprovider(this, sendbuf, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:,:)   :: sendbuf
    double complex, dimension(:,:,:,:) :: recvbuf
    integer, dimension(3) :: lbound_send, ubound_send
    integer, dimension(4) :: lbound_recv, ubound_recv
    integer :: sendcount, recvcount
    integer :: ierr

    lbound_send = lbound(sendbuf)
    ubound_send = ubound(sendbuf)
    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    sendcount = (ubound_send(3)-lbound_send(3)+1) * (ubound_send(2)-lbound_send(2)+1) * (ubound_send(1)-lbound_send(1)+1)
    recvcount = sendcount

    !write (*,*) ubound_send, lbound_send

    call MPI_ALLGATHER(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_ALLGATHER", ierr
      stop
    end if
  end subroutine allgather_complex_4_mpiprovider

!> BEGIN wrappers for allgather with MPI_IN_PLACE

  subroutine allgather_inplace_double_1_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:) :: recvbuf
    integer :: recvcount
    integer :: ierr

    recvcount = 1

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_double_1_mpiprovider

  subroutine allgather_inplace_double_2_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:) :: recvbuf
    integer, dimension(2) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = ubound_recv(1)-lbound_recv(1)+1

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_double_2_mpiprovider

  subroutine allgather_inplace_double_3_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:,:) :: recvbuf
    integer, dimension(3) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_double_3_mpiprovider

  subroutine allgather_inplace_double_4_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:,:,:) :: recvbuf
    integer, dimension(4) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(3)-lbound_recv(3)+1) * &
                (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_double_4_mpiprovider

  subroutine allgather_inplace_double_5_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double precision, dimension(:,:,:,:,:) :: recvbuf
    integer, dimension(5) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(4)-lbound_recv(4)+1) * &
                (ubound_recv(3)-lbound_recv(3)+1) * &
                (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_double_5_mpiprovider

  subroutine allgather_inplace_complex_1_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:) :: recvbuf
    integer :: recvcount
    integer :: ierr

    recvcount = 1

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_complex_1_mpiprovider

  subroutine allgather_inplace_complex_2_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:) :: recvbuf
    integer, dimension(2) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = ubound_recv(1)-lbound_recv(1)+1

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_complex_2_mpiprovider

  subroutine allgather_inplace_complex_3_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:,:) :: recvbuf
    integer, dimension(3) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_complex_3_mpiprovider

  subroutine allgather_inplace_complex_4_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:,:,:) :: recvbuf
    integer, dimension(4) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(3)-lbound_recv(3)+1) * &
                (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_complex_4_mpiprovider

  subroutine allgather_inplace_complex_5_mpiprovider(this, recvbuf)
    class(mpiprovider) :: this
    double complex, dimension(:,:,:,:,:) :: recvbuf
    integer, dimension(5) :: lbound_recv, ubound_recv
    integer :: recvcount
    integer :: ierr

    lbound_recv = lbound(recvbuf)
    ubound_recv = ubound(recvbuf)

    recvcount = (ubound_recv(4)-lbound_recv(4)+1) * &
                (ubound_recv(3)-lbound_recv(3)+1) * &
                (ubound_recv(2)-lbound_recv(2)+1) * &
                (ubound_recv(1)-lbound_recv(1)+1)

    call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, &
                       MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_ALLGATHER", ierr
       stop
    end if
  end subroutine allgather_inplace_complex_5_mpiprovider

!> END wrappers for allgather with MPI_IN_PLACE


  subroutine barrier_mpiprovider(this)
    class(mpiprovider) :: this
    integer :: ierr

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write (*,*) "Error in MPI_BARRIER", ierr
       stop
    end if
  end subroutine barrier_mpiprovider

  !> Wrapper function for send
  subroutine send_mpiprovider(this, dest, tag, buffer)
    class(mpiprovider) :: this
    integer, intent(in) :: dest, tag
    character(len=*) :: buffer
    integer :: ierr

    call mlog%logEvent(mpe_e_sendA)
    call MPI_SEND(buffer, len(buffer), MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_SEND", ierr
      stop
    end if
    call mlog%logEvent(mpe_e_sendB, dest, tag)
  end subroutine send_mpiprovider

  !> Wrapper function for receiving
  subroutine recv_mpiprovider(this, source, tag, buffer)
    class(mpiprovider) :: this
    integer, intent(in) :: tag
    integer, intent(out) :: source
    character(len=4), intent(out) :: buffer
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    call mlog%logEvent(mpe_e_recvA)
    call MPI_RECV(buffer, 4, MPI_CHARACTER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, rstatus, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,*) "Error in MPI_RECV", ierr
      stop
    end if
    source = rstatus(MPI_SOURCE)
    call mlog%logEvent(mpe_e_recvB, source, tag)

  end subroutine recv_mpiprovider

  !> Get name of current process
  function getProcName_mpiprovider(this)
    class(mpiprovider) :: this
    character(len=MPI_MAX_PROCESSOR_NAME) :: getProcName_mpiprovider

    getProcName_mpiprovider = trim(this%procname)
  end function getProcName_mpiprovider

  !> Get current rank
  function getRank_mpiprovider(this)
    class(mpiprovider) :: this
    integer :: getRank_mpiprovider

    getRank_mpiprovider = this%rank
  end function getRank_mpiprovider

  !> Get number of all processes
  function getNumProcs_mpiprovider(this)
    class(mpiprovider) :: this
    integer :: getNumProcs_mpiprovider

    getNumProcs_mpiprovider = this%numprocs
  end function getNumProcs_mpiprovider

  function isParallel_mpiprovider(this)
    class(mpiprovider) :: this
    logical :: isParallel_mpiprovider

    isParallel_mpiprovider = this%numprocs > 1
  end function isParallel_mpiprovider

  function paraMode_mpiprovider(this) result(res)
    class(mpiprovider) :: this
    logical :: res

    res = (this%numprocs > 1)
  end function paraMode_mpiprovider

  function isMaster_mpiprovider(this) result(res)
    class(mpiprovider) :: this
    logical :: res

    res = (this%getRank() == 0)
  end function isMaster_mpiprovider

  !> Setting a different random seed for each rank
  subroutine initrandomseed_mpiprovider(this)
    class(mpiprovider) :: this
    integer :: n
    integer, dimension(:), allocatable :: seed
    real :: r

    call random_seed(size = n)
    allocate(seed(n))

    ! Read a random seed from /dev/urandom
    open(50, file="/dev/urandom", access='stream',form='UNFORMATTED')
    read(50) seed
    close(50)

    ! Set new seed
    call random_seed(put = seed)

    deallocate(seed)

    call random_number(r)
    if (this%getNumProcs() > 2) then
      ! Used for 1. Benchmark
      this%clientSpeed = 1.0*r + 0.5
    else
      ! If running in linear mode, then always have the normal speed
      this%clientSpeed = 1.0
    end if
    !write (*,*) "Client performance of", mpro%getRank(), "=", this%clientSpeed
  end subroutine initrandomseed_mpiprovider

  !> Simulates longer calculations for benchmarks
  subroutine randomSleep_mpiprovider(this, meanTime)
    class(mpiprovider) :: this
    integer :: meanTime
    integer :: waitTime
    real :: r

    if (meanTime > 0) then
      call random_number(r)
      !write (*,*) r
      waitTime = nint(2*r*meanTime * this%clientSpeed)
      call sleep(waitTime)

    end if
  end subroutine randomSleep_mpiprovider

end module mpiprovider_module
