  
  include "listener_header.f90"
  include "comListener_header.f90"
  include "wuListener_header.f90"

  include "dispatcher_header.f90"
  include "comDispatcher_header.f90"
  include "wuDispatcher_header.f90"

  !> Settings for the generic scheduler (hardcoded defaults - namelist removed for nvfortran bug)
  logical :: loadBalancing  = .false.
  integer :: bufferSize     = 32    !Megabytes
  integer :: verbose        = 0
  logical :: activateMPE         = .false.

  !> Basic scheduler class with matmul performance testing
  type :: scheduler
     type(comListener) :: tComListener    !< Listener for commands
     type(comDispatcher) :: tComDisp      !< Dispatcher for commands
     type(wuListener)  :: tWuListener     !< Listener for workunits
     type(wuDispatcher) :: tWuDisp        !< Dispatcher for workunits

     character(len=maxStrLen) :: configFilename = "./schedConfig.txt"!< Sets the filename for the namelist

     logical :: isRunning                 !< .true. when scheduling is running
     type(clientStatus), dimension(:), allocatable :: clientStats   !< This array stores informations for every client

     double precision :: proTime = 0      !< For performance analysis, marks the time spent in workunits

     integer :: lastWuId = -1             !< Stores the last workunit processed by the client, used for printing the result
     integer :: maxSpool = 1              !< Scheduler 1 is implemented with spool support to send workunits to clients, even when they are not ready
     logical :: balance = .false.         !< Scheduler 2 and 3 are implemented with loadBalancing support
     
     integer, allocatable, dimension(:) :: workunits_per_client
     
     
     ! Main lists for workunit processing
     type(workunitlist), allocatable :: waitingWorkunits   !< Workunits which have not been sent to a client yet
     type(workunitlist), allocatable :: pendingWorkunits   !< Workunits which are already on clients
     type(workunitlist), allocatable :: processedWorkunits !< Workunits which are completely done

     !> UID for next packable object
     integer    :: nextUID = 1

     contains
       procedure :: init => init_scheduler
       procedure :: deinit => deinit_scheduler
       procedure :: terminate => terminate_scheduler

       procedure :: addWorkunit => addWorkunit_scheduler
       procedure :: getPendingNodesCount => getPendingNodesCount_scheduler
       
       procedure :: allocateWU => allocateWU_scheduler
       procedure :: allocateMergeWU => allocateMergeWU_scheduler

       procedure :: initPhase => initPhase_scheduler
       procedure :: runInitWU => runInitWU_scheduler
       procedure :: prepare => prepare_scheduler
       procedure :: summarize => summarize_scheduler
       procedure :: schedule => schedule_scheduler
       procedure :: processWaitingWorkunits => processWaitingWorkunits_scheduler
       procedure :: rebuildWU => rebuildWU_scheduler
       procedure :: setClientStatus => setClientStatus_scheduler
       procedure :: cleanup => cleanup_scheduler
       procedure :: partNearlyFair => partNearlyFair_scheduler
       procedure :: repairNeighbors => repairNeighbors_scheduler

       procedure :: printLast => printLast_scheduler
       procedure :: printPerformanceAnalysis => printPerformanceAnalysis_scheduler

       ! This subroutines CAN BE overwritten for special behaviour
       procedure, private :: initMaster => initMaster_scheduler
       procedure, private :: loadSettings => loadSettings_scheduler
       procedure, private :: initClient => initClient_scheduler
       procedure, private :: initLinear => initLinear_scheduler
       procedure, private :: checkIfClientDone => checkIfClientDone_scheduler ! Corresponds to dynamicWorkunits()
       procedure, private :: loadBalancing => loadBalancing_scheduler

       ! This subroutines HAVE TO be inherited when creating a new scheduler with new workunits
       procedure :: allocateSpecific => allocateSpecific_scheduler
       procedure :: allocateSpecificMergeWU => allocateSpecificMergeWU_scheduler
  end type scheduler
