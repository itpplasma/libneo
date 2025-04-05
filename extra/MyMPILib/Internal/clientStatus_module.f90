!> Module for class clientStatus
module clientStatus_module
  use workunit_module
  use intlist_module
  use wulist_module

  implicit none

  !> Class clientStatus:
  !> Stores information per client
  type clientStatus

     !> Defines the rank of the process
     integer :: rank = -1

     !> True, if client has no Workunits pending
     logical :: isReady = .true.

     !> For blocking client when it should recveive from a dataRequester
     logical :: isBlocked = .false.

     logical :: isDone = .false.

     !> If working with spools, the number of workunits pending on the client
     integer :: wuSpool = 0

     !> For merge process of Sched1
     integer :: nextNeighbor

     !> For merge process of Sched1: Defines which nodes have to be locally merged
     type(workunitList) :: pendingNodes
    
     !> For performance analysis
     double precision :: lastTime = -1

     !> Mean time, the client needed for a workunit
     double precision :: meanWorkTime = 0

     !> Units per client, used to calculate mergetrees
     integer :: upc

     !> Defines the workunits, the client has in local memory
     type(intList) :: myWorkunits

     !> Workunits, which have been sent to the client, but not processed, yet
     type(intList) :: lastWUs

     integer :: workunitsLeft = -1

     !> Required for Scheduler 2 and 3 to indicate that the client has local merge possibilities, before requesting workunits from other clients
     logical :: localMerge = .false.

   contains

     !> Destructor
     procedure :: free => free_clientStatus
  end type clientStatus

contains

   !> Destructor
  subroutine free_clientStatus(this)
    class(clientStatus) :: this
    
    call this%myWorkunits%free();
    call this%pendingNodes%free(.false.);
    call this%lastWUs%free();
  end subroutine free_clientStatus

end module clientStatus_module
