!> Module for abstract class workunit, for documentation see genericWorkunit
module workunit_module

  use packable_module
  use intList_module

  implicit none
  
  type, extends(packable), abstract :: workunit
    real    :: fracIndex        !< Required for special purposes (correct sort order for matrix multiplication)
    character(len=1024) :: type = ''        !< Defines the type of the workunit to allocate the correct class. Has to be unique!
    integer :: client = -1      !< Defines the client on which the workunit will run
    integer :: oldClient = -1   !< Defines the client which was suspected to run the workunit, but has maybe changed by loadbalancing
    integer :: druid = -1       !< UID of datarequester


    logical :: balance = .false.  !< Indicates, if work unit can be balanced automatically
    logical :: isProcessed = .false.    !< Indicates if workunit has already been processed
    logical :: sendBack = .false.       !< If true, the result of the workunit will be send back to the scheduler
    type(intList) :: neededWUs          !< Dependencies

  contains

    procedure(iinit), deferred :: init
    procedure(iprocess), deferred :: process

    procedure(isetClient), deferred :: setClient
    procedure(isetOldClient), deferred :: setOldClient
    procedure(isetDrUID), deferred :: setDrUID

  end type workunit

  abstract interface

    subroutine iinit(this)
      import workunit
      class(workunit) :: this

    end subroutine

    subroutine isetDrUID(this, druid)
      import workunit
      class(workunit) :: this
      integer :: druid

    end subroutine isetDrUID

    subroutine isetOldClient(this, cl)
      import workunit

      class(workunit) :: this
      integer :: cl

    end subroutine isetOldClient

    subroutine isetClient(this, cl)
      import workunit
    
      class(workunit) :: this
      integer :: cl

    end subroutine isetClient

    subroutine iprocess(this)
      import workunit
    
      class(workunit) :: this

    end subroutine iprocess

  end interface

end module workunit_module
