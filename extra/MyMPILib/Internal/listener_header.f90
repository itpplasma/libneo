  !> Class listener
  type :: listener
     integer :: iTag    !< Defines the tag for listening
     class(scheduler), pointer :: myScheduler => null()   !< Defines the scheduler for callback functions
   contains
     procedure :: init => init_listener
     procedure :: deinit => deinit_listener

     procedure :: doListen => doListen_listener
  end type listener
