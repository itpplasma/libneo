!> Class dispatcher, see body for documentation
  type :: dispatcher
     integer :: iTag    !< Defines the tag for sending
     integer :: iErr    !< Indicates errors
   contains
     procedure :: init => init_dispatcher
     procedure :: deinit => deinit_dispatcher
  end type dispatcher
