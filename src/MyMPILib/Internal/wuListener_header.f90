  !> Child class of listener
  type, extends(listener) :: wuListener

   contains
     procedure :: init => init_wuListener
     procedure :: doListen => doListen_wuListener
  end type wuListener
