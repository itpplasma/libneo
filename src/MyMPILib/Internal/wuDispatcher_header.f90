!> Child class of dispatcher
type, extends(dispatcher) :: wuDispatcher

   contains
     procedure :: init  => init_wuDispatcher
     procedure :: send  => send_wuDispatcher
     procedure :: isend => isend_wuDispatcher
   
end type wuDispatcher
