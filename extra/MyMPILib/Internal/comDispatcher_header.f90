  !> Class definition for comDispatcher, see body file for documentation of the functions
  type, extends(dispatcher) :: comDispatcher
     
   contains
     procedure :: init => init_comDispatcher
     procedure, private :: send => send_comDispatcher
     procedure, private :: sendBroadcast => sendBroadcast_comDispatcher

     procedure :: sendTermSignal => sendTermSignal_comDispatcher
     procedure :: sendInitSignal => sendInitSignal_comDispatcher
     procedure :: sendPrintSignal=> sendPrintSignal_comDispatcher
     procedure :: sendSignal => sendSignal_comDispatcher
  end type comDispatcher
