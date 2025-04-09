  !> Class definition for comListener
  type, extends(listener) :: comListener
  
   contains
       procedure :: init => init_comListener
       procedure :: doListen => doListen_comListener
  end type comListener
 
