!> Module for class initWorkunit
module initWorkunit_module

  use mpiprovider_module
  use genericworkunit_module
  
  implicit none

  !> Used to inherit a specific initial workunit
  type, extends(genericworkunit) :: initWorkunit

   contains
     procedure :: init    => init_InitWorkunit
     procedure :: process => process_InitWorkunit
     procedure :: free    => free_InitWorkunit
     procedure :: print   => print_InitWorkunit
     procedure :: pack    => pack_InitWorkunit
     procedure :: unpack  => unpack_InitWorkunit
     procedure :: get     => get_InitWorkunit
  end type initWorkunit

contains

  !> Constructor
  subroutine init_InitWorkunit(this)
    class(initWorkunit) :: this
    
  end subroutine init_InitWorkunit

  !> This subroutine has to be overwritten for defining special init-processes
  subroutine process_initWorkunit(this)
    class(initWorkunit) :: this

  end subroutine process_initWorkunit
  
  !> Can be used for getting data from the init-workunit
  function get_InitWorkunit(this, uid) result(res)
    class(initWorkunit) :: this
    integer, intent(in) :: uid
    class(packable), pointer :: res
    
  end function get_InitWorkunit
  
  !> Destructor
  subroutine free_initWorkunit(this)
    class(initWorkunit) :: this

  end subroutine free_initWorkunit

  subroutine pack_initWorkunit(this)
    class(initWorkunit) :: this

    write (*,*) "InitWU: Pack here not needed"
  end subroutine pack_initWorkunit

  subroutine unpack_initWorkunit(this)
    class(initWorkunit) :: this

    write (*,*) "InitWU: Unpack here not needed"
  end subroutine unpack_initWorkunit

  subroutine print_initWorkunit(this)
    class(initWorkunit) :: this

    write(*,*) "InitWU: Print not yet implemented"
  end subroutine print_initWorkunit

end module initWorkunit_module
