  !> Constructor of listener
  subroutine init_listener(this, mySched)
    class(listener) :: this
    class(scheduler), target :: mySched

  end subroutine init_listener

  !> Listen routine, should be called after MPI_PROBE
  subroutine doListen_listener(this, source,  tag)
    class(listener) :: this
    integer :: source, tag
    
  end subroutine doListen_listener
  
  !> Destructor
  subroutine deinit_listener(this)
    class(listener) :: this

  end subroutine deinit_listener
