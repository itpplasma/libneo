!> Constructor
subroutine init_comListener(this, mySched)
  class(comListener) :: this
  class(scheduler), target :: mySched      !< Pointer to the parent scheduler

  !Define scheduler for callbacks
  this%myScheduler => mySched
  
  this%iTag = comTag
end subroutine init_comListener

!> Execute listening process, function checks if tag is allowed for this listener
subroutine doListen_comListener(this, source, tag)
  class(comListener) :: this
  integer :: tag          !< Listening tag
  integer :: source       !< Source rank
  character(len=4) :: buffer

  ! Only receive, if tags are corresponding
  if (tag == this%iTag) then

     ! Wrapper for MPI receive command
     call mpro%recv(source, tag, buffer)

     if (buffer == "TERM") then
        call terminate_scheduler(this%myScheduler)
        !call this%myScheduler%terminate() ! Works since gfortran-4.7.0
     end if
     if (buffer == "INIT") then
        call runInitWU_scheduler(this%myScheduler)
        !call this%myScheduler%runInitWU()
     end if
     if (buffer == "REDY") then
        call setClientStatus_scheduler(this%myScheduler, source)
        !call this%myScheduler%setClientStatus(source)
     end if
     if (buffer == "PRNT") then
        call printLast_scheduler(this%myScheduler)
        !call this%myScheduler%printLast()
     end if
 
  end if

end subroutine doListen_comListener
