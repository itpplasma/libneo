!> Constructor sets the tag for listining for workunits
subroutine init_wuListener(this, mySched)
  class(wuListener) :: this
  class(scheduler), target :: mySched

  this%myScheduler => mySched
  this%iTag = wuTag
end subroutine init_wuListener

!> Listens for workunits and calls the rebuildWU-routine fo scheduler
subroutine doListen_wuListener(this, source, tag)
  class(wuListener) :: this
  integer :: tag, source

  if (tag == this%iTag) then

     call mpro%packBuffer%receiveFrom(source, this%iTag)
     if (mpro%getRank() > 0) then
       call rebuildWU_scheduler(this%myScheduler)
       !call this%myScheduler%rebuildWU()!rebuildWU_scheduler(this%myScheduler)
     else

       call setClientStatus_scheduler(this%myScheduler, source)
       call rebuildWU_scheduler(this%myScheduler)

       !call this%myScheduler%setClientStatus(source)
       !call this%myScheduler%rebuildWU()!rebuildWU_scheduler(this%myScheduler)
     end if
  end if

end subroutine doListen_wuListener
