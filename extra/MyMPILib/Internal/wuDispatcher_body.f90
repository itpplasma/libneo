!> Constructor of wuDispatcher sets the tag for sending workunits
subroutine init_wuDispatcher(this)
  class(wuDispatcher) :: this
  this%iTag = wuTag
end subroutine init_wuDispatcher

!> Calls the send routine from mpiprovider
subroutine send_wuDispatcher(this, destRank)
  class(wuDispatcher) :: this
  integer :: destRank

  !write (*,*) "Sending buffer to ", destRank
  call mpro%packBuffer%sendTo(destRank, this%iTag)
      
end subroutine send_wuDispatcher

!> Calls the isend routine from mpiprovider
function isend_wuDispatcher(this, destRank) result(res)
  class(wuDispatcher) :: this
  integer :: destRank
  integer :: res
  	
  res = mpro%packBuffer%iSendTo(destRank, this%iTag)
end function isend_wuDispatcher
