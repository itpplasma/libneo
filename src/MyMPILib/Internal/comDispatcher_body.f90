!> Constructor
subroutine init_comDispatcher(this)
  class(comDispatcher) :: this

  ! Defines the tag for sending and listening
  this%iTag = comTag
end subroutine init_comDispatcher

!> Stops all clients from listening
subroutine sendTermSignal_comDispatcher(this, dest)
  class(comDispatcher) :: this
  integer, optional :: dest         !< Destination rank

  if (mpro%getNumProcs() > 1) then

    ! Decision if broadcasting or not
    if (present(dest)) then
       call this%send(dest, "TERM")
    else
       call this%sendBroadcast("TERM")
    end if

  end if
end subroutine sendTermSignal_comDispatcher

!> Starts the InitWU on the clients
subroutine sendInitSignal_comDispatcher(this, dest)
  class(comDispatcher) :: this
  integer, optional :: dest        !< Destination rank

  ! Decision if broadcasting or not
  if (present(dest)) then
     call this%send(dest, "INIT")
  else
     call this%sendBroadcast("INIT")
  end if
end subroutine sendInitSignal_comDispatcher

!> Sends the command to print the last processed workunit to one client
subroutine sendPrintSignal_comDispatcher(this, dest)
  class(comDispatcher) :: this
  integer, optional :: dest       !< Destination rank

  ! Decision if broadcasting or not  
  if (present(dest)) then
     call this%send(dest, "PRNT")
  else
     call this%sendBroadcast("PRNT")
  end if
end subroutine sendPrintSignal_comDispatcher

!> Generic function to send a command string to a client
subroutine sendSignal_comDispatcher(this, buf, dest)
  class(comDispatcher) :: this
  character(len=4) :: buf       !< Command string (max. length 4)
  integer, optional :: dest     !< Destination rank, leave empty to send a broadcast

  ! Decision if broadcasting or not
  if (present(dest)) then
     call this%send(dest, buf)
  else
     call this%sendBroadcast(buf)
  end if
end subroutine sendSignal_comDispatcher

!> Wrapper function for sending the command
subroutine send_comDispatcher(this, dest, command)
  class(comDispatcher) :: this
  integer :: dest               !< Destination rank
  character(len=4) :: command   !< Command string (max. length 4)

  ! Send command
  call mpro%send(dest, this%iTag, command)

end subroutine send_comDispatcher

!> Wrapper function for sending a broadcast command
subroutine sendBroadcast_comDispatcher(this, command)
  class(comDispatcher) :: this
  character(len=4) :: command   !< Command string
  integer :: i

  ! Broadcasting
  do i = 1, mpro%getNumProcs()-1
     call this%send(i, command)
  end do
end subroutine sendBroadcast_comDispatcher
