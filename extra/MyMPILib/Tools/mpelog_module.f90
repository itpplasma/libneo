!> Module for mpelog-class
module mpelog_module

  ! Include header file for MPE
#if defined(MPE_SUPPORT)
  include "mpe_logf.h"
#endif
  
  !> Class mpelog
  type :: mpelog
     logical :: active = .false.
     character(len = 8+4+1) :: startDateTime
   contains
     procedure :: init => init_mpelog
     procedure :: deinit => deinit_mpelog
     
     procedure :: logEvent => logEvent_mpelog
  end type mpelog
  
  integer :: mpe_e_sendA, mpe_e_sendB, mpe_e_recvA, mpe_e_recvB
  integer :: mpe_e_compA, mpe_e_compB, mpe_e_createA, mpe_e_createB
  integer :: mpe_e_busyA, mpe_e_busyB
  
  ! Singleton
  type(mpelog) :: mlog  
  
contains
  
  !> Init the mpe mechanism and create event numbers
  subroutine init_mpelog(this, myid)
    class(mpelog) :: this
    integer :: myid
    character(8) :: date
    character(10) :: time

#if defined(MPE_SUPPORT)
    if (this%active) then
       ierr = MPE_Init_Log()
       mpe_e_sendA = MPE_Log_get_event_number()
       mpe_e_sendB = MPE_Log_get_event_number()
       mpe_e_recvA = MPE_Log_get_event_number()
       mpe_e_recvB = MPE_Log_get_event_number()
       mpe_e_compA = MPE_Log_get_event_number()
       mpe_e_compB = MPE_Log_get_event_number()
       mpe_e_createA = MPE_Log_get_event_number()
       mpe_e_createB = MPE_Log_get_event_number()
       mpe_e_busyA = MPE_Log_get_event_number()
       mpe_e_busyB = MPE_Log_get_event_number()               
       
       if (myid == 0) then
          ierr = MPE_Describe_state(mpe_e_sendA, mpe_e_sendB, "Sending", "red")
          ierr = MPE_Describe_State(mpe_e_recvA, mpe_e_recvB, "Recvieving", "blue")
          ierr = MPE_Describe_State(mpe_e_compA, mpe_e_compB, "Computing", "green")
          ierr = MPE_Describe_State(mpe_e_createA, mpe_e_createB, "Creating", "lightgreen")
          ierr = MPE_Describe_State(mpe_e_busyA, mpe_e_busyB, "Scheduler", "grey")
       end if
       
       ierr = MPE_Start_log()

       call date_and_time(date, time)
       this%startDateTime = date // "_" // time(1:4)
    end if
#endif
  end subroutine init_mpelog
  
  !> Stop logging and write file
  subroutine deinit_mpelog(this)
    class(mpelog) :: this
#if defined(MPE_SUPPORT)
    if (this%active) then
       ierr = MPE_Finish_Log("mpeProfile_" // this%startDateTime)
    end if
#endif
  end subroutine deinit_mpelog
  
  !> Log special event
  subroutine logEvent_mpelog(this, ltype, opt1, opt2)
    class(mpelog) :: this
    integer, intent(in) :: ltype
    integer, optional :: opt1, opt2
    
    integer :: ierr
#if defined(MPE_SUPPORT)
    
    if (this%active) then
       ierr = MPE_Log_Event(ltype, 0, "")
       if (ltype == mpe_e_sendB) then
          ierr = MPE_Log_Send(opt1, opt2, 0)
       end if
       if (ltype == mpe_e_recvB) then
          ierr = MPE_Log_Receive(opt1, opt2, 0)
       end if
    end if
#endif
  end subroutine logEvent_mpelog
  
end module mpelog_module
