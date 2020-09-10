!> Module for writing log-files
module myLog_module
  use commandline_parser_module

  implicit none
  
  type myLogClass
     
    integer :: verbose = 0
  contains

    procedure :: init => init_mylogClass
    procedure :: deinit => deinit_myLogClass

    procedure :: logRelinkWaitToPend => logRelinkWaitToPend_mylogClass
    procedure :: logRelinkPendToProc => logRelinkPendToProc_mylogClass
    procedure :: logSending          => logSending_mylogClass
    procedure :: logCreateDR         => logCreateDR_myLogClass

    procedure, private :: getDateTime => getDateTime_myLogClass
       
  end type myLogClass

  ! Singleton
  type(myLogClass) :: myLog
  
  
contains

  function getDateTime_myLogClass(this) result(res)
    class(myLogClass) :: this
    character(8) :: date
    character(10) :: time
    character(len=8+4+7) :: res

    call date_and_time(date, time)
    res = date(7:8) // "." // date(5:6) // "." // date(1:4) // " "&
    // time(1:2) // ":" // time(3:4) // ":" // time(5:6)
  end function getDateTime_myLogClass

  subroutine init_mylogClass(this)
    class(myLogClass) :: this
    integer :: stat
    character(8) :: date
    character(10) :: time
    call date_and_time(date, time)

    !this%verbose = comlineParser%getInt("-verbose=", 0)

    if (this%verbose > 0) then
      open(10, file="myLog_" // date // "_" // time(1:4) // ".txt" , action='write', iostat=stat)

      write (10, *) "--- Begin logfile at", this%getDateTime(),  "---"

    end if
  end subroutine init_mylogClass
    
  subroutine deinit_myLogClass(this)
    class(myLogClass) :: this

    if (this%verbose > 0) then
      write (10,*) "--- End logfile at ", this%getDateTime(), "---"
      close(10)
    end if
  end subroutine deinit_myLogClass

  subroutine logSending_myLogClass(this, source, dest, uid)
    class(mylogClass) :: this
    integer :: source, dest, uid

    if (this%verbose > 0) then
      write (10, *) this%getDateTime(), ": Sending ", uid, " from ", source, " to ", dest
    end if
  end subroutine logSending_myLogClass

  subroutine logRelinkWaitToPend_mylogClass(this, uid)
    class(mylogClass) :: this
    integer :: uid

    if (this%verbose > 0) write (10, *) this%getDateTime(), ": Scheduler linking", uid, "to pending list"
  end subroutine logRelinkWaitToPend_mylogClass

  subroutine logRelinkPendToProc_mylogClass(this, uid)
    class(mylogClass) :: this
    integer :: uid
      
    if (this%verbose > 0) write (10, *) this%getDateTime(), ": Scheduler linking", uid, "to processed list"
  end subroutine logRelinkPendToProc_mylogClass
    
  subroutine logCreateDR_myLogClass(this, uid, source, dest, whichUID)
    class(myLogClass) :: this
    integer :: uid, source, dest, whichUID

    if (this%verbose > 0) write (10, *) this%getDateTime(), ": Create DataRequester (", uid, &
    ") from client", source, " to client ", dest, " for WU " , whichUID

  end subroutine logCreateDR_myLogClass


end module myLog_module
