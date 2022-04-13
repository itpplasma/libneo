!> Module for class datarequester
module wuDataRequester_module
  use genericWorkunit_module
  use mpiprovider_module
  use list_module

  use initWorkunit_module
  use wuMergeChunk_module
  
    
  implicit none
  
  !> Class wuDataRequester, used for sending workunits between clients
  type, extends(genericWorkunit) :: wuDataRequester

    integer :: whichUID    !< Which workunit
    integer :: dest        !< Target client
 
  contains
    procedure :: process => process_wuDataRequester
    procedure :: print   => print_wuDataRequester
    procedure :: pack    => pack_wuDataRequester
    procedure :: unpack  => unpack_wuDataRequester
    procedure :: free    => free_wuDataRequester
    procedure :: init    => init_wuDataRequester
  end type wuDataRequester

contains
    
  !> Constructor sets the type and the next UID
  subroutine init_wuDataRequester(this)
    class(wuDataRequester) :: this

    this%type = "wuDataRequester"

    ! Call the scheduler of parent class
    call this%genericworkunit%init()
      
  end subroutine init_wuDataRequester

  subroutine print_wuDataRequester(this)
    class(wuDataRequester) :: this
      
    write (*,*) "DataRequester ", this%client, this%uid, this%whichUID, this%dest
  end subroutine print_wuDataRequester
    
  !> This routine packs and sends a workunit from one client to another
  subroutine process_wuDataRequester(this)
    class(wuDataRequester) :: this
    class(packable), pointer :: wu
    !integer, dimension(MPI_STATUS_SIZE) :: stat
    character(len=4) :: buffer
    double precision :: stime, etime, spacktime, epacktime

    stime = MPI_WTime()
    !write (*,*) "Requesting ", this%whichUID, this%dest

    wu => mpro%storage%waitingWorkunits%get(this%whichUID)
    !if (wu%type == 3) then
    select type (q => wu)
      class is (wuMergeChunk)
      this%whichUID = q%resultUID

      !call mpro%storage%waitingWorkunits%print()
      if (mpro%storage%waitingWorkunits%hasElement(this%whichUID)) then
        wu => mpro%storage%waitingWorkunits%get(this%whichUID)
      else
        write (*,*) "Error in dataRequester, element not found", this%whichUID
        stop
      end if
    end select
    !end if
    call mpro%packBuffer%clear()
    spacktime = MPI_WTime()
    call wu%pack()
    epacktime = MPI_WTime()
    !write (*,*) "DataRequester on", this%client, this%uid ,"sending", wu%uid, "to", this%dest

    call mpro%packBuffer%ssendTo(this%dest, wuTag)
       
    ! Another possibility:

    !!$       req = mpro%packBuffer%isendTo(this%dest, 111)
    !!$
    !!$       flag = .false.
    !!$       errormsg = .false.
    !!$       emergCounter = 0
    !!$       st = MPI_WTIME()
    !!$       do while (.not. flag)
    !!$          call MPI_TEST(req, flag, stat, ierr)
    !!$          emergCounter = emergCounter + 1
    !!$          !write (*,*) "waiting..."
    !!$          if (MPI_WTIME() - st > 10) then
    !!$             write (*,*) mpro%getRank(), "maybe in a deadlock"
    !!$            ! emergCounter = 0
    !!$          end if
    !!$       end do
    !!$       !write (*,*) "DONE.", mpro%getRank(), this%uid
    !!$
    !!$       call mlog%logEvent(mpe_e_recvA)
    !call MPI_WAIT(req, stat, ierr)
    !!$
    !!$       !write (*,*) "DataReq is waiting for", this%dest
    !!     call mpro%recv(this%dest, 77, buffer)
    !!$       !write (*,*) "OK"
    !!$       call mlog%logEvent(mpe_e_recvB, this%dest, 77)
    this%isProcessed = .true.

    etime = MPI_WTime()
    mpro%drcount = mpro%drcount + 1
    mpro%drTime = mpro%drTime + (etime - stime)
    mpro%drPackTime = mpro%drPackTime + (epacktime - spacktime)
  end subroutine process_wuDataRequester

  !> Unpack envelope of datarequester
  subroutine unpack_wuDataRequester(this)
    class(wuDataRequester) :: this

    call this%genericWorkunit%unpack()

    call mpro%packBuffer%get_int(this%whichUID)
    call mpro%packBuffer%get_int(this%dest)
  end subroutine unpack_wuDataRequester
    
  !> Pack envelope of datarequester
  subroutine pack_wuDataRequester(this)
    class(wuDataRequester) :: this
      
    call this%genericWorkunit%pack()

    call mpro%packBuffer%add_int(this%whichUID)
    call mpro%packBuffer%add_int(this%dest)
        
    if (this%isProcessed) then
         
    end if
  end subroutine pack_wuDataRequester

  !> Free memory
  subroutine free_wuDataRequester(this)
    class(wuDataRequester) :: this
      
    call this%genericWorkunit%free()
  end subroutine free_wuDataRequester
  
end module wuDataRequester_module
