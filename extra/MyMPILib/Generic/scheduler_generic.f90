! Body file of generic scheduler

  include "listener_body.f90"
  include "comListener_body.f90"
  include "wuListener_body.f90"

  include "dispatcher_body.f90"
  include "comDispatcher_body.f90"
  include "wuDispatcher_body.f90"

!> Constructor of generic scheduler
subroutine init_scheduler(this)
    class(scheduler) :: this
    integer :: stat
    integer :: f = 500
    character(len=100) :: msg

    mpro%schedInitTime = MPI_WTime()
    
    ! Namelist reading removed due to nvfortran 25.11 bug - using defaults
    write (*,*) "MyMPILib: Using default scheduler settings (namelist disabled for nvfortran compatibility)"
    ! Set values of nameList
    this%balance = loadBalancing
    myLog%verbose = verbose
    mlog%active = activateMPE

    call mpro%allocateBuffers(buffersize)

    this%isRunning = .false.

    ! Allocate Array for Clients
    allocate(this%clientStats(1:mpro%getNumProcs()-1))

    ! Allocate workunit lists
    allocate(mpro%storage%waitingWorkunits)
    allocate(mpro%storage%pendingWorkunits)
    allocate(mpro%storage%processedWorkunits)

    mpro%storage%processedWorkunits%sortList = .false.

    ! Initialize Listeners and Dispatchers
    call this%tComDisp%init()
    call this%tComListener%init(this)
    call this%tWuDisp%init()
    call this%tWuListener%init(this)

    if (mpro%getRank() .eq. 0 ) then
        !Initialize the logging on the master process
        call mylog%init()
    end if

    call mlog%init(mpro%getRank())

end subroutine init_scheduler

!> Wrapper for adding workunits to waitingWorkunits list
subroutine addWorkunit_scheduler(this, wu)
    class(scheduler) :: this
    class(workunit) :: wu

    call mpro%storage%waitingWorkunits%add(wu)

end subroutine addWorkunit_scheduler

!> Destructor of generic scheduler
subroutine deinit_scheduler(this)
    class(scheduler) :: this

    if (allocated(this%workunits_per_client)) deallocate(this%workunits_per_client)

    if (mpro%getRank() .eq. 0) then

        ! Performance output
        call this%printPerformanceAnalysis()

        ! Send TERM-Signal to all clients
        call this%tComDisp%sendTermSignal()

        ! Stop the logging
        call mylog%deinit()
    end if

    ! Stop MPE logger
    call mlog%deinit()

    ! Do cleaning up activities
    call this%cleanup()

    mpro%schedDeinitTime = MPI_WTime()
    
end subroutine deinit_scheduler

!> Prints a performance analysis section
subroutine printPerformanceAnalysis_scheduler(this)
    class(scheduler) :: this
    integer :: i
    write (*,*) "Scheduler: All jobs are done!"

    write (*,*) "----- PERFORMANCE ANALYSIS (Scheduler) -----"
    write (*,*) "Client    MeanTime"
    do i = 1, mpro%getNumProcs()-1
        write (*,'(I7, F12.3)') i, this%clientStats(i)%meanWorkTime
        call this%clientStats(i)%free()
    end do
    write (*,*) ""

    write (*,*) "----- PERFORMANCE ANALYSIS (Client) -----"
    write (*,*) "                                        Workunits   Time[s]   Total[s]    DataRequesters   Time[s]    Packtime[s]"

    ! Wait for flushing output
    call sleep(1)
    
    ! Print all processed work units
    ! call mpro%processedWorkunits%print()

end subroutine printPerformanceAnalysis_scheduler

!> Routine for the partition of a number of work units depending on the number of clients
subroutine partNearlyFair_scheduler(this, workunits_count)
    class(scheduler) :: this
    integer :: workunits_count

    integer :: client_count
    integer :: workunits_per_client
    integer :: workunits_offset
    integer :: workunits_offset_left
    integer :: current_client
    integer :: workunits_left_client

    client_count  = mpro%getNumProcs() - 1
    allocate(this%workunits_per_client(1:client_count))

    workunits_per_client        = workunits_count / client_count
    workunits_offset = mod(workunits_count, client_count)

    do current_client = 1, client_count
        workunits_left_client = workunits_per_client
        if (workunits_offset > 0) then
            workunits_left_client = workunits_left_client + 1
            workunits_offset = workunits_offset - 1
        end if

        this%workunits_per_client(current_client) = workunits_left_client
    end do
end subroutine partNearlyFair_scheduler

!> Internal routine for controlling the merge-process of merge-able objects, e.g., matrices
subroutine repairNeighbors_scheduler(this, p_uid, ln, wu_uid, rn)
    class(scheduler) :: this
    integer :: p_uid, ln, wu_uid, rn

    class(workunit), pointer :: selectWU => null()

    selectWU => mpro%storage%processedWorkunits%get(p_uid)
    select type (last => selectWU)
        class is (wuMergeWorkunit)

        call last%setMerged(.true.)
    end select


    ! Repair right Neighbors
    if (ln /= -1) then
        selectWU => mpro%storage%processedWorkunits%get(ln, .false.)
        if (associated(selectWU)) then
            select type (q1 => selectWU)
                class is (wuMergeWorkunit)
                call q1%setNeighbors(q1%leftNeighbor, wu_uid)
            end select
        else
            selectWU => mpro%storage%pendingWorkunits%get(ln, .false.)
            if (associated(selectWU)) then
                select type (q2 => selectWU)
                    class is (wuMergeWorkunit)
                    call q2%setNeighbors(q2%leftNeighbor, wu_uid)
                end select
            else
                selectWU => mpro%storage%waitingWorkunits%get(ln)
                if (associated(selectWU)) then
                    select type (q3 => selectWU)
                        class is (wuMergeWorkunit)
                        call q3%setNeighbors(q3%leftNeighbor, wu_uid)
                    end select
                end if
            end if
        end if
    end if

    ! Repair left Neighbors
    if (rn /= -1) then
        selectWU => mpro%storage%processedWorkunits%get(rn, .false.)
        if (associated(selectWU)) then
            select type (q1 => selectWU)
                class is (wuMergeWorkunit)
                call q1%setNeighbors(wu_uid, q1%rightNeighbor)
            end select
        else
            selectWU => mpro%storage%pendingWorkunits%get(rn, .false.)
            if (associated(selectWU)) then
                select type (q2 => selectWU)
                    class is (wuMergeWorkunit)
                    call q2%setNeighbors(wu_uid, q2%rightNeighbor)
                end select
            else
                selectWU => mpro%storage%waitingWorkunits%get(rn)
                if (associated(selectWU)) then
                    select type (q3 => selectWU)
                        class is (wuMergeWorkunit)
                        call q3%setNeighbors(wu_uid, q3%rightNeighbor)
                    end select
                end if
            end if
        end if
    end if
end subroutine repairNeighbors_scheduler

!> Function to allocate a work unit
function allocateWU_scheduler(this, wuType) result(res)
    class(scheduler) :: this
    character(len=maxStrLen) :: wuType
    class(workunit), pointer :: res

    nullify(res)
    select case (wutype)
        case ("wuDataRequester")
            allocate(wuDataRequester :: res)
        case default
            ! If no matching work unit has been found, then the specific routine of the inherited scheduler is called
            ! Therefore, the generic scheduler has to be inherited if new work units are introduced.
            res => this%allocateSpecific(wuType)
    end select

    if (.not. associated(res)) then
        write (*,*) "Recieved unknown workunit!! (type = ", wuType, "). Maybe you have forgotten to pack something?"
        stop
    end if

end function allocateWU_scheduler

!> Allocates merge-able work units
function allocateMergeWU_scheduler(this) result(res)
    class(scheduler) :: this
    character(len=maxStrLen) :: wuType
    class(wuMergeWorkunit), pointer :: res

    nullify(res)
    res => this%allocateSpecificMergeWu()

    if (.not. associated(res)) then
        write (*,*) "Received unknown work unit! Maybe you have forgotten to pack something?"
        stop
    end if

end function allocateMergeWU_scheduler

!> Routine to rebuild a received workunit and to process it
subroutine rebuildWU_scheduler(this)
    class(scheduler) :: this
    class(workunit), pointer :: wu
    integer :: uid
    character(len=maxStrLen) ::  wutype
    double precision :: stime, etime

    ! Check which Workunit
    call mpro%packBuffer%resetPos()
    call mpro%packBuffer%get_string(wutype)
    !write (*,*) mpro%getRank(), "The type is: ", type

    if (mpro%getRank() > 0) then

        ! Allocate object
        wu => this%allocateWU(wuType)
    
        ! Initialize and Unpack the received workunit
        call wu%init()

        call wu%unpack()
        ! Add workunit to list
        call mpro%storage%waitingWorkunits%add(wu)
        !call mpro%waitingWorkunits%print()
    
        stime = MPI_WTime()
        ! Decide if processing is needed
        if (.not. wu%isProcessed) then

            !write (*,*) mpro%getRank(), "Processing ", wu%uid
            call mlog%logEvent(mpe_e_compA)
            call wu%process()
            call mlog%logEvent(mpe_e_compB)
       
        else
           !write (*,*) mpro%getRank(), "Workunit ", wu%uid, "already processed, nothing to do!"
           !write (*,*) mpro%getRank(), "Received a processed workunit", wu%uid
        end if

        this%lastWuID = wu%uid

        ! Performance recordings
        etime = MPI_WTime()
        this%proTime = this%proTime + (etime - stime)
    
        mpro%wuTime = mpro%wuTime + (etime - stime)
        mpro%wuCount = mpro%wuCount + 1

    else

        ! Receive a sendBack work unit
        call mpro%packBuffer%get_int(uid)
        wu => mpro%storage%processedWorkunits%get(uid)
        call wu%unpack()

      !write (*,*) "Scheduler recieved sendBack-workunit"
    end if

end subroutine rebuildWU_scheduler

!> Runs the init workunit on the client
subroutine runInitWU_scheduler(this)
    class(scheduler) :: this
    double precision :: stime, etime

    ! Run the init-workunit on every client, when the INIT-command was received
    write (*,"(A, I2, A, A, A)") "This is ", mpro%getRank(), " on ", trim(mpro%getProcName()), " running the Init-WU!"
   
    if (associated(mpro%initWU)) then
        call mpro%initWU%init()
       
        stime = MPI_WTime()
        call mpro%initWU%process()
        etime = MPI_WTime()
       
        select type (q => mpro%initWU)
            class is (initWorkunit)

            class default
            write (*,*) "FATAL Error: No matching INIT-WU found!"
        end select
    else
       write (*,*) mpro%getRank(), "No Init-WU defined!, Skipping..."
    end if

end subroutine runInitWU_scheduler

!> Import routine which is called, when scheduler receives a Ready-signal from a client.
!> Re-links objects of the work unit - lists, sets clients states and more
subroutine setClientStatus_scheduler(this, source, status)
    class(scheduler) :: this
    integer :: source
    logical, optional :: status
    integer :: wuid
    class(workunit), pointer :: wu


    if (.not. present(status)) then

        ! Check, if client has pendingWorkunits
        if (this%clientStats(source)%lastWUs%getCount() > 0) then
       
            call this%clientStats(source)%lastWUs%rewind()
            wuid = this%clientStats(source)%lastWUs%getCurrent()

            wu => mpro%storage%pendingWorkunits%get(wuid)

            select type (q => wu)
                class is (wuDataRequester)
                if (this%clientStats(q%dest)%isBlocked) then
                    this%clientStats(q%dest)%isBlocked = .false.
                  !write (*,*) "Unblocking client", q%dest
                else
                    write (*,*) "FATAL Warning: Requesting to unblock Client ", q%dest, " while it is not blocked!"
                    write (*,*) "Additional information: DR (", q%uid, ")"
                    stop
                end if
                !this%clientStats(q%dest)%isReady = .true.
            end select

            ! Relink last pending workunit (FIFO)
            call mpro%storage%pendingWorkunits%relinkElementTo(mpro%storage%processedWorkunits, wu)
            call mylog%logRelinkPendToProc(wuid)

            ! Reduce spool count
            this%clientStats(source)%wuSpool = this%clientStats(source)%wuSpool - 1
            call this%clientStats(source)%lastWUs%del(wuid)
            call this%clientStats(source)%myWorkunits%add(wuid)
            this%clientStats(source)%workunitsLeft = this%clientStats(source)%workunitsLeft - 1
          
            this%clientStats(source)%isReady = .true.
            if (this%clientStats(source)%lastWUs%getCount() > 0) then
                this%clientStats(source)%isReady = .false.
            end if
           !write (*,*) "Setting status of", source, "to", this%clientStats(source)%isReady, wuid

        else
            if (mpro%getRank() == 0) then
                write (*,*) "FATAL ERROR: Client response not expected"
                write (*,*) "Information: Client ", source
                stop
            else
               !write (*,*) mpro%getRank(), "gets a response from", source
            end if
        end if
    else
        write (*,*) "WARNING: no stat"
        this%clientStats(source)%isReady = status
    end if
    
end subroutine setClientStatus_scheduler

!> Callback function for listener, when receiving a stop-signal
subroutine terminate_scheduler(this)
    class(scheduler) :: this
    
    this%isRunning = .false.
end subroutine terminate_scheduler
  
!> Initiates the init-workunit on the clients and waits for them to finish
subroutine initPhase_scheduler(this)
    class(scheduler) :: this
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: ierr, tag, i, source
    character(len=4) :: buffer
    
    if (associated(mpro%initWU)) then
        write (*,*) "Init-Phase"
        ! Send INIT-Signal to all clients to start the Init-Workunit
        call this%tComDisp%sendInitSignal()

        call mpro%initWU%init()

        ! Init-Phase - Wait for REDY-Signal from every client
        do i=1, mpro%getNumProcs()-1
            call MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            tag = status(MPI_TAG)
            source = status(MPI_SOURCE)
            !write (*,*) mpro%getRank(), "Received a message with tag = ", tag, "from ", status(MPI_SOURCE)

            call mpro%recv(source, tag, buffer)
            if (buffer == "REDY") this%clientStats(source)%isReady = .true.
        end do

        do i = 1, mpro%getNumProcs()-1
            if (.not. this%clientStats(source)%isReady) then
                write (*,*) "Not all clients responded, aborting."
                stop
            end if
        end do

    else
        write (*,*) "No Init-Workunit defined, skipping prepare()"
    end if
  
end subroutine initPhase_scheduler

!> Does some cleaning-ups after running the scheduling process and prints benchmark information.
subroutine cleanup_scheduler(this)
    class(scheduler) :: this
    integer :: i
	
    ! Deallocate array for clients
    deallocate(this%clientStats)

    mpro%storage%nextUID = 1
    this%lastWuId = -1

    call mpro%storage%waitingWorkunits%free()
    call mpro%storage%pendingWorkunits%free()
    call mpro%storage%processedWorkunits%free()

    deallocate(mpro%storage%waitingWorkunits)
    deallocate(mpro%storage%pendingWorkunits)
    deallocate(mpro%storage%processedWorkunits)
	  
end subroutine cleanup_scheduler

!> Needed for  scheduler Version 1
function getPendingNodesCount_scheduler(this) result(res)
    class(scheduler) :: this
    integer :: res
    integer :: i

    res = 0
    do i = 1, mpro%getNumProcs()-1
        res = res + this%clientStats(i)%pendingNodes%getCount()
    end do
    
end function getPendingNodesCount_scheduler

!> The core routine of the scheduler. Iterates over all waiting workunits and dispatches them to the clients
subroutine processWaitingWorkunits_scheduler(this)
    class(scheduler) :: this
    integer :: ierr, source, tag, i, uid, wuCount, req, cl, nextAutoClient
    real :: wuPercent, wuPercentBefore
    logical :: depCheck, probeFlag
    integer, dimension(MPI_STATUS_SIZE) :: status
 
    double precision :: time, stime, etime, timeDiff, sCalctime, eCalctime, sLoop, eLoop, tLoop, sLoop1, eLoop1, tLoop1
    character(8) :: date
    character(10) :: strtime
    character(len=8+4) :: res

    class(workunit), pointer :: currentWU => null()

    call date_and_time(date, strtime)

    write (*,*) "Scheduler: All workunits are prepared! Count = ", mpro%storage%waitingWorkunits%getCount()
    write (*,*) "Starting framework at ", date, " ", strtime(1:4)

    sCalctime = MPI_WTime()

    time = 0
    nextAutoClient = 1

    ! Set client neighbors for merging of Scheduler Version 1
    do i = 1, mpro%getNumProcs()-1
        this%clientStats(i)%nextNeighbor = i + 1
    end do
    this%clientStats(mpro%getNumProcs()-1)%nextNeighbor = 0
  
    wuCount = mpro%storage%waitingWorkunits%getCount()
    wuPercent = 0
    tloop = 0
    tloop1 = 0
    
    req = MPI_REQUEST_NULL
    
    ! --- Distribute workunits to clients ---
    ! Loop until every workunit has been send to a client
    do while (((mpro%storage%waitingWorkunits%getCount() > 0) .or. (mpro%storage%pendingWorkunits%getCount() > 0)) &
        .or. (this%getPendingNodesCount() > 1))
                     
        ! Rewind list-iterator to first element
        call mpro%storage%waitingWorkunits%rewind()
       
        call mlog%logEvent(mpe_e_busyA)

        ! Iterate through list until the last element is reached
        tloop = 0
        do while (associated(mpro%storage%waitingWorkunits%currentElement))

            currentWU => mpro%storage%waitingWorkunits%getCurrent()

            if (currentWU%client == -1) then    ! Client auto detection
                call currentWU%setClient(nextAutoClient)
                !write (*,*) "AutoClient: ", nextAutoClient, q%uid
                nextAutoClient = nextAutoClient + 1
                if (nextAutoClient >= mpro%getNumProcs()) nextAutoClient = 1
            end if

            !if (this%clientStats(q%client)%wuSpool < 1) then
            !	write (*,*) "Warning, Spool of ", q%client, " gets empty!", this%clientStats(q%client)%wuSpool
            !end if
             
            ! If the particular client is ready then prepare to sen
            if (this%clientStats(currentWU%client)%isReady .and. &
                (.not. this%clientStats(currentWU%client)%isBlocked)) then! .or. this%clientStats(q%client)%wuSpool < this%maxSpool) then
                
                ! Check dependencies
                sloop = MPI_WTime()
                depCheck = .true.
                call currentWU%neededWUs%rewind()
                do while (associated(currentWU%neededWUs%currentElement))
                    uid = currentWU%neededWUs%getCurrent()
                    if (.not. mpro%storage%processedWorkunits%hasElement(uid) .and. &
                        (.not. this%clientStats(currentWU%client)%lastWUs%hasElement(uid))) then
                        depCheck = .false.
                    end if
                   
                    call currentWU%neededWUs%gotoNext()
                end do
                
                eloop = MPI_WTime()
                tloop = tloop + (eloop - sloop)

                if (depCheck) then
                   
                    sloop1 = MPI_WTime()

                    ! Pack and Send
                    call MPI_WAIT(req, status, ierr)
                    call mpro%packBuffer%clear()

                    call currentWU%pack()
                   
                    call myLog%logSending(0, currentWU%client, currentWU%uid)

                    req = this%tWuDisp%isend(currentWU%client)
                                     
                    if (this%clientStats(currentWU%client)%lastTime > 0) then
                        timeDiff = MPI_WTime() - this%clientStats(currentWU%client)%lastTime
                        this%clientStats(currentWU%client)%meanWorkTime =  &
                            0.5*(this%clientStats(currentWU%client)%meanWorkTime + timeDiff)
                    end if
                   
                    this%clientStats(currentWU%client)%lastTime = MPI_WTime()
                    this%clientStats(currentWU%client)%isReady = .false.
                    this%clientStats(currentWU%client)%isDone  = .false.
                   
                    call this%clientStats(currentWU%client)%lastWUs%add(currentWU%uid)
                   
                    this%clientStats(currentWU%client)%wuSpool = this%clientStats(currentWU%client)%wuSpool + 1
                   
                    !write (*,*) "Scheduler sending WU", q%uid, " to ", q%client
                   
                    ! Relink element to pending-workunits - list
                    call myLog%logRelinkWaitToPend(currentWU%uid)
                    call mpro%storage%waitingWorkunits%relinkElementTo(mpro%storage%pendingWorkunits, currentWU)
                    cl = currentWU%client
                    call mpro%storage%waitingWorkunits%rewind()
                   
                    eloop1 = MPI_WTime()
                    tloop1 = tloop1 + (eloop1 - sloop1)

                else
                    !write (*,*) "Dependencies of workunit ",  q%uid, " for client ", q%client," not yet fullfilled!"
                    !call q%neededWUs%print()
                   
                    ! Iterate to next workunit
                    call mpro%storage%waitingWorkunits%gotoNext()
                end if !DepCheck
                
            else
                call mpro%storage%waitingWorkunits%gotoNext()
            end if  !Client ready?
          
        end do
              
        call MPI_WAIT(req, status, ierr)
        call mlog%logEvent(mpe_e_busyB)

        ! --- ---

        ! Percentage output
        wuPercentBefore = wuPercent;
        wuPercent = 100.0 - (real(mpro%storage%waitingWorkunits%getCount() + &
                                  mpro%storage%pendingWorkunits%getCount()) / wuCount)*100;
        if ((int(wuPercent / 10) /= int(wuPercentBefore / 10)) .and. (wuPercentBefore <= 99)) then
            if (wuPercent > wuPercentBefore) then
                write(*,"(A, A, A, F6.2, A, A, A)") achar(27), '[32m', ' ----- ', wuPercent, "% done -----", achar(27), '[0m'
            end if
        end if
       
        stime = MPI_WTime()

        ! ------ Listening --------
        if (mpro%storage%pendingWorkunits%getCount() > 0) then
            !--- Listening ---
          
            ! write (*,*) "Probing...", mpro%storage%pendingWorkunits%getCount(), mpro%storage%waitingWorkunits%getCount()
            !call mpro%storage%pendingWorkunits%print()
            call MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            !write (*,*) "Probing Done"
          
            probeFlag = .true.
            do while (probeFlag)
             
                call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, probeFlag, status, ierr)
                if (probeFlag) then
                    tag = status(MPI_TAG)
                    source = status(MPI_SOURCE)
                    call this%tComListener%doListen(source, tag)
                    call this%tWuListener%doListen(source, tag)
                end if
             
            end do
          
            ! Check, if one client is completely finished and can send his result to the next neighbor
            call this%checkIfClientDone()
            call this%loadBalancing()
          
           !--- End Listening ---
        else
            write (*,*) "The scheduler detected a possible deadlock! Aborting program!"
            call mpro%storage%waitingWorkunits%print()
            call mpro%storage%processedWorkunits%print()
            stop
        end if
       
        etime = MPI_WTime()
        time = time + (etime - stime)
       
    end do

    eCalctime = MPI_WTime()              
    write (*,"(A, F12.2, A)") " Scheduler needed  ", (eCalctime - sCalctime), " s for processing all work units."
    write (*,"(A, F12.2, A)") " Scheduler spent   ", time, " s for waiting for responses."
    write (*,"(A, I6)")      " Balanced work units: ", mpro%balanceCount

  !call this%tcomDisp%sendPrintSignal(1)
    
  !write (*,*) "Scheduling done, printing processed work units:"
  !call mpro%storage%processedWorkunits%print()

end subroutine processWaitingWorkunits_scheduler

!> Used for printing the last processed workunit on a client
subroutine printLast_scheduler(this)
    class(scheduler) :: this
    class(workunit), pointer :: wu
    
    write (*,*) "Printing last element:"
    wu => mpro%storage%waitingWorkunits%get(this%lastWuId)
    call wu%print()
    
end subroutine printLast_scheduler
  
!> This functions runs on a client, when starting the scheduling
subroutine initClient_scheduler(this)
    class(scheduler) :: this
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: stat
    integer :: tag, source
    double precision :: stime, etime, time
    integer :: countsOfRedy
    class(workunit), pointer :: wu
    logical :: sendBack
    integer :: i
    double precision :: wuTimeOverWuCount, drTimeOverDrCount, drPacktimeOverDrCount

    time = 0
    
    this%isRunning = .true.
    countsOfRedy = 1
  
    ! Run until TERM-Signal
    do while (this%isRunning)
        stime = MPI_WTime()
        ! Check for new messages
  
        call MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)
        tag = stat(MPI_TAG)
        source = stat(MPI_SOURCE)

        etime = MPI_WTime()
        time = time + (etime - stime)
        !write (*,*) mpro%getRank(), "Received a message with tag = ", tag, " from ", source

        ! Call listeners
        call this%tComListener%doListen(source, tag)
        call this%tWuListener%doListen(source, tag)

        sendBack = .false.
        if (tag == wuTag) then
            wu => mpro%storage%waitingWorkunits%get(this%lastWUid)
            sendBack = wu%sendBack
        end if

        ! Send REDY to scheduler
        if (source == 0) then
            if (this%isRunning) then
                !write (*,*) mpro%getRank(), " sends ready because of ", this%lastWUid
                if (.not. sendBack) then
                    call this%tComDisp%sendSignal("REDY", source)
                else
                    !write (*,*) "Client ", mpro%getRank(), " sending Workunit to scheduler"

                    call wu%pack()
                    call this%tWuDisp%send(0)
                end if
            end if
        end if
    end do

    ! Deinitialize everything
    !call mpro%initWU%free()
    call mpro%storage%waitingWorkunits%free()
    call this%tComListener%deinit()

    ! Synchronize timing information output
    do i = 1, mpro%getNumProcs()
       
       call MPI_BARRIER(mpro%mpi_comm_clients, ierr)
       if (i == mpro%getRank()) then

          wuTimeOverWuCount = 0
          drTimeOverDrCount = 0
          drPacktimeOverDrCount = 0
          
          if (mpro%wuCount .ne. 0) wuTimeOverWuCount = mpro%wuTime / mpro%wuCount
          if (mpro%drCount .ne. 0) then
             drTimeOverDrCount = mpro%drTime / mpro%drCount
             drPackTimeOverDrCount =  mpro%drPackTime / mpro%drCount
          end if
          
          write (*,"(A, I6, A4, A10, I13, F10.3, F11.3, I18, F10.3, F15.3)") & 
               " Details for node", mpro%getRank(), " on ", mpro%getProcName(), &
               mpro%wuCount, wuTimeOverWuCount, mpro%wuTime, &
               mpro%drCount,drTimeOverDrCount, drPacktimeOverDrCount

       end if
    end do
     
end subroutine initClient_scheduler

!> Initias the initPhase on the clients
subroutine prepare_scheduler(this)
    class(scheduler) :: this

    call this%loadSettings()

    if (mpro%getNumProcs() > 1) then
        ! Decide if master or client
        if (mpro%getRank() .eq. 0) then
            ! Call initWorkunit on every client
            call this%initPhase()
        end if
    end if

end subroutine prepare_scheduler

!> Starts the scheduling process on master and client
subroutine schedule_scheduler(this)
    class(scheduler) :: this

    if (mpro%getNumProcs() > 1) then
    
        ! Decide if master or client
        if (mpro%getRank() .eq. 0) then

            ! Prepare all workunits                -> This will be specific for the problem
            call this%initMaster()

            ! Process all workunits
            call this%processWaitingWorkunits()

            ! Call summarize-routine (can be overwritten by the user)
            call this%summarize()

        else
            ! Set client to ready-mode for recvieving workunits
            call this%initClient()
        end if
    else
        ! No parallel execution,
        ! Run in linear mode
        call this%initLinear()
    end if
    
end subroutine schedule_scheduler
