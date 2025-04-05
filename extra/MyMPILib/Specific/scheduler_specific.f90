! Specific part of the scheduler. This has to be overwritten for a specific problem.

!> Needs to be overwritten to allocate specific workunits
function allocateSpecific_scheduler(this, wuType) result(res)
    class(scheduler) :: this
    character(len=maxStrLen) :: wuType
    class(workunit), pointer :: res

    write (*,*) "ERROR, you have forgotten to set your work units in the function allocateSpecific ", wuType
    stop
end function allocateSpecific_scheduler

!> Needs to be overwritten when there is need for special mergeworkunits
function allocateSpecificMergeWU_scheduler(this) result(res)
    class(scheduler) :: this
    class(wuMergeWorkunit), pointer :: res

    write (*,*) "ERROR, you have forgotten to set your work units in the function allocateSpecificMergeWu "
    stop
end function allocateSpecificMergeWU_scheduler

!> Can be overwritten for special type of loadbalancing
subroutine loadBalancing_scheduler(this)
    class(scheduler) :: this
    integer :: k

    class(workunit), pointer :: selectWU

    ! If balancing is activated in config file
    if (this%balance) then

        ! Iterate through all clients
        do k = 1, mpro%getNumProcs()-1

            ! Check if the client has no workunit at the moment
            if (this%clientStats(k)%isReady) then

                ! Check if the client has nothing more to do
                if (this%clientStats(k)%isDone) then
                    !write (*,*) "Client", k, "has nothing more to do."
                    call mpro%storage%waitingWorkunits%rewind()

                    ! Search unprocessed workunit
                    do while (associated(mpro%storage%waitingWorkunits%currentElement))
                        ! Tell the fortran compiler that the list contains workunits :-)
                        selectWU => mpro%storage%waitingWorkunits%getCurrent()
                        if (selectWU%balance) then
                            ! Search next free workunit
                            if ((selectWU%client /= k) .and. (.not. this%clientStats(selectWU%client)%isReady)) then
                                write (*,*) "------- Balancing", selectWU%uid, "from", selectWU%client, "to", k, " -------"

                                ! Set client of the workunit to the client, which is done
                                call selectWU%setClient(k)
                                mpro%balanceCount = mpro%balanceCount + 1
                                exit
                            end if
                        end if

                        call mpro%storage%waitingWorkunits%gotoNext()
                    end do
                else
                    ! If a client is ready, but not done yet, set it to done
                    ! If it will be done the next iteration of scheduler, it will get a balanced workunit
                    this%clientStats(k)%isDone = .true.
                end if
            end if
        end do
    end if


end subroutine loadBalancing_scheduler

!> Check if workunits have to be created during runtime
subroutine checkIfClientDone_scheduler(this)
  class(scheduler) :: this
  integer :: k
  type(wuDataRequester), pointer :: dr
  class(wuMergeWorkunit), pointer :: wu
  class(wuMergeWorkunit), pointer :: p
  integer :: temp
  
  class(workunit), pointer :: selectWU => null()
  class(workunit), pointer :: selectWU2 => null()
  
  do k = 1, mpro%getRank() -1
     this%clientStats(k)%localMerge = .false.
  end do
  
  ! In the first iteration, local merges will be done, in the second iteration merges over different clients will be done
  ! So, local merge has always higher priority
  do k = 1, 2
     
     call mpro%storage%processedWorkunits%rewind()
     !if ((mpro%waitingWorkunits%getCount() == 0 .and. mpro%pendingWorkunits%getCount() == 0) .or. this%readyToMerge) then
     
     do while (associated(mpro%storage%processedWorkunits%currentElement))
        
        selectWU => mpro%storage%processedWorkunits%getCurrent()
        select type (q => selectWU)
        class is (wuMergeWorkunit)
           
           !Do start with the second element
           call mpro%storage%processedWorkunits%gotoNext()
           
           !Check if element was not already involved in a merge-process
           if ((.not. q%isMerged) .and. (.not. q%doNotMerge)) then
              nullify(p)
              
              ! Get needed object to merge
              selectWU2 => mpro%storage%processedWorkunits%get(q%leftNeighbor, .false.)
              if (associated(selectWU2)) then
                 
                 select type (selectWU2)
                 class is (wuMergeWorkunit)
                    p => selectWU2
                 end select
                 
              end if
              
              ! If found (object was already processed and is in processedWorkunits list)
              if (associated(p)) then
                 
                 !write (*,*) "I would merge ", p%uid, q%uid, "on client", p%client
                 
                 !Check if clients have nothing to do
                 if (this%clientStats(p%client)%isReady .and. this%clientStats(q%client)%isReady .and. &
                      .not. this%clientStats(p%client)%isBlocked .and. .not. this%clientStats(q%client)%isBlocked) then
                    
                    ! Local merge has first priority !
                    if ((k == 1) .and. (p%client == q%client)) then

                       this%clientStats(p%client)%localMerge = .true.
                       
                       
                       wu => this%allocateMergeWU()!("wuExternalJoin")
                       !allocate(wu)
                       call wu%init()
                       wu%client = p%client
                       wu%resultUID = wu%uid
                       !wu%fracIndex = 0
                       call wu%setNeighbors(p%leftNeighbor, q%rightNeighbor)
                       call wu%setMergeInfo(p%resultUID, q%resultUID)

                       write (*,*) "Local merge on client", p%client, p%uid, q%uid

                       if (p%druid /= -1) then
                          call wu%neededWUs%add(p%druid)
                       end if

                       if (q%druid /= -1) then
                          call wu%neededWUs%add(q%druid)
                       end if

                       call mpro%storage%waitingWorkunits%add(wu)

                       call q%setMerged(.true.)
                       call this%repairNeighbors(p%uid, wu%leftNeighbor, wu%uid, wu%rightNeighbor)
                       call mpro%storage%processedWorkunits%rewind()


                       !call mpro%storage%waitingWorkunits%print()
                    else
                       if ((k == 2)) then
                          if (q%oldClient == -1) then
                             nullify(dr)
                             allocate(dr)
                             call dr%init()
                             dr%client   = q%client
                             dr%dest     = p%client !Right to left
                             dr%whichUID = q%uid
                             dr%fracIndex = 0
                             call q%setOldClient(q%client)
                             !call dr%neededWUs%add(dr%whichUID)
                             call q%setClient(p%client)

                             !this%clientStats(dr%dest)%isReady = .false.
                             this%clientStats(dr%dest)%isBlocked = .true.
                             !write (*,*) "Blocking client", dr%dest

                             temp = dr%uid
                             call mpro%storage%waitingWorkunits%add(dr)
                             call p%setDrUID(temp)
                             call q%setDrUID(temp)

                             call myLog%logCreateDR(dr%uid, dr%client, dr%dest, dr%whichUID)
                             write (*,*) "Creating dataRequester", dr%uid,"for client", dr%client, "UID=", &
                                  dr%whichUID, "Dest=", dr%dest
                             !write (*,*) "Workunit status: ", q%doNotMerge, q%type, q%client

                          end if
                       end if
                    end if

                    !exit
                 else

                 end if !Local merge ?

              end if ! Client are ready?

           else
              cycle
           end if ! workunit already merged

           ! Skip other workunits
           class default
           call mpro%storage%processedWorkunits%gotoNext()

        end select !Get only mergeWorkunits

     end do !Loop over processedWorkunits

  end do ! loop over k

end subroutine checkIfClientDone_scheduler


!> Can be overwritten for special initial routines before starting the scheduling
subroutine initMaster_scheduler(this)
  class(scheduler) :: this

end subroutine initMaster_scheduler

!> Can be overwritten for routines after the scheduling
subroutine summarize_scheduler(this)
    class(scheduler) :: this

end subroutine summarize_scheduler

!> Linear mode of scheduler
subroutine initLinear_scheduler(this)
    class(scheduler) :: this

end subroutine initLinear_scheduler

subroutine loadSettings_scheduler(this)
    class(scheduler) :: this

end subroutine loadSettings_scheduler
