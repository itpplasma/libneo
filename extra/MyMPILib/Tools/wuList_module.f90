!> List of work units
module wuList_module

  use list_module
  use workunit_module
  
  implicit none
  
  type, extends(node) ::  workunitNode
     class(workunit), pointer :: container => null()
   contains
     procedure :: getUID => getUID_workunitNode
     procedure :: set => set_workunitNode
     procedure :: print => print_workunitNode

     procedure :: free => free_workunitNode

     procedure :: getSortIndex => getSortIndex_workunitNode
  end type workunitNode

  type, extends(list) :: workunitList
     !class(workunitnode), pointer :: first
     !class(workunitnode), pointer :: currentElement
   contains
     procedure :: add => add_workunitlist
     procedure :: get => get_workunitList
     procedure :: getCurrent => getCurrent_workunitList
     procedure :: relinkElementTo => relinkElementTo_workunitlist
     procedure :: moveListTo      => moveListTo_workunitlist
  end type workunitList
  
contains

  function getSortIndex_workunitNode(this) result(res)
    class(workunitNode) :: this
    real :: res

    if (associated(this%container)) then
        res = this%container%fracIndex
    else
        write (*,*) "ERROR in getSortIndex, Container of listnode is not associated"
        stop
    end if
  end function getSortIndex_workunitNode
  
  subroutine free_workunitNode(this)
    class(workunitNode) :: this

    if (associated(this%container)) then
        call this%container%free()
    end if
    if (associated(this%container)) deallocate(this%container)   !Not working with ifort
  end subroutine free_workunitNode
  
  function get_workunitList(this, index, throwExcp)
    class(workunitList) :: this
    integer :: index
    logical, optional :: throwExcp
    class(node), pointer :: element
    class(workunit), pointer :: get_workunitList
    
    get_workunitList => null()
    if (present(throwExcp)) then
       element => this%getnode(index, throwExcp)
    else
       element => this%getnode(index)
    end if
    if (associated(element)) then
       select type (q => element)
       class is (workunitnode)
          get_workunitList => q%container
       end select
    end if
  end function get_workunitList

  subroutine add_workunitList(this, data)
    class(workunitList) :: this
    class(workunit) :: data
    type(workunitnode), pointer :: newNode

    this%sortList = .true.
    allocate(newNode)
    call newNode%set(this%count + 1, data)
    call this%addNode(newNode)
  end subroutine add_workunitList
  
  subroutine relinkElementTo_workunitlist(this, newList, element)
    class(workunitlist) :: this
    type(workunitlist)  :: newList
    class(workunit) :: element
    
    !call element%print()
    call newList%add(element)
    call this%del(element%uid, .false.)
  end subroutine relinkElementTo_workunitlist

  subroutine moveListTo_workunitlist(this, newList)
    class(workunitlist) :: this
    type(workunitlist) :: newList
   
    class(node), pointer :: arrow

    arrow => this%first
    do while (associated(arrow))
       arrow => this%first
       select type (q => arrow)
       type is (workunitNode)
          call newList%add(q%container)
          call this%del(q%container%uid, .false.)           
       end select
    
    end do

    call this%free()
    !call newList%add(element)
    !call this%del(element%uid, .false.)
  end subroutine moveListTo_workunitlist
  
  function getCurrent_workunitList(this) result(res)
    class(workunitList) :: this
    class(workunit), pointer :: res

    !if (associated(this%currentElement)) then
    
       select type (q => this%currentElement)
       type is (workunitNode)
          if (associated(q%container)) then
             res => q%container
          else
             write (*,*) "FATAL ERROR: Container nil!"
             stop
          end if
       end select
       
  end function getCurrent_workunitList
  
  subroutine set_workunitnode(this, idx, val)
    class(workunitnode) :: this
    integer, intent(in) :: idx
    class(workunit), target :: val

    if (val%uid == -1) then
       val%uid = idx
       !write (*,*) "Element has no ID"
    end if
    this%container => val

  end subroutine set_workunitnode

  function getUID_workunitNode(this) result(res)
    class(workunitNode) :: this
    integer :: res

    res = -1
    if (associated(this%container)) then
       res = this%container%uid
    else
       write (*,*) "FATAL ERROR: Container not associated!"
    end if
  end function getUID_workunitNode

  subroutine print_workunitnode(this)
    class(workunitnode) :: this

    !write (*,*) this%container%uid
    call this%container%print()
  end subroutine print_workunitnode
  
end module wuList_module
