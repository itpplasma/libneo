!> List of packable objects
module packableList_module

  use list_module
  use packable_module
  
  implicit none
  
  type, extends(node) ::  packableNode
     class(packable), pointer :: container => null()
   contains
     procedure :: getUID => getUID_packableNode
     procedure :: set => set_packableNode
     procedure :: print => print_packableNode

     procedure :: free => free_packableNode
  end type packableNode

  type, extends(list) :: packableList
     !class(packablenode), pointer :: first
     !class(packablenode), pointer :: currentElement
   contains
     procedure :: add => add_packablelist
     procedure :: get => get_packableList
   
     procedure :: getCurrent => getCurrent_packableList
     procedure :: relinkElementTo => relinkElementTo_packablelist
  end type packableList
  
contains

  subroutine free_packableNode(this)
    class(packableNode) :: this

    call this%container%free()
    !if (associated(this%container)) deallocate(this%container)   !Not working with ifort
  end subroutine free_packableNode
  
  function get_packableList(this, index)
    class(packableList) :: this
    integer :: index
    class(node), pointer :: element
    class(packable), pointer :: get_packableList
    
    get_packableList => null()
    element => this%getnode(index)
    if (associated(element)) then
       select type (q => element)
       class is (packablenode)
          get_packableList => q%container
       end select
    end if
  end function get_packableList

  subroutine add_packableList(this, data)
    class(packableList) :: this
    class(packable) :: data

    class(packablenode), pointer :: newNode
    
    allocate(packableNode::newNode)
    call newNode%set(this%count + 1, data)
    call this%addNode(newNode)
  end subroutine add_packableList
  
  subroutine relinkElementTo_packablelist(this, newList, element)
    class(packablelist) :: this
    type(packablelist)  :: newList
    class(packable) :: element
    
    !call element%print()
    call newList%add(element)
    call this%del(element%uid, .false.)
  end subroutine relinkElementTo_packablelist
  
  function getCurrent_packableList(this) result(res)
    class(packableList) :: this
    class(packable), pointer :: res
    
    select type (q => this%currentElement)
    type is (packableNode)
       res => q%container
    end select
  end function getCurrent_packableList
  
  subroutine set_packablenode(this, idx, val)
    class(packablenode) :: this
    integer, intent(in) :: idx
    class(packable), target :: val

    if (val%uid == -1) then
       val%uid = idx
       !write (*,*) "Element has no ID"
    end if
    this%container => val

  end subroutine set_packablenode

  function getUID_packableNode(this)
    class(packableNode) :: this
    integer :: getUID_packableNode

    getUID_packableNode = this%container%uid
  end function getUID_packableNode

  subroutine print_packablenode(this)
    class(packablenode) :: this

    !write (*,*) this%container%uid
    call this%container%print()
  end subroutine print_packablenode
  
end module packableList_module
