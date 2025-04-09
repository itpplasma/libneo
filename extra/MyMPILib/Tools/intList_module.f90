!> List of integer values
module intlist_module
  use list_module
  implicit none

  type, extends(node) :: intnode
     integer :: value
   contains
     procedure :: getUID => getUID_intnode
     procedure :: set    => set_intnode
     procedure :: print  => print_intnode
     procedure :: free   => free_intnode
  end type intnode

  type, extends(list) :: intlist

   contains
     procedure :: add => add_intlist
     procedure :: get => get_intlist
     procedure :: getCurrent => getCurrent_intlist
     procedure :: relinkElementTo => relinkElementTo_intlist
  end type intlist

contains

  function getUID_intnode(this) result(res)
    class(intnode) :: this
    integer :: res

    res = this%value
  end function getUID_intnode

  subroutine print_intnode(this)
    class(intnode) :: this
  
    write (*,*) this%value
  end subroutine print_intnode
  
  subroutine free_intnode(this)
    class(intnode) :: this

  end subroutine free_intnode

  subroutine set_intnode(this, val)
    class(intnode) :: this
    integer :: val

    this%value = val
  end subroutine set_intnode

  function get_intlist(this, index) result(res)
    class(intList) :: this
    integer :: index
    class(node), pointer :: element
    integer :: res

    res = 0
    element => this%getnode(index)
    
    if (associated(element)) then
       select type (q => element)
          class is (intnode)
          res = q%value
       end select
    end if
  end function get_intlist

  subroutine add_intlist(this, data)
    class(intlist) :: this
    integer :: data
    class(intnode), pointer:: newNode

    allocate(intnode ::newNode)
    call newNode%set(data)
    call this%addNode(newNode)
  end subroutine add_intlist

  subroutine relinkElementTo_intlist(this, newList, element)
    class(intlist) :: this
    type(intlist)  :: newList
    integer :: element

    !call element%print()
    call newList%add(element)
    call this%del(element, .false.)
  end subroutine relinkElementTo_intlist

  function getCurrent_intList(this) result(res)
    class(intList) :: this
    integer :: res
    
    select type (q => this%currentElement)
    type is (intnode)
       res = q%value
    class default
        write (*,*) "FATAL ERROR: Integer-List: Something else than an integer-node in the list!"
       res = -1
    end select
  end function getCurrent_intList
  
end module intlist_module
