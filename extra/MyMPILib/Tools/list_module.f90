!> Module for classes list and node
module list_module

    implicit none

    !> Class node, for a node of the list
    type :: node
        class(node), pointer :: next => null()
        class(node), pointer :: prev => null()
    contains
        procedure :: add => add_node
        procedure :: free => free_node

        procedure :: print => print_node
        procedure :: getNext => getNext_node
        procedure :: getPrev => getPrev_node

        procedure :: getUID => getUID_node
        procedure :: getSortIndex => getSortIndex_node
    end type node
  
    !> Class list for storing all packable-objects
    type list
        class(node), pointer :: first => null()
        class(node), pointer :: last  => null()
        class(node), pointer :: currentElement => null()
        integer :: count = 0

        !> If true, list will be sorted by sortIndex
        logical :: sortList = .false.
    contains
        procedure, private :: changeRoot => changeRoot_list

        procedure :: init => init_list
        procedure :: addnode => addnode_list
        procedure :: del => del_list
     
        procedure :: print => print_list
        procedure :: free => free_list
        procedure :: getCount => getCount_list

        procedure :: getnode => getnode_list
        procedure :: hasElement => hasElement_list

        procedure :: gotoNext => gotoNext_list
        procedure :: gotoPrev => gotoPrev_list
        procedure :: rewind   => rewind_list
    end type list

contains

    !> ---------- NODE -------------

    subroutine free_node(this)
        class(node) :: this

        write (*,*) "List: Call of free_node should not happen, because it is designed to be inherited"
    end subroutine free_node

    function getUID_node(this)
        class(node) :: this
        integer :: getUID_node

        getUID_node = -1
        write (*,*) "List: Call of getUID_node should not happen, because it is designed to be inherited"
    end function getUID_node

    subroutine print_node(this)
        class(node) :: this
        write (*,*) "List: Call of print_node should not happen, because it is designed to be inherited"
    end subroutine print_node
    
    function getSortIndex_node(this) result(res)
        class(node) :: this
        real :: res

        res = -1
        write (*,*) "List: Call of getSortIndex_node should not happen, because it is designed to be inherited"
    end function getSortIndex_node

    function getNext_node(this)
        class(node) :: this
        class(node), pointer :: getNext_node
    
        getNext_node => null()
        if (associated(this%next)) then
            getNext_node => this%next
        end if
    end function getNext_node

    function getPrev_node(this)
        class(node) :: this
        class(node), pointer :: getPrev_node

        getPrev_node => null()
        if (associated(this%prev)) then
            getPrev_node => this%prev
        end if
    end function getPrev_node

    subroutine add_node(this, nextNode)
        class(node) :: this
        class(node), target :: nextNode
    
        this%next => nextNode
    end subroutine add_node
  
    !-------- LIST ---------

    subroutine init_list(this)
        class(list) :: this

    end subroutine init_list
  
    subroutine rewind_list(this)
        class(list) :: this
    
        if (associated(this%first)) then
            this%currentElement => this%first
        else
            this%currentElement => null()
        end if
    end subroutine rewind_list

    subroutine forward_list(this)
        class(list) :: this

        if (associated(this%last)) then
            this%currentElement => this%last
        else
            this%currentElement => null()
        end if
    end subroutine forward_list
  
    subroutine gotoNext_list(this)
        class(list) :: this

        if (associated(this%currentElement)) then
            this%currentElement => this%currentElement%getNext()
        else
            this%currentElement => this%first
        end if
    
    end subroutine gotoNext_list
  
    subroutine gotoPrev_list(this)
        class(list) :: this
    
        if (associated(this%currentElement)) then
            this%currentElement => this%currentElement%getPrev()
            if (.not. associated(this%currentElement)) then
                this%currentElement => this%first
            end if
        else
            this%currentElement => this%first
        end if
    end subroutine gotoPrev_list

    function hasElement_list(this, index)
        class(list) :: this
        integer :: index
    
        logical :: hasElement_list
        class(node), pointer :: arrow
        logical :: found

        found = .false.

        if (associated(this%first)) then
            arrow => this%first
            do while (associated(arrow))
                if (arrow%getUID() == index) then
                    found = .true.
                    exit
                end if
                arrow => arrow%getNext()
            end do
       
        end if
    
        hasElement_list = found
    end function hasElement_list

    function getnode_list(this, index, throwExcp) result(res)
        class(list) :: this
        integer :: index
        logical, intent(in), optional :: throwExcp
    
        class(node), pointer :: arrow
        class(node), pointer :: res
        logical :: found, stopIfNotFound

        ! Decide if program should abort, if element is not in list
        if (.not. present(throwExcp)) then
            stopIfNotFound = .true.
        else
            stopIfNotFound = throwExcp
        end if
    
        found = .false.
        nullify(res)
        arrow => this%first
  
        do while (associated(arrow))

            if (arrow%getUID() == index) then
                found = .true.
                res => arrow
                exit
            end if
            arrow => arrow%getNext()
        end do

        if (.not. found) then
            if (stopIfNotFound) then
                write (*,*) "FATAL ERROR: Element ", index, "in list not found!"
                call this%print()
                stop
            else
                !write (*,*) "WARNING! Element ", index, " not found"
                nullify(res)
            end if
        end if
    
    end function getnode_list
    
    subroutine changeRoot_list(this)
        class(list) :: this

        this%first => getNext_node(this%first)!this%first%getNext()

    end subroutine changeRoot_list

    subroutine addnode_list(this, newNode)
        class(list) :: this
        class(node), target :: newNode
    
        class(node), pointer :: arrow, arrow2
        integer :: i
    
        if (associated(this%first)) then
            if (.not. this%sortList) then
                if (associated(this%last)) then
                    arrow => this%last
                else
                    write (*,*) "An error occurred, last element of list not defined"
                    stop
                end if

                !arrow => this%first
                !do i = 1, this%count-1
                !    arrow => arrow%getNext()
                !end do
                arrow%next => newNode
                newNode%prev => arrow
                this%last => newNode
            else
                arrow => this%first
                if (arrow%getSortIndex() > newNode%getSortIndex()) then
                    ! New Root
                    newNode%prev => null()
                    newNode%next => arrow
                    arrow%prev   => newNode
                    this%first   => newNode
                else
                    do while (associated(arrow))
                        if (arrow%getSortIndex() <= newNode%getSortIndex()) then
                            arrow2 => arrow
                            arrow => arrow%getNext()
                        else
                            arrow2 => arrow
                            arrow  => arrow%getPrev()
                   
                            newNode%prev => arrow
                            newNode%next => arrow2

                            arrow%next   => newNode
                            arrow2%prev  => newNode
                            exit
                        end if
                    end do
                    if (.not. associated(arrow)) then
                        !New element should be the last element
                        newNode%prev => arrow2
                        arrow2%next  => newNode
                    end if
                end if
            end if
            this%count = this%count + 1

        else
            this%first => newNode
            this%last  => newnode
            this%count = 1
        end if

    end subroutine addnode_list

    subroutine del_list(this, index, freeMem)
        class(list) :: this
        integer :: index
        class(node), pointer :: arrow, prev, next
        logical, optional :: freeMem

        logical :: dealloc

        ! Decide if container element should be freed or not, after deletion of node
        dealloc = .true.
        if (present(freeMem)) then
            dealloc = freeMem
        end if
    
        if (associated(this%first)) then
            arrow => this%first
            do while (associated(arrow))
                if (arrow%getUID() /= index) then
                    arrow => arrow%getNext()
                else
                    exit
                end if
            end do

            ! Decide behaviour of deletion of current element
            !if (arrow%getUID() == this%currentElement%getUID()) then
            !   this%currentElement
            !end if

            prev => arrow%getPrev()
            next => arrow%getNext()

            if (this%count > 1) then
                if (associated(prev) .and. associated(next)) then
                    prev%next => next
                    next%prev => prev
                else
                    if (associated(next)) then
                        next%prev => null()
                        this%first => next
                    else
                        prev%next => null()
                    end if

                    if (associated(prev)) then
                        prev%next => null()
                        this%last => prev
                    else
                        next%prev => null()
                    end if
                end if
            else
                this%first => null()
                this%last  => null()
               !list empty
            end if

            this%count = this%count - 1

            if (associated(this%currentElement)) then
                if (associated(this%currentElement%prev)) then
                    this%currentElement => this%currentElement%prev
                else
                    if (associated(this%first)) then
                        this%currentElement => this%first
                    else
                        this%currentElement => null()
                    end if
                end if
            end if
       
            if (dealloc) then
                call arrow%free()
            end if
            if (associated(arrow)) deallocate(arrow)

        end if
    end subroutine del_list

    subroutine print_list(this)
        class(list) :: this
        class(node), pointer :: arrow

        if (associated(this%first)) then
            arrow => this%first
    
            write (*,*) "Output of list:"
            do while (associated(arrow))

                call arrow%print()
                arrow => arrow%getNext()

            end do
        else
            write (*,*) "<Empty list>"
        end if

    end subroutine print_list

    subroutine free_list(this, freeMem)
        class(list) :: this
        logical, optional :: freeMem
        logical :: dealloc
        class(node), pointer :: arrow
        integer :: i, id

        dealloc = .true.
        if (present(freeMem)) then
            dealloc = freeMem
        end if

        do i = 1, this%count
            arrow => this%first
            id = arrow%getUID()
            call this%del(id, dealloc)
  
        end do
    
    end subroutine free_list

    function getCount_list(this)
        class(list) :: this
        integer :: getCount_list

        getCount_list = this%count
    end function getCount_list

end module list_module
