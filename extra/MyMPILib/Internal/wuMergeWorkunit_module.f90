!> Module for class wuMergeWorkunit
module wuMergeWorkunit_module

  use genericworkunit_module
  use mpiprovider_module
  
  implicit none

  !> This workunit is special for mergeWorkunits, which means that two workunits get merged to one
  !> For example matrix multiplication
  type, extends(genericWorkunit) :: wuMergeWorkunit
    integer :: idxM1, idxM2 = -1
    integer :: sourceM1, sourceM2 = 1 !0... initWU, 1... processedList

    integer :: leftNeighbor = -1
    integer :: rightNeighbor = -1

    integer :: resultUID = -1

    logical :: isMerged = .false.
    logical :: doNotMerge = .false.

  contains
    procedure :: process => process_wuMergeWorkunit
    procedure :: print   => print_wuMergeWorkunit
    procedure :: pack    => pack_wuMergeWorkunit
    procedure :: unpack  => unpack_wuMergeWorkunit
    procedure :: free    => free_wuMergeWorkunit
    procedure :: init    => init_wuMergeWorkunit

    procedure :: setNeighbors => setNeighbors_wuMergeWorkunit
    procedure :: setMerged    => setMerged_wuMergeWorkunit
    procedure :: setMergeInfo => setMergeInfo_wuMergeWorkunit
  end type wuMergeWorkunit

contains

  subroutine init_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this
    this%uid = mpro%storage%nextUID
    this%resultuid = this%uid
    mpro%storage%nextUId = mpro%storage%nextUID + 1
    !this%type = 1
    call this%neededWUs%init()
      
  end subroutine init_wuMergeWorkunit

  subroutine setMergeInfo_wuMergeWorkunit(this, left_uid, right_uid)
    class(wuMergeWorkunit) :: this
    integer :: left_uid, right_uid

  end subroutine setMergeInfo_wuMergeWorkunit

  subroutine setMerged_wuMergeWorkunit(this, m)
    class(wuMergeWorkunit) :: this
    logical, intent(in) :: m
    
    this%isMerged = m
  end subroutine setMerged_wuMergeWorkunit
    
  subroutine setNeighbors_wuMergeWorkunit(this, ln, rn)
    class(wuMergeWorkunit) :: this
    integer, intent(in) :: ln, rn
    !write (*,*) "This is MergeWorkunit", this%uid, "and someone sets my neighbors to", ln, rn
    this%leftNeighbor = ln
    this%rightNeighbor = rn
  end subroutine setNeighbors_wuMergeWorkunit

  subroutine print_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this
      
    write (*,*) "Warning, a call of this print()-method should not happen!"
  end subroutine print_wuMergeWorkunit
    
  subroutine process_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this
      
  end subroutine process_wuMergeWorkunit

  subroutine unpack_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this

  end subroutine unpack_wuMergeWorkunit
    
  subroutine pack_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this

  end subroutine pack_wuMergeWorkunit

  subroutine free_wuMergeWorkunit(this)
    class(wuMergeWorkunit) :: this
      
    call this%neededWUs%free()
  end subroutine free_wuMergeWorkunit
    

end module wuMergeWorkunit_module
