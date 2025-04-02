!> Module for class wuMergeChunk
module wuMergeChunk_module

    use wuMergeWorkunit_module

    implicit none

    !> Child class of wuMergeWorkunit, for collection of workunits
    type, extends(wuMergeWorkunit) :: wuMergeChunk
      integer :: startUID
      integer :: rangeUID
      integer :: uidOffset
    contains

    end type wuMergeChunk
end module wuMergeChunk_module
