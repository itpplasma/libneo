!***********************************************************************
!
! AUTHOR: Bernhard Seiwald
!
! DATE:   18.07.2001
!
!***********************************************************************

!> Contains simple wrapper for solvers for real system of linear
!> equations  A * X = B
!>
!> \todo Useful as own module? (shift into other module?)
module solve_systems

  use libneo_kinds, only : integer_kind, real_kind

  implicit none

  public :: solve_eqsys

! --------------------------------------------------------------------
contains

  !> simple wrapper for solvers for real system of linear
  !> equations  A * X = B
  subroutine solve_eqsys(a, b, info)

    implicit none

    real(real_kind), dimension(:,:), intent(inout) :: a
    real(real_kind), dimension(:),   intent(inout) :: b
    integer(integer_kind),             intent(out)   :: info
    integer(integer_kind) :: i_alloc
    integer(integer_kind) :: n, nrhs, lda, ldb
    integer(integer_kind), dimension(:), allocatable :: ipiv
    ! ------------------------------------------------------------------

    lda  = size(a,1)
    n    = size(a,2)
    ldb  = size(b,1)
    nrhs = 1

    allocate(ipiv(n),  stat = i_alloc)
    if (i_alloc /= 0) stop 'solve_eqsys: Allocation for array failed!'

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    info = 0

    deallocate(ipiv,  stat = i_alloc)
    if (i_alloc /= 0) stop 'solve_eqsys: Deallocation for array failed!'

  end subroutine solve_eqsys

end module solve_systems
