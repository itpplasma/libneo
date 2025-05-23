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

  use libneo_kinds, only : integer_kind, dp

  implicit none

  public :: solve_eqsys

  interface
  subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    integer, intent(in) :: n, nrhs, lda, ldb
    double precision, intent(in out) :: a(lda,*), b(ldb,*)
    integer, intent(out) :: ipiv(*), info
  end subroutine dgesv
  end interface

! --------------------------------------------------------------------
contains

  !> simple wrapper for solvers for real system of linear
  !> equations  A * X = B
  subroutine solve_eqsys(a, b, info)

    implicit none

    real(dp), dimension(:,:), intent(inout) :: a
    real(dp), dimension(:),   intent(inout) :: b
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
