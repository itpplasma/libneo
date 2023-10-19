subroutine chamb(y,phi,ierr)
  !
  ! checks whether the point is inside the vacuum chamber
  !  Input parameters:
  !            formal: y(i) - coordinates on the poloidal cross section
  !                    phi  - toroidal angle
  ! Outout parameters:
  !            formal: ierr -error code (0 if the point is inside 1 - othervice)

  use libneo_kinds, only : real_kind

  implicit none


  integer :: ierr,ir,ip,iz
  real(kind=real_kind) :: phi
  real(kind=real_kind), dimension(2) :: y

  ierr = 0  ! dummy routine here. TODO: implement
end subroutine chamb
