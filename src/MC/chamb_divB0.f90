subroutine chamb(y,phi,ierr)
  !
  ! checks whether the point is inside the vacuum chamber
  !  Input parameters:
  !            formal: y(i) - coordinates on the poloidal cross section
  !                    phi  - toroidal angle
  ! Output parameters:
  !            formal: ierr -error code (0 if the point is inside 1 - othervice)

  use libneo_kinds, only : dp

  implicit none


  real(dp), dimension(2), intent(in) :: y
  real(dp), intent(in) :: phi
  integer, intent(out) :: ierr

  ! dummy routine here. TODO: implement
  associate(dummy => y)
  end associate
  associate(dummy => phi)
  end associate
  ierr = 0
end subroutine chamb
