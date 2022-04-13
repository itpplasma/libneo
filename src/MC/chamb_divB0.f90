subroutine chamb(y,phi,ierr)
    !
    ! checks whether the point is inside the vacuum chamber
    !  Input parameters:
    !            formal: y(i) - coordinates on the poloidal cross section
    !                    phi  - toroidal angle
    ! Outout parameters:
    !            formal: ierr -error code (0 if the point is inside 1 - othervice)
    !
    !
          implicit none
    !
    !
          integer :: ierr,ir,ip,iz
          double precision :: phi
          double precision, dimension(2) :: y

          ierr = 0  ! dummy routine here. TODO: implement
end subroutine chamb
