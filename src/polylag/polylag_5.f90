module neo_polylag_5
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

integer, parameter :: mp=6

contains

subroutine indef(u,umin,dum1,nup,indu)
! defines interval for 1D interpolation on uniform mesh, normally 
! looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    u - coordinate of a point (to be interpolated)
!    umin - minimal value of u
!    dum1 = 1./h reciprocal to length of mesh interval
!    nup - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points
!
! the power 5 of polinomial is fixed strictly:
real(dp), intent(in) :: u, umin, dum1
integer, intent(in) :: nup 
integer, dimension(mp), intent(out) :: indu

integer :: i
                        
indu(1) = int((u-umin)*dum1)+1
if( indu(1) .le. 0 ) indu(1) = 1
indu(mp) = indu(1) + mp - 1
if( indu(mp) .gt. nup ) then
    indu(mp) = nup
    indu(1) = indu(mp) - mp + 1
endif
do i=2,mp-1
    indu(i) = indu(i-1) + 1
enddo
end subroutine indef


subroutine plag1d(x,fp,dxm1,xp,polyl1d,p1x1d)
! 1D interpolation by means of Lagrange polynomial
! the power 5 is fixed strictly:
! uniform mesh (increasingly ordered) in all dimensions is implied
!
! Input parameters:
!   x - coordinate of the point for interpolation
!   dxm1 - 1./h_x
!   xp - vertices of stencil
!
! Output parameters:
! polyl1d - polynomial itself
! poly1x - its derivative
real(dp), intent(in) :: x, dxm1
real(dp), dimension(mp), intent(in) :: xp, fp
real(dp), intent(out) :: polyl1d, p1x1d

real(dp), dimension(mp) :: cx, cx1
integer :: i

call coefs(x,xp,dxm1,cx)
polyl1d = 0.d0
do i=1,mp
    polyl1d = polyl1d + fp(i)*cx(i)
enddo

call coefs1(x,xp,dxm1,cx1)
p1x1d = 0.d0
do i=1,mp
    p1x1d = p1x1d + fp(i)*cx1(i)
enddo
end subroutine plag1d


subroutine plag3d(x,y,z,fp,dxm1,dym1,dzm1,xp,yp,zp, &
                    polyl3d,poly1x,poly1y,poly1z)
! 3D interpolation by means of Lagrange polynomial
! the power 5 is fixed strictly:
! uniform mesh (increasingly ordered) in all dimensions is implied
!
! Input parameters:
!   x,y,z - coordinates of the point for interpolation
!   dxm1,dym1,dzm1 - steps in each direction
!   xp,yp,zp - vertices of stencil
!
! Output parameters:
! polyl3d - polynomial itself
! poly1x - its x-derivative
! poly1y - its y-derivative
! poly1z - its z-derivative
real(dp), intent(in) :: x, y, z, dxm1, dym1, dzm1
real(dp), dimension(mp,mp,mp), intent(in) :: fp
real(dp), dimension(mp), intent(in) :: xp, yp, zp
real(dp), intent(out) :: polyl3d, poly1x, poly1y, poly1z

real(dp), dimension(mp) :: cx, cy, cz, cx1, cy1, cz1
integer :: i, j, k

call coefs(x,xp,dxm1,cx)
call coefs(y,yp,dym1,cy)
call coefs(z,zp,dzm1,cz)

polyl3d = 0.d0
do k=1,mp
    do j=1,mp
    do i=1,mp
        polyl3d = polyl3d + fp(i,j,k)*cx(i)*cy(j)*cz(k)
    enddo
    enddo
enddo

call coefs1(x,xp,dxm1,cx1)
call coefs1(y,yp,dym1,cy1)
call coefs1(z,zp,dzm1,cz1)

poly1x = 0.d0
do k=1,mp
    do j=1,mp
    do i=1,mp
        poly1x = poly1x + fp(i,j,k)*cx1(i)*cy(j)*cz(k)
    enddo
    enddo
enddo

poly1y = 0.d0
do k=1,mp
    do j=1,mp
    do i=1,mp
        poly1y = poly1y + fp(i,j,k)*cx(i)*cy1(j)*cz(k)
    enddo
    enddo
enddo

poly1z = 0.d0
do k=1,mp
    do j=1,mp
    do i=1,mp
        poly1z = poly1z + fp(i,j,k)*cx(i)*cy(j)*cz1(k)
    enddo
    enddo
enddo
end subroutine plag3d


subroutine coefs(u,up,dum1,cu)
real(dp), intent(in) :: u, dum1
real(dp), dimension(mp), intent(in) :: up

real(dp) :: du5
real(dp), dimension(mp) :: cu

du5 = dum1**5
cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * &
        (u - up(5)) * (u - up(6)) * (-du5)/120.d0
cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * &
        (u - up(5)) * (u - up(6)) * (du5)/24.d0
cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * &
        (u - up(5)) * (u - up(6)) * (-du5)/12.d0
cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
        (u - up(5)) * (u - up(6)) * (du5)/12.d0
cu(5) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
        (u - up(4)) * (u - up(6)) * (-du5)/24.d0
cu(6) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
        (u - up(4)) * (u - up(5)) * (du5)/120.d0
end subroutine coefs


subroutine coefs1(u,up,dum1,cu1)
real(dp), intent(in) :: u, dum1
real(dp), dimension(mp), intent(in) :: up

real(dp) :: du5
real(dp), dimension(mp) :: cu1

du5 = dum1**5
cu1(1) = ((u - up(3)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(2)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(2)) * (u - up(3)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(2)) * (u - up(3)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(2)) * (u - up(3)) * &
          (u - up(4)) * (u - up(5)))*(-du5)/120.d0
cu1(2) = ((u - up(3)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(4)) * (u - up(5)))*(du5)/24.d0
cu1(3) = ((u - up(2)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(4)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(4)) * (u - up(5)))*(-du5)/12.d0
cu1(4) = ((u - up(2)) * (u - up(3)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(5)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(5)))*(du5)/12.d0
cu1(5) = ((u - up(2)) * (u - up(3)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(4)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(6)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(4)))*(-du5)/24.d0
cu1(6) = ((u - up(2)) * (u - up(3)) * &
          (u - up(4)) * (u - up(5)) + &
          (u - up(1)) * (u - up(3)) * &
          (u - up(4)) * (u - up(5)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(4)) * (u - up(5)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(5)) + &
          (u - up(1)) * (u - up(2)) * &
          (u - up(3)) * (u - up(4)))*(du5)/120.d0
end subroutine coefs1

end module neo_polylag_5
