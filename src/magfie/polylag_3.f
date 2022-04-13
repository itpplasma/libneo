      module polylag_3
        integer, parameter :: mp=4
      contains
!
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
! the power 3 of polinomial is fixed strictly:
!
      implicit double precision (a-h,o-z)
!
      integer indu(mp)  
                             
      indu(1) = int((u-umin)*dum1)
      if( indu(1) .le. 0 ) indu(1) = 1
      indu(mp) = indu(1) + mp - 1
      if( indu(mp) .gt. nup ) then
         indu(mp) = nup
         indu(1) = indu(mp) - mp + 1
      endif
      do i=2,mp-1
         indu(i) = indu(i-1) + 1
      enddo

      return 
      end subroutine indef
!---------------------------------------------------------------------
      subroutine indsmp(index,nup,indu)
! defines interval for 1D interpolation on uniform mesh
! by known index.
! Normally looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    index - number of a cell on the mesh
!    nup - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points

! the power 3 of polinomial is fixed strictly:
      integer indu(mp)  
                             
      indu(1) = index - 1
      if( indu(1) .le. 0 ) indu(1) = 1
      indu(mp) = indu(1) + mp - 1
      if( indu(mp) .gt. nup ) then
         indu(mp) = nup
         indu(1) = indu(mp) - mp + 1
      endif
      do i=2,mp-1
         indu(i) = indu(i-1) + 1
      enddo

      return 
      end subroutine indsmp
!---------------------------------------------------------------------
      subroutine plag1d(x,fp,dxm1,xp,polyl1d,p1x1d)
! 1D interpolation by means of Lagrange polynomial
! the power 3 is fixed strictly:
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
!
      implicit double precision (a-h,o-z)
!
      dimension cx(mp),xp(mp),cx1(mp),fp(mp)
!
      call coefs(x,xp,dxm1,cx)
      polyl1d = 0.d0
      do i=1,mp
         polyl1d = polyl1d + fp(i)*cx(i)
      enddo
!
      call coefs1(x,xp,dxm1,cx1)
      p1x1d = 0.d0
      do i=1,mp
         p1x1d = p1x1d + fp(i)*cx1(i)
      enddo
!
      return
      end subroutine plag1d
!---------------------------------------------------------------------
      subroutine plag3d(x,y,z,fp,dxm1,dym1,dzm1,xp,yp,zp
     $     ,polyl3d,poly1x,poly1y,poly1z)
!
      implicit double precision (a-h,o-z)
!
! 3D interpolation by means of Lagrange polynomial
! the power 3 is fixed strictly:
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
      dimension cx(mp),cy(mp),cz(mp),fp(mp,mp,mp),xp(mp)
      dimension yp(mp),zp(mp),cx1(mp),cy1(mp),cz1(mp)
!
      call coefs(x,xp,dxm1,cx)
      call coefs(y,yp,dym1,cy)
      call coefs(z,zp,dzm1,cz)
!
      polyl3d = 0.d0
      do k=1,mp
         do j=1,mp
            do i=1,mp
               polyl3d = polyl3d + fp(i,j,k)*cx(i)*cy(j)*cz(k)
            enddo
         enddo
      enddo
!
      call coefs1(x,xp,dxm1,cx1)
      call coefs1(y,yp,dym1,cy1)
      call coefs1(z,zp,dzm1,cz1)
!
      poly1x = 0.d0
      do k=1,mp
         do j=1,mp
            do i=1,mp
               poly1x = poly1x + fp(i,j,k)*cx1(i)*cy(j)*cz(k)
            enddo
         enddo
      enddo
!
      poly1y = 0.d0
      do k=1,mp
         do j=1,mp
            do i=1,mp
               poly1y = poly1y + fp(i,j,k)*cx(i)*cy1(j)*cz(k)
            enddo
         enddo
      enddo
!
      poly1z = 0.d0
      do k=1,mp
         do j=1,mp
            do i=1,mp
               poly1z = poly1z + fp(i,j,k)*cx(i)*cy(j)*cz1(k)
            enddo
         enddo
      enddo
!
      return
      end subroutine plag3d
!---------------------------------------------------------------------
      subroutine coefs(u,up,dum1,cu)
!
      implicit double precision (a-h,o-z)
!
      dimension up(mp),cu(mp)
      data one6/0.16666666666667d0/
      du3 = dum1**3
      cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-one6*du3)
      cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5d0*du3)
      cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5d0*du3)
      cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (one6*du3)
      return
      end subroutine coefs
!---------------------------------------------------------------------
      subroutine coefs1(u,up,dum1,cu1)
!
      implicit double precision (a-h,o-z)
!
      dimension up(mp), cu1(mp)
      data one6/0.16666666666667d0/
      du3 = dum1**3
      cu1(1) = (  (u - up(3))*(u - up(4)) 
     $          + (u - up(2))*(u - up(4))
     $          + (u - up(2))*(u - up(3))
     $         )* (-one6*du3)
      cu1(2) = (  (u - up(3))*(u - up(4))
     $          + (u - up(1))*(u - up(4))
     $          + (u - up(1))*(u - up(3)) 
     $         ) * (0.5d0*du3)
      cu1(3) = (  (u - up(2))*(u - up(4))
     $          + (u - up(1))*(u - up(4))
     $          + (u - up(1))*(u - up(2))
     $         )* (-0.5d0*du3)
      cu1(4) = (  (u - up(2))*(u - up(3))
     $          + (u - up(1))*(u - up(3))
     $          + (u - up(1))*(u - up(2))
     $         )* (one6*du3)
      return
      end subroutine coefs1
!
      end module
