!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine s2dcut(nx,ny,hx,hy,f,imi,ima,jmi,jma,icount,spl,ipoint)
  ! Calculates coefficients of a 2D spline for a convex domain
  ! (for a non-rectangular domain the interpolation is not continious)
  ! equidistant mesh, but hx must not be = hy
  !
  !  Input parameters:
  !                    nx           - horizontal size of the mesh (over x)
  !                    ny           - vertical size of the mesh (over y)
  !                    hx           - step of the mesh over x
  !                    hy           - step of the mesh over y
  !                    f(i,j)       - array of values to be interpolated
  !                                   (i = 1, ..., nx; j = 1, ..., ny).
  !
  !                                   For the case of non-rectangular domain:
  !                                   numbers of the mesh points which
  !                                   correspond to the boundaries of the
  !                                   interpolation region:
  !                    imi(j)       - left boundary of the row j (j=1,...,ny)
  !                    ima(j)       - right boundary of the row j (j=1,...,ny)
  !                    jmi(i)       - lower boundary of the column i (i=1,...,nx)
  !                    jma(i)       - upper boundary of the column i (i=1,...,nx)
  !                                   in a rectangle should be:
  !                                   imi(:) = 1
  !                                   ima(:) = nx
  !                                   jmi(:) = 1
  !                                   jma(:) = ny
  !
  !                    icount       - maximum number of entries in spl
  !
  ! Output parameters:
  !                    spl(l,m,k)   - spline coefficients (i,j = 1, ... , n;
  !                    ipoint(i,j)    l,m = 1, ..., 4; i,j - numbers of the
  !                                   mesh point in horizontal and vertical
  !                                   direction (over x and over y), l,m -
  !                                   the numbers of expansion power over x
  !                                   and y (~ dx**(l-1)*dy**(m-1) ))
  !                                   ipoint(i,j) contains the pointer to k
  use spl_three_to_five_sub, only: spl_five_reg
  implicit double precision (a-h,o-z)

  dimension f(nx,ny),spl(6,6,icount),ipoint(nx,ny)
  dimension imi(ny),ima(ny),jmi(nx),jma(nx)

  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi

  nmax=max(nx,ny)

  allocate( ai(nmax),bi(nmax),ci(nmax),di(nmax),ei(nmax),fi(nmax) )

  spl=0.d0
  ipoint=-1

  !  spline along Y-axis

  ic = 0
  do i=1,nx
    if(jmi(i).gt.0) then
      nsi=jma(i)-jmi(i)+1
      do j=jmi(i),jma(i)
        ai(j-jmi(i)+1)=f(i,j)
      end do
      call spl_five_reg(nsi,hy,ai,bi,ci,di,ei,fi)
      do j=jmi(i),jma(i)
        jj=j-jmi(i)+1
        ic = ic+1
        ipoint(i,j)=ic
        spl(1,1,ic)=ai(jj)
        spl(1,2,ic)=bi(jj)
        spl(1,3,ic)=ci(jj)
        spl(1,4,ic)=di(jj)
        spl(1,5,ic)=ei(jj)
        spl(1,6,ic)=fi(jj)
      end do
    end if
  end do

  if (ic .ne. icount) then
    write (6,*) 'Warning, ic, icount:  ',ic,icount
  endif

  !  spline along X-axis

  do j=1,ny
    if(imi(j).gt.0) then
      nsi=ima(j)-imi(j)+1
      do l=1,6
        do i=imi(j),ima(j)
          ai(i-imi(j)+1)=spl(1,l,ipoint(i,j))
        end do
        call spl_five_reg(nsi,hx,ai,bi,ci,di,ei,fi)
        do i=imi(j),ima(j)
          ii=i-imi(j)+1
          spl(2,l,ipoint(i,j))=bi(ii)
          spl(3,l,ipoint(i,j))=ci(ii)
          spl(4,l,ipoint(i,j))=di(ii)
          spl(5,l,ipoint(i,j))=ei(ii)
          spl(6,l,ipoint(i,j))=fi(ii)
        end do
      end do
    end if
  end do

  deallocate( ai,bi,ci,di,ei,fi )

end subroutine s2dcut


! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine spline(nx,ny,x,y,hx,hy,icount,spl,ipoint,xb,yb,u,ux,uy, &
    uxx,uxy,uyy,ierr)
  ! Evaluates interpolated value u(x,y) and its derivatives ux,uy,uxx,uyy,uxy
  ! using arrays calculated by  s2dcut
  ! see also comments in subroutine s2dcut
  !  ierr = 1 - point out of the domain

  implicit double precision (a-h,o-z)

  dimension spl(6,6,icount),x(nx),y(ny),ipoint(nx,ny)
  dimension a(6),ax(6),axx(6)

  xk=(xb-x(1))/hx
  kx=int(xk)+1
  kx = min(nx,max(1,kx))
  yk=(yb-y(1))/hy
  ky=int(yk)+1
  ky = min(ny,max(1,ky))

  ierr=0
  dx=xb-x(kx)
  dy=yb-y(ky)
  do l=1,6
    a(l) =     spl(1,l,ipoint(kx,ky)) &
         + dx*(spl(2,l,ipoint(kx,ky)) &
         + dx*(spl(3,l,ipoint(kx,ky)) &
         + dx*(spl(4,l,ipoint(kx,ky)) &
         + dx*(spl(5,l,ipoint(kx,ky)) &
         + dx* spl(6,l,ipoint(kx,ky))))))
    ax(l) =          spl(2,l,ipoint(kx,ky)) &
          + dx*(2.d0*spl(3,l,ipoint(kx,ky)) &
          + dx*(3.d0*spl(4,l,ipoint(kx,ky)) &
          + dx*(4.d0*spl(5,l,ipoint(kx,ky)) &
          + dx* 5.d0*spl(6,l,ipoint(kx,ky)))))
    axx(l) = 2.d0*spl(3,l,ipoint(kx,ky)) &
           + dx*(6.d0*spl(4,l,ipoint(kx,ky)) &
           + dx*(12.d0*spl(5,l,ipoint(kx,ky)) &
           + dx*(20.d0*spl(6,l,ipoint(kx,ky)))))
  end do

  u = a(1) + dy*(a(2) + dy*(a(3) + dy*(a(4) + dy*(a(5) + dy*a(6)))))
  ux = ax(1) + dy*(ax(2) + dy*(ax(3) + dy*(ax(4) + dy*(ax(5) + dy*ax(6)))))
  uy = a(2) + dy*(2.d0*a(3) + dy*(3.d0*a(4) + dy*(4.d0*a(5) + dy*5.d0*a(6))))
  uxx =  axx(1) + dy*(axx(2) + dy*(axx(3) + dy*(axx(4) + dy*(axx(5) + dy*axx(6)))))
  uxy= ax(2) + dy*(2.d0*ax(3) + dy*(3.d0*ax(4) + dy*(4.d0*ax(5) + dy*5.d0*ax(6))))
  uyy = 2.d0*a(3) + dy*(6.d0*a(4) + dy*(12.d0*a(5) + dy*20.d0*a(6)))

end subroutine spline
