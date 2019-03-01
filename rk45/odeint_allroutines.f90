module odeint_mod
  integer :: kmax=0, kount=0, kmaxx=200, ialloc
  double precision :: dxsav=0.d0
  double precision, dimension(:),   allocatable :: dydx,xp,y,yscal
  double precision, dimension(:,:), allocatable :: yp
  double precision, dimension(:),   allocatable :: ak2,ak3,ak4,ak5
  double precision, dimension(:),   allocatable :: ak6,ytemp
  double precision, dimension(:),   allocatable :: yerr,ytemp1

!$omp threadprivate(kmax, kount, kmaxx, ialloc, dxsav, dydx, xp, y)
!$omp threadprivate(yscal, yp, ak2, ak3, ak4, ak5, ak6, ytemp, yerr)
!$omp threadprivate(ytemp1)

contains

  subroutine odeint_allroutines(y,nvar,x1,x2,eps,derivs)

    implicit none

    external :: derivs,rkqs
    integer :: nvar,nok,nbad
    double precision :: x1,x2,eps,h1,hmin
    double precision, dimension(nvar) :: y

    h1=x2-x1
    hmin=0.d0

    call odeint(y,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

    return
  end subroutine odeint_allroutines

  subroutine alloc_odeint(nvar)

    if (ialloc.eq.1) then
      allocate(dydx(nvar),xp(kmaxx),y(nvar))
      allocate(yp(nvar,kmaxx),yscal(nvar))
      allocate(ak2(nvar),ak3(nvar),ak4(nvar),ak5(nvar))
      allocate(ak6(nvar),ytemp(nvar))
      allocate(yerr(nvar),ytemp1(nvar))
    else
      deallocate(dydx,xp,y,yp,yscal)
      deallocate(ak2,ak3,ak4,ak5,ak6,ytemp)
      deallocate(yerr,ytemp1)
    end if

  end subroutine alloc_odeint

  subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

    use gbpi_mod, only : ierrfield

    implicit double precision (a-h,o-z)

    integer nbad,nok,nvar,KMAXX,MAXSTP
    double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
    external derivs,rkqs
    parameter (MAXSTP=1000000,TINY=1.e-30)
    integer i,nstp
    double precision h,hdid,hnext,x,xsav

    ialloc=1
    call alloc_odeint(nvar)

    x=x1
    h=sign(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    do i=1,nvar
      y(i)=ystart(i)
    end do

    if (kmax.gt.0) xsav=x-2.*dxsav

    do nstp=1,MAXSTP
      call derivs(x,y,dydx)

      if (ierrfield.ne.0) then
        print*,ierrfield,'=ierrfield, derivs, odeint'
        ialloc=0
        call alloc_odeint(nvar)
        return
      end if

      do i=1,nvar
        yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
      end do
      if (kmax.gt.0) then
        if (abs(x-xsav).gt.abs(dxsav)) then
          if (kount.lt.kmax-1) then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
            xsav=x
          end if
        end if
      end if
      if ((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)

      if (ierrfield.ne.0) then
        print*,ierrfield,'=ierrfield, rkqs, odeint'
        ialloc=0
        call alloc_odeint(nvar)
        return
      end if

      if (hdid.eq.h) then
        nok=nok+1
      else
        nbad=nbad+1
      end if
      if ((x-x2)*(x2-x1).ge.0.) then
        do i=1,nvar
          ystart(i)=y(i)
        end do
        if (kmax.ne.0) then
          kount=kount+1
          xp(kount)=x
          do i=1,nvar
            yp(i,kount)=y(i)
          end do
        end if
        ialloc=0
        call alloc_odeint(nvar)
        return
      end if

      if (abs(hnext).lt.hmin) then
        print*,'stepsize smaller than minimum in odeint'
        write(99,*)'stepsize smaller than minimum in odeint, stop'
        stop
      end if

      h=hnext
    end do

    print*,'too many steps in odeint'
    write(99,*)'to many steps in odeint, stop'
    stop


    ialloc=0
    call alloc_odeint(nvar)
    return
  end subroutine odeint

  subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)

    use gbpi_mod

    implicit double precision (a-h,o-z)

    integer n
    double precision h,x,dydx(n),y(n),yerr(n),yout(n)
    external derivs

    integer i
    double precision A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
      & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,
     *B21=.2d0,B31=3.d0/40.d0,
     *B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,
     *B51=-11.d0/54.d0,B52=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,
     *B62=175.d0/512.d0,
     *B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,
     *C1=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,
     *DC1=C1-2825.d0/27648.d0,
     *DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,
     *DC5=-277.d0/14336.d0,
     *DC6=C6-.25d0)
    do i=1,n
      ytemp(i)=y(i)+B21*h*dydx(i)
    end do
    call derivs(x+A2*h,ytemp,ak2)

    if (ierrfield.ne.0) return

    do i=1,n
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
    end do
    call derivs(x+A3*h,ytemp,ak3)

    if (ierrfield.ne.0) return

    do i=1,n
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    end do
    call derivs(x+A4*h,ytemp,ak4)

    if (ierrfield.ne.0) return

    do i=1,n
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
    end do
    call derivs(x+A5*h,ytemp,ak5)

    if (ierrfield.ne.0) return

    do i=1,n
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do
    call derivs(x+A6*h,ytemp,ak6)

    if (ierrfield.ne.0) return

    do i=1,n
      yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
    end do
    do i=1,n
      yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
    end do

    return
  end subroutine

  subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)

    implicit double precision (a-h,o-z)

    integer n
    double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
    external derivs

    integer i
    double precision :: errmax,h,htemp,xnew,SAFETY,PGROW, PSHRNK,ERRCON
    parameter (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,  ERRCON=1.89d-4)
    h=htry
1   call rkck(y,dydx,n,x,h,ytemp1,yerr,derivs)
    errmax=0.
    do i=1,n
      errmax=max(errmax,abs(yerr(i)/yscal(i)))
    end do
    errmax=errmax/eps
    if (errmax.gt.1.) then
      htemp=SAFETY*h*(errmax**PSHRNK)
      h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h

      if (xnew.eq.x) then
        print *,'stepsize underflow in rkqs, x,y = ',x,y

        write(99,*)'stepsize underflow in rkqs, x,y = ',x,y

        stop
      end if
      goto 1
    else
      if(errmax.gt.ERRCON)then
        hnext=SAFETY*h*(errmax**PGROW)
      else
        hnext=5.*h
      endif
      hdid=h
      x=x+h
      do i=1,n
        y(i)=ytemp1(i)
      end do
      return
    end if
  end subroutine rkqs


  subroutine RK4b(Y,N,X,H,DERIVS)

    use gbpi_mod

    implicit double precision (a-h,o-z)
    external derivs
    parameter (NMAX=12)
    dimension Y(N),DYDX(NMAX),YT(NMAX),DYT(NMAX),DYM(NMAX)

    HH=H*0.5d0

    H6=H/6.0d0
    XH=X+HH

    call DERIVS(X,Y,DYDX)
    if(ierrfield.ne.0) return

    do I=1,N
      YT(I)=Y(I)+HH*DYDX(I)
    end do

    call DERIVS(XH,YT,DYT)
    if (ierrfield.ne.0) return

    do I=1,N
      YT(I)=Y(I)+HH*DYT(I)
    end do

    call DERIVS(XH,YT,DYM)
    if (ierrfield.ne.0) return

    do I=1,N
      YT(I)=Y(I)+H*DYM(I)
      DYM(I)=DYT(I)+DYM(I)
    end do

    call DERIVS(X+H,YT,DYT)
    if (ierrfield.ne.0) return

    do I=1,N
      Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.d0*DYM(I))
    end do

    return
  end subroutine RK4b

end module odeint_mod
