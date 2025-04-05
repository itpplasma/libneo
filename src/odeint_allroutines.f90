!
!------------------------------------------------------------------------------
module odeint_mod
    integer :: kmax=0, kount=0, kmaxx=200, ialloc
    double precision :: dxsav=0.d0
    double precision, dimension(:),   allocatable :: dydx,xp,y,yscal
    double precision, dimension(:,:), allocatable :: yp
    double precision, dimension(:),   allocatable :: ak2,ak3,ak4,ak5
    double precision, dimension(:),   allocatable :: ak6,ytemp
    double precision, dimension(:),   allocatable :: yerr,ytemp1

    !$omp threadprivate(kount, ialloc, dydx, xp, y)
    !$omp threadprivate(yscal, yp, ak2, ak3, ak4, ak5, ak6, ytemp, yerr)
    !$omp threadprivate(ytemp1)
end module odeint_mod

module odeint_allroutines_sub
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    abstract interface
    subroutine compute_derivative(x, y, dydx)
        import :: dp
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
    end subroutine compute_derivative
  end interface

    contains
!
      SUBROUTINE odeint_allroutines(y,nvar,x1,x2,eps,derivs,initial_stepsize)
!
      implicit none
!
      procedure(compute_derivative) :: derivs
      integer :: nvar,nok,nbad
      double precision :: x1,x2,eps,h1,hmin
      double precision, optional :: initial_stepsize
      double precision, dimension(nvar) :: y
!
      if (present(initial_stepsize)) then
        h1 = sign(initial_stepsize,x2-x1)
      else
        h1 = x2 - x1
      end if
      hmin=0.d0
!
      call odeint(y,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs)
!
      return
      END SUBROUTINE odeint_allroutines
!
!------------------------------------------------------------------------------
      SUBROUTINE alloc_odeint(nvar)
!
      use odeint_mod
!
      integer, intent(in) :: nvar
!
      if(ialloc.eq.1) then
        allocate(dydx(nvar),xp(kmaxx),y(nvar))
        allocate(yp(nvar,kmaxx),yscal(nvar))
        allocate(ak2(nvar),ak3(nvar),ak4(nvar),ak5(nvar))
        allocate(ak6(nvar),ytemp(nvar))
        allocate(yerr(nvar),ytemp1(nvar))
      else
        deallocate(dydx,xp,y,yp,yscal)
        deallocate(ak2,ak3,ak4,ak5,ak6,ytemp)
        deallocate(yerr,ytemp1)
      endif
!
      END SUBROUTINE alloc_odeint
!------------------------------------------------------------------------------
!
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs)
!
      use odeint_mod, only : kmax,kount,ialloc,dxsav,dydx,xp,y,yscal,yp
!
      implicit double precision (a-h,o-z)
!
      INTEGER nbad,nok,nvar,MAXSTP
      double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
      procedure(compute_derivative) :: derivs
      PARAMETER (MAXSTP=1000000,TINY=1.e-30)
      INTEGER i,nstp
      double precision h,hdid,hnext,x,xsav
!
      xsav = 1.234e5
!
      ialloc=1
      call alloc_odeint(nvar)
!
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(abs(hdid-h).le.0.0d0)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          ialloc=0
          call alloc_odeint(nvar)
          return
        endif
        if(abs(hnext).lt.hmin) then
          error stop 'stepsize smaller than minimum in odeint'
        endif
        h=hnext
16    continue
      error stop 'stepsize smaller than minimum in odeint'
      ialloc=0
      call alloc_odeint(nvar)
      return
      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
!
      use odeint_mod, only : ak2,ak3,ak4,ak5,ak6,ytemp
!
      implicit double precision (a-h,o-z)
!
      INTEGER n
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      procedure(compute_derivative) :: derivs
!U    USES derivs
      INTEGER i
      double precision &
     A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,&
     B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0, &
     B21=.2d0,B31=3.d0/40.d0, &
     B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0, &
     B51=-11.d0/54.d0,B52=2.5d0, &
     B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0, &
     B62=175.d0/512.d0, &
     B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0, &
     C1=37.d0/378.d0, &
     C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0, &
     DC1=C1-2825.d0/27648.d0, &
     DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0, &
     DC5=-277.d0/14336.d0, &
     DC6=C6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+&
     B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*&
     ak6(i))
17    continue
      return
      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
!
      use odeint_mod, only : yerr,ytemp1
!
      implicit double precision (a-h,o-z)
!
      INTEGER n
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      procedure(compute_derivative) :: derivs
!U    USES derivs,rkck
      INTEGER i
      double precision errmax,h,htemp,xnew,SAFETY,PGROW, PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp1,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
!        if(xnew.eq.x) pause 'stepsize underflow in rkqs'
        if(abs(xnew-x).le.0.0d0) then
          print *,'stepsize underflow in rkqs, x,y = ',x,y
          stop
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp1(i)
12      continue
        return
      endif
      END
end module odeint_allroutines_sub
