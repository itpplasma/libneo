module ode_integration
  implicit none

  integer, parameter :: dp = kind(1.d0)

  private

  integer, parameter :: MAXSTP=1000000, kmax=0, kmaxx=200
  real(dp), parameter :: TINY_LOCAL=1.e-30
  real(dp), parameter :: SAFETY=0.9d0, PGROW=-.2d0, &
      & PSHRNK=-.25d0, ERRCON=1.89d-4, dxsav=0.d0
  real(dp), dimension(:),   allocatable :: dydx,xp,y,yscal
  real(dp), dimension(:,:), allocatable :: yp
  real(dp), dimension(:),   allocatable :: ak2,ak3,ak4,ak5
  real(dp), dimension(:),   allocatable :: ak6,ytemp
  real(dp), dimension(:),   allocatable :: yerr,ytemp1

  integer, save :: kount=0, ialloc

  abstract interface
    subroutine compute_derivative(x, y, dydx, ierr)
      use libneo_kinds, only : dp
      integer, intent(out) :: ierr
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydx(:)
    end subroutine compute_derivative
  end interface

  !$omp threadprivate(kount, ialloc, dydx, xp, y)
  !$omp threadprivate(yscal, yp, ak2, ak3, ak4, ak5, ak6, ytemp, yerr)
  !$omp threadprivate(ytemp1)

  public :: odeint_allroutines

contains

  subroutine odeint_allroutines(y, nvar, x1, x2, eps, derivs, initial_stepsize)
    use libneo_kinds, only : dp

    implicit none

    procedure(compute_derivative) :: derivs
    integer, intent(in) :: nvar
    real(dp), intent(in) :: x1, x2, eps
    real(dp), intent(in), optional :: initial_stepsize
    real(dp), intent(inout) :: y(nvar)

    integer :: nok,nbad
    real(dp) :: h1,hmin

    if (present(initial_stepsize)) then
      h1 = sign(initial_stepsize,x2-x1)
    else
        h1 = x2 - x1
    end if
    hmin=0.d0

    call odeint(y,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs)

  end subroutine odeint_allroutines

  subroutine alloc_odeint(nvar)
    implicit none

    integer, intent(in) :: nvar

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

  subroutine odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, &
    & derivs)

    use libneo_kinds, only : dp

    implicit none

    procedure(compute_derivative) :: derivs

    integer, intent(in) :: nvar
    integer, intent(out) :: nok, nbad
    real(dp), intent(in) :: x1, x2, eps, h1, hmin
    real(dp), intent(inout) :: ystart(nvar)

    integer :: i,nstp
    integer :: ierr

    real(dp) :: h,hdid,hnext,x,xsav

    xsav = 1.234e5

    ialloc = 1
    call alloc_odeint(nvar)

    x = x1
    h = sign(h1,x2-x1)
    nok = 0
    nbad = 0
    kount = 0

    do i=1,nvar
      y(i) = ystart(i)
    end do

    if (kmax .gt. 0) xsav = x-2.*dxsav

    do nstp=1,MAXSTP
      call derivs(x, y, dydx, ierr)
      if (ierr.ne.0) then
        write(*,*) ierr,'=ierrfield, derivs, odeint'
        ialloc = 0
        call alloc_odeint(nvar)
        return
      end if

      do i=1,nvar
        yscal(i) = abs(y(i)) + abs(h*dydx(i)) + TINY_LOCAL
      end do
      if (kmax .gt. 0) then
        if (abs(x-xsav) .gt. abs(dxsav)) then
          if (kount .lt. kmax-1) then
            kount = kount+1
            xp(kount) = x
            do i=1,nvar
              yp(i,kount) = y(i)
            end do
            xsav = x
          end if
        end if
      end if
      if ((x+h-x2)*(x+h-x1) .gt. 0.0) h=x2-x
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)

      if (abs(hdid-h) <= 0.0d0) then
        nok = nok+1
      else
        nbad = nbad+1
      end if
      if ((x-x2)*(x2-x1) .ge. 0.0) then
        do i=1,nvar
          ystart(i) = y(i)
        end do
        if (kmax .ne. 0) then
          kount = kount+1
          xp(kount) = x
          do i=1,nvar
            yp(i,kount) = y(i)
          end do
        end if
        ialloc = 0
        call alloc_odeint(nvar)
        return
      end if

      if (abs(hnext) .lt. hmin) write(*,*) 'stepsize smaller than minimum in odeint'
      h=hnext
    end do

    write(*,*) 'too many steps in odeint'
    ialloc=0
    call alloc_odeint(nvar)
  end subroutine odeint

  subroutine rkck(y, dydx, n, x, h, yout, yerr, derivs)

    use libneo_kinds, only : dp

    implicit none

    procedure(compute_derivative) :: derivs

    integer, intent(in) :: n
    real(dp), intent(in) :: y(n), dydx(n), x, h
    real(dp), intent(out) :: yout(n), yerr(n)

    integer :: i, ierr
    real(dp), parameter :: A2 = 0.2d0, A3 = 0.3d0, A4 = 0.6d0, &
      & A5 = 1.0d0, A6 = 0.875d0, &
      & B21 = 0.2d0, B31 = 3.0d0/40.0d0, B32 = 9.0d0/40.0d0, B41 = 0.3d0, B42 = -0.9d0, &
      & B43 = 1.2d0, B51 = -11.0d0/54.0d0, B52 = 2.5d0, B53 = -70.0d0/27.0d0, &
      & B54 = 35.0d0/27.0d0, B61 = 1631.0d0/55296.0d0, B62 = 175.0d0/512.0d0, &
      & B63 = 575.0d0/13824.0d0, B64 = 44275.0d0/110592.0d0, B65 = 253.0d0/4096.0d0, &
      & C1 = 37.0d0/378.0d0, C3 = 250.0d0/621.0d0, C4 = 125.0d0/594.0d0, &
      & C6 = 512.0d0/1771.0d0, &
      & DC1 = C1-2825.0d0/27648.0d0, DC3 = C3-18575.0d0/48384.0d0, &
      & DC4 = C4-13525.0d0/55296.0d0, DC5 = -277.0d0/14336.0d0, DC6 = C6-0.25d0

    do i=1,n
      ytemp(i) = y(i) + B21*h*dydx(i)
    end do

    call derivs(x+A2*h, ytemp, ak2, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      ytemp(i) = y(i) + h*(B31*dydx(i) + B32*ak2(i))
    end do

    call derivs(x+A3*h, ytemp, ak3, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      ytemp(i) = y(i) + h*(B41*dydx(i) + B42*ak2(i) + B43*ak3(i))
    end do
    call derivs(x+A4*h, ytemp, ak4, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      ytemp(i) = y(i) + h*(B51*dydx(i) + B52*ak2(i) + B53*ak3(i) + B54*ak4(i))
    end do
    call derivs(x+A5*h, ytemp, ak5, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      ytemp(i) = y(i) + h*(B61*dydx(i) + B62*ak2(i) + B63*ak3(i) &
        & + B64*ak4(i) + B65*ak5(i))
    end do
    call derivs(x+A6*h, ytemp, ak6, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      yout(i) = y(i) + h*(C1*dydx(i) + C3*ak3(i) + C4*ak4(i) + C6*ak6(i))
    end do
    do i=1,n
      yerr(i) = h*(DC1*dydx(i) + DC3*ak3(i) + DC4*ak4(i) + DC5*ak5(i) &
        & + DC6*ak6(i))
    end do

  end subroutine rkck

  subroutine rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
    use libneo_kinds, only : dp

    implicit none

    procedure(compute_derivative) :: derivs

    integer, intent(in) :: n
    real(dp), intent(in) :: dydx(n), htry, eps, yscal(n)
    real(dp), intent(inout) :: y(n), x
    real(dp), intent(out) :: hdid, hnext

    integer :: i

    real(dp) :: errmax, h, htemp, xnew

    h = htry

    do

      call rkck(y, dydx, n, x, h, ytemp1, yerr, derivs)

      errmax=0.0
      do i=1,n
        errmax = max(errmax, abs(yerr(i)/yscal(i)))
      end do
      errmax=errmax/eps

      if (errmax .le. 1.0) exit

      htemp = SAFETY*h*(errmax**PSHRNK)
      h = sign(max(abs(htemp),0.1*abs(h)),h)
      xnew = x+h

      if (abs(xnew - x) <= 0.0d0) then
        print *,'stepsize underflow in rkqs, x,y = ', x, y
        stop
      end if

    end do

    if (errmax .gt. ERRCON) then
      hnext = SAFETY*h*(errmax**PGROW)
    else
      hnext = 5.0*h
    end if
    hdid = h
    x = x+h
    do i=1,n
      y(i)=ytemp1(i)
    end do

  end subroutine rkqs

  subroutine rk4b(y, n, x, h, derivs)

    use libneo_kinds, only : dp

    implicit none

    procedure(compute_derivative) :: derivs

    integer, intent(in) :: n
    real(dp), intent(in) :: x, h
    real(dp), intent(inout) :: y(n)

    integer, parameter :: NMAX=12
    integer :: i, ierr
    real(dp) :: dydx(NMAX), yt(NMAX), dyt(NMAX), dym(NMAX)
    real(dp) :: xh, hh, h6

    hh = h*0.5d0

    h6 = h/6.0d0
    xh = x + hh

    call derivs(x, y, dydx, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      yt(i) = y(i) + hh*dydx(i)
    end do

    call derivs(xh, yt, dyt, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      yt(i) = y(i) + hh*dyt(i)
    end do

    call derivs(xh, yt, dym, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      yt(i) = y(i) + h*dym(i)
      dym(i)= dyt(i) + dym(i)
    end do

    call derivs(x+h, yt, dyt, ierr)
    if (ierr .ne. 0) return

    do i=1,n
      y(i) = y(i) + h6*(dydx(i) + dyt(i) + 2.0d0*dym(i))
    end do

  end subroutine rk4b

end module ode_integration
