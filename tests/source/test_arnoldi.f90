program test_arnoldi
  use mpiprovider_module, only : mpro

  use arnoldi_mod, only: arnoldi, leigen,ngrow,tol,eigvecs
  use libneo_kinds, only : complex_kind

  integer :: info
  integer, parameter :: nsize = 10
  integer ipiv(nsize)
  complex(kind=complex_kind) :: amat(nsize,nsize), mmat(nsize,nsize)
  complex(kind=complex_kind) :: bvec(nsize), yvec(nsize), xsol(nsize), xold(nsize), xnew(nsize)
  integer :: kit, maxit = 20

  ! for Arnoldi
  integer, parameter :: nritz = 8  ! for Arnoldi iterations
  complex(kind=complex_kind), dimension(nritz) :: ritznum
  complex(kind=complex_kind), allocatable, dimension(:) :: coefren
  complex(kind=complex_kind), allocatable, dimension(:,:) :: amat2, bvec2

#ifdef PARALLEL
  call mpro%init()
#endif

  open(1,file='amat.dat')
  read(1,*) amat
  close(1)

  open(1,file='bvec.dat')
  read(1,*) bvec
  close(1)

  open(1,file='mmat.dat')
  read(1,*) mmat
  close(1)

  open(1,file='yvec.dat')
  read(1,*) yvec
  close(1)

  call zgesv(nsize,1,amat,nsize,ipiv,bvec,nsize,info)

  if (info.ne.0) then
    if (info.gt.0) then
      print *,'iterator: singular matrix in zgesv'
    else
      print *,'iterator: argument ',-info,' has illegal value in zgesv'
    end if
    stop
  end if

  xsol = bvec ! reference solution

  print *, xsol
  print *
  print *

  ! direct iteration
  xold = (0d0,0d0)
  do kit = 1,maxit
    call next_iteration(nsize, xold, xnew)
    xold = xnew
    print *, NORM2([NORM2(real(xnew-xsol)),NORM2(aimag(xnew-xsol))])
  end do

  ! Arnoldi
  print *, "Finding eigenvalues"
  leigen= .true.
  tol = 0.7d0
  call arnoldi(nsize, nritz, ritznum, next_iteration)
  do kit = 1, ngrow
    print *, ritznum(kit)
  end do

  print *, "Solving subsystem in eigenbasis"
  allocate(amat2(ngrow,ngrow),bvec2(ngrow,ngrow),coefren(ngrow))
  bvec2=(0.d0,0.d0)
  do i=1,ngrow
    bvec2(i,i)=(1.d0,0.d0)
    do j=1,ngrow
      amat2(i,j)=sum(conjg(eigvecs(:,i))*eigvecs(:,j))*(ritznum(j)-(1.d0,0.d0))
    end do
  end do

  call zgesv(ngrow,ngrow,amat2,ngrow,ipiv,bvec2,ngrow,info)

  if (info.ne.0) then
    if (info.gt.0) then
      print *,'iterator: singular matrix in zgesv'
    else
      print *,'iterator: argument ',-info,' has illegal value in zgesv'
    end if

    stop
  end if

  xold = (0d0,0d0)
  do kit = 1,maxit
    call next_iteration(nsize, xold, xnew)

    do j=1,ngrow
      coefren(j)=ritznum(j)*sum(bvec2(j,:)                           &
           *matmul(transpose(conjg(eigvecs(:,1:ngrow))),xnew-xold))
    end do
    xnew = xnew - matmul(eigvecs(:,1:ngrow),coefren(1:ngrow))
    xold = xnew
    print *, NORM2([NORM2(real(xnew-xsol)),NORM2(aimag(xnew-xsol))])
  end do

contains

  subroutine next_iteration(n, hold, hnew)
    use libneo_kinds, only : real_kind, complex_kind

    implicit none

    integer, intent(in) :: n
    complex(kind=complex_kind), intent(in) :: hold(n)
    complex(kind=complex_kind), intent(out) :: hnew(n)

    complex(kind=complex_kind) :: alpha, beta, x(n), y(n)

    x = hold+yvec
    y = cmplx(0d0,0d0)
    alpha = cmplx(1d0,0d0)
    beta = cmplx(0d0,0d0)
    call zgemv('N',n,n,alpha,mmat,n,x,1,beta,y,1)
    hnew = y
  end subroutine next_iteration

end program test_arnoldi
