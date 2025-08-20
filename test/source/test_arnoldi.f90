module test_arnoldi_mod
#ifdef PARALLEL
  use mpiprovider_module, only : mpro
#endif
  use arnoldi, only: calc_ritz_eigenvalues, leigen,ngrow,tol,eigvecs
  use libneo_kinds, only : cdp
  implicit none

  integer :: info
  integer :: ios
  integer, parameter :: nsize = 10
  integer :: ipiv(nsize)
  complex(cdp) :: amat(nsize,nsize), mmat(nsize,nsize)
  complex(cdp) :: bvec(nsize), yvec(nsize), xsol(nsize), xold(nsize), xnew(nsize)
  integer :: kit, maxit = 20

  ! for Arnoldi
  integer, parameter :: nritz = 8  ! for Arnoldi iterations
  complex(cdp), dimension(nritz) :: ritznum
  complex(cdp), allocatable, dimension(:) :: coefren
  complex(cdp), allocatable, dimension(:,:) :: amat2, bvec2

contains

  subroutine next_iteration(n, hold, hnew)
    implicit none

    integer, intent(in) :: n
    complex(cdp), intent(in) :: hold(n)
    complex(cdp), intent(out) :: hnew(n)

    complex(cdp) :: alpha, beta, x(n), y(n)

    x = hold + yvec(1:n)
    y = cmplx(0d0,0d0,cdp)
    alpha = cmplx(1d0,0d0,cdp)
    beta = cmplx(0d0,0d0,cdp)
    call zgemv('N',n,n,alpha,mmat,n,x,1,beta,y,1)
    hnew = y
  end subroutine next_iteration

  subroutine run_test()
    integer :: i, j

#ifdef PARALLEL
    call mpro%init()
#endif

    open(1,file='amat.dat', action='read', status='old', iostat=ios)
    if (ios .ne. 0) stop 'Error while trying to open amat.dat'
    read(1,*) amat
    close(1)

    open(1,file='mmat.dat', action='read', status='old', iostat=ios)
    if (ios .ne. 0) stop 'Error while trying to open mmat.dat'
    read(1,*) mmat
    close(1)

    open(1,file='yvec.dat', action='read', status='old', iostat=ios)
    if (ios .ne. 0) stop 'Error while trying to open yvec.dat'
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
    call calc_ritz_eigenvalues(nsize, nritz, ritznum, next_iteration)
    do kit = 1, ngrow
      print *, ritznum(kit)
    end do

    print *, "Solving subsystem in eigenbasis"
    allocate(amat2(ngrow,ngrow),bvec2(ngrow,ngrow),coefren(ngrow))
    bvec2=(0.d0,0.d0)
    do i=1,ngrow
      do j=1,ngrow
        bvec2(i,j)=sum(bvec*conjg(eigvecs(:,i)))
      end do
    end do

    do i=1,ngrow
      do j=1,ngrow
        amat2(i,j)=ritznum(i)*bvec2(i,j)
      end do
    end do

    xold = (0d0,0d0)
    do kit = 1,maxit
      call next_iteration(nsize, xold, xnew)

      do j=1,ngrow
        coefren(j)=ritznum(j)*sum(bvec2(j,:) &
             *matmul(transpose(conjg(eigvecs(:,1:ngrow))),xnew-xold))
      end do
      xnew = xnew - matmul(eigvecs(:,1:ngrow),coefren(1:ngrow))
      xold = xnew
      print *, NORM2([NORM2(real(xnew-xsol)),NORM2(aimag(xnew-xsol))])
    end do

  end subroutine run_test

end module test_arnoldi_mod

program test_arnoldi
  use test_arnoldi_mod
  call run_test()
end program test_arnoldi