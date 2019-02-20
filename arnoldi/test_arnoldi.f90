program test_arnoldi
  use arnoldi_mod, only: ieigen,ngrow,tol,eigvecs,arnoldi_ierr=>ierr
  use for_mpi, only: mype, npes, mpi_p_root

  integer :: ierr, info
  integer, parameter :: nsize = 10
  integer ipiv(nsize)
  double complex :: amat(nsize,nsize), mmat(nsize,nsize)
  double complex :: bvec(nsize), yvec(nsize), xsol(nsize), xold(nsize), xnew(nsize)
  double complex :: cdummy(nsize)
  integer :: kit, maxit = 20

  ! for Arnoldi
  integer, parameter :: nritz = 8  ! for Arnoldi iterations
  double complex, dimension(nritz) :: ritznum
  double complex, allocatable, dimension(:) :: coefren
  double complex, allocatable, dimension(:,:) :: amat2, bvec2
  
#ifdef PARALLEL
  call MPI_INIT(ierr)
  if(ierr .ne. MPI_SUCCESS) then 
     write(*,*)'MPI init. err.'
     stop 1
  endif
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
#else
  npes = 1
  mype = 0
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
        !
  if(info.ne.0) then
     if(info.gt.0) then
        print *,'iterator: singular matrix in zgesv'
     else
        print *,'iterator: argument ',-info,' has illegal value in zgesv'
     endif
     return
  endif

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
  ieigen=1
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
     enddo
  enddo
  !
  call zgesv(ngrow,ngrow,amat2,ngrow,ipiv,bvec2,ngrow,info)
  !
  if(info.ne.0) then
     if(info.gt.0) then
        print *,'iterator: singular matrix in zgesv'
     else
        print *,'iterator: argument ',-info,' has illigal value in zgesv'
     endif
     !deallocate(coefren,amat,bvec,ipiv)
     !deallocate(eigvecs)
     return
  endif
  
  xold = (0d0,0d0)
  do kit = 1,maxit
     call next_iteration(nsize, xold, xnew)
     
     do j=1,ngrow
        coefren(j)=ritznum(j)*sum(bvec2(j,:)                           &
             *matmul(transpose(conjg(eigvecs(:,1:ngrow))),xnew-xold))
     enddo
     xnew = xnew - matmul(eigvecs(:,1:ngrow),coefren(1:ngrow))
     xold = xnew
     print *, NORM2([NORM2(real(xnew-xsol)),NORM2(aimag(xnew-xsol))]) 
  end do
contains

  
  subroutine next_iteration(n, hold, hnew)
    integer :: n
    double complex :: hold(n), hnew(n), x(n), y(n)
    double complex :: alpha, beta
    x = hold+yvec
    y = cmplx(0d0,0d0)
    alpha = cmplx(1d0,0d0)
    beta = cmplx(0d0,0d0)
    call zgemv('N',n,n,alpha,mmat,n,x,1,beta,y,1)
    hnew = y
  end subroutine next_iteration
  
  
end program test_arnoldi
