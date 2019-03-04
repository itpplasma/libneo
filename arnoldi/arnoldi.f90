
module arnoldi_mod
  use libneo_kinds, only : real_kind, complex_kind

  integer :: ieigen=0
  integer :: ngrow,ierr
  real(kind=real_kind) :: tol
  complex(kind=complex_kind), dimension(:,:), allocatable :: eigvecs

contains

  !> Computes m Ritz eigenvalues (approximations to extreme eigenvalues)
  !> of the iteration procedure of the vector with dimension n.
  !> Eigenvalues are computed by means of Arnoldi iterations.
  !> Optionally computes Ritz vectors (approximation to eigenvectors).
  !>
  !> Input  parameters:
  !> Formal:             n              - system dimension
  !>                     m              - number of Ritz eigenvalues
  !> Module arnoldi_mod: ieigen         - flag to compute eigenvectors (1 - yes,
  !>                                      0 - no), module is initialized with 0
  !>                     tol            - eigenvectors are not computed for
  !>                                      eigenvalues smaller than this numer
  !> External:           next_iteration - routine computing next iteration
  !>                                      of the solution from the previous
  !> Output parameters:
  !> Formal:             ritznum        - Ritz eigenvalues
  !> Module arnoldi_mod: ngrow          - number of eigenvalues larger or equal
  !>                                      to TOL
  !>                     eigvecs        - array of eigenvectors, size - (m,ngrow)
  !>                     ierr           - error code (0 - normal work, 1 - error)
  subroutine arnoldi(n,m,ritznum,next_iteration)

    use for_mpi, only : mype, mpi_p_root
    use libneo_kinds, only : real_kind, complex_kind
#ifdef PARALLEL
    use mpi
#endif

    implicit none

    external :: next_iteration
    integer                                       :: n,m,k,j,lwork,info
    complex(kind=complex_kind)                    :: tmp
    complex(kind=complex_kind), dimension(m)      :: ritznum
    logical,          dimension(:),   allocatable :: select
    integer,          dimension(:),   allocatable :: ifailr
    real(kind=real_kind), dimension(:),   allocatable :: rwork

    complex(kind=complex_kind), dimension(:),   allocatable :: fold,fnew,fzero,work,rnum
    complex(kind=complex_kind), dimension(:,:), allocatable :: qvecs,hmat,eigh
    complex(kind=complex_kind), dimension(:,:), allocatable :: qvecs_save

    allocate(fold(n),fnew(n),fzero(n))
    fold=(0.d0,0.d0)
    call next_iteration(n,fold,fnew)
    fzero=fnew

    if (mype == mpi_p_root) then
      ierr=0
      allocate(qvecs(n,m),hmat(m,m))
      hmat=(0.d0,0.d0)
      qvecs(:,1)=fnew/sqrt(sum(conjg(fnew)*fnew))
    end if

    do k=2,m
      if (mype == mpi_p_root) fold=qvecs(:,k-1)
#ifdef PARALLEL
      call MPI_BCAST(fold, n, MPI_DOUBLE_COMPLEX, mpi_p_root, MPI_COMM_WORLD, ierr)
#endif
      !print *, mype, sum(fold)
      call next_iteration(n,fold,fnew)
      !print *, mype, sum(fnew)
      if (mype == mpi_p_root) then
        qvecs(:,k)=fnew-fzero
        do j=1,k-1
          hmat(j,k-1)=sum(conjg(qvecs(:,j))*qvecs(:,k))
          qvecs(:,k)=qvecs(:,k)-hmat(j,k-1)*qvecs(:,j)
        end do
        hmat(k,k-1)=sqrt(sum(conjg(qvecs(:,k))*qvecs(:,k)))
        qvecs(:,k)=qvecs(:,k)/hmat(k,k-1)
      end if
    end do

    if (mype == mpi_p_root) then
      tol=0.7d0
      allocate(eigh(m,m))
      print *,'in'
      print *,m,size(hmat,1),size(hmat,2),size(ritznum),size(eigh,1),size(eigh,2)
      call try_eigvecvals(m,tol,hmat,ngrow,ritznum,eigh,ierr)
      print *,'out',m,ngrow

      if(allocated(eigvecs)) deallocate(eigvecs)
      allocate(eigvecs(n,ngrow))

      eigvecs=matmul(qvecs,eigh(:,1:ngrow))

      deallocate(qvecs,hmat,eigh)
    end if
    deallocate(fold,fnew,fzero)

  end subroutine arnoldi

  !> Computes eigenvalues, ritznum, of the upper Hessenberg matrix hmat
  !> of the dimension (m,m), orders eigenvelues into the decreasing by module
  !> sequence and computes the eigenvectors, eigh, for eigenvalues exceeding
  !> the tolerance tol (number of these eigenvalues is ngrow)
  !>
  !> Input arguments:
  !>          Formal: m        - matrix size
  !>                  tol      - tolerance
  !>                  hmat     - upper Hessenberg matrix
  !> Output arguments:
  !>          Formal: ngrow    - number of exceeding the tolerance
  !>                  ritznum  - eigenvalues
  !>                  eigh     - eigenvectors
  !>                  ierr     - error code (0 - normal work)
  subroutine try_eigvecvals(m,tol,hmat,ngrow,ritznum,eigh,ierr)
    use libneo_kinds, only : real_kind, complex_kind

    implicit none

    integer :: m,ngrow,ierr,k,j,lwork,info

    real(kind=real_kind) :: tol
    complex(kind=complex_kind)   :: tmp

    complex(kind=complex_kind), dimension(m)   :: ritznum
    complex(kind=complex_kind), dimension(m,m) :: hmat,eigh

    logical,          dimension(:),   allocatable :: selec
    integer,          dimension(:),   allocatable :: ifailr
    real(kind=real_kind), dimension(:),   allocatable :: rwork
    complex(kind=complex_kind),   dimension(:),   allocatable :: work,rnum
    complex(kind=complex_kind),   dimension(:,:), allocatable :: hmat_work

    print *,size(hmat)
    print *,size(ritznum)
    print *,size(eigh)
    ierr=0

    allocate(hmat_work(m,m))

    hmat_work=hmat

    allocate(work(1))
    lwork=-1

    call zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    if (info.ne.0) then
      if (info.gt.0) then
        print *,'arnoldi: zhseqr failed to compute all eigenvalues'
      else
        print *,'arnoldi: argument ',-info,' has illigal value in zhseqr'
      endif
      deallocate(hmat_work,work)
      ierr=1
      return
    endif

    lwork=work(1)
    deallocate(work)
    allocate(work(lwork))
    print *,'lwork = ',lwork

    call zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    if(info.ne.0) then
      if(info.gt.0) then
        print *,'arnoldi: zhseqr failed to compute all eigenvalues'
      else
        print *,'arnoldi: argument ',-info,' has illigal value in zhseqr'
      endif
      deallocate(hmat_work,work)
      ierr=1
      return
    endif

    do k=1,m
      info=0
      do j=2,m
        if (abs(ritznum(j)).gt.abs(ritznum(j-1))) then
          info=1
          tmp=ritznum(j-1)
          ritznum(j-1)=ritznum(j)
          ritznum(j)=tmp
        end if
      end do
      if(info.eq.0) exit
    end do


    ! compute how many eigenvalues exceed the tolerance (TOL):

    allocate(selec(m),rnum(m))
    selec=.false.
    ngrow=0
    do j=1,m
      if(abs(ritznum(j)).lt.tol) exit
      ngrow=ngrow+1
      selec(j)=.true.
    end do
    rnum=ritznum
    hmat_work=hmat
    deallocate(work)
    allocate(work(m*m),rwork(m),ifailr(m))
    eigh=(0.d0,0.d0)

    call zhsein('R','Q','N',selec,m,hmat_work,m,rnum,rnum,1,eigh(:,1:ngrow),m,  &
                ngrow,ngrow,work,rwork,ifailr,ifailr,info)

    if (info.ne.0) then
      if (info.gt.0) then
        print *,'arnoldi: ',info,' eigenvectors not converged in zhsein'
      else
        print *,'arnoldi: argument ',-info,' has illigal value in zhsein'
      end if
      ierr=1
    end if

    print *,'ierr = ',ierr

    deallocate(hmat_work,work,rwork,selec,rnum,ifailr)

  end subroutine try_eigvecvals

end module arnoldi_mod
