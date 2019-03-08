!> \brief Module for arnoldi algorithm for computation of Ritz eigenvalues.
!>
!> This module implements an arnoldi algorithm for the computation of
!> Ritz eigenvlues.
!> Parameters can be set via a namelist 'arnoldi', see subroutine
!> arnoldi_namelist. This routine is also responsible for initializing
!> the module.
!> Main work is done by the subroutine arnoldi.
module arnoldi_mod
  use libneo_kinds, only : real_kind, complex_kind

  logical, save :: leigen=.false. !< Determines if eigenvectors should be calculated.
  integer, save :: ngrow,ierr
  real(kind=real_kind), save :: tol !< Tolerance below which eigenvectors are no longer calculated.
  complex(kind=complex_kind), dimension(:,:), allocatable :: eigvecs

  abstract interface
    subroutine interface_next_iteration(n, old, new)
      use libneo_kinds, only : complex_kind
      integer, intent(in) :: n
      complex(kind=complex_kind), dimension(n), intent(in) :: old
      complex(kind=complex_kind), dimension(n), intent(out) :: new
    end subroutine interface_next_iteration
  end interface
contains

  !> \brief Init module, read/write namelist.
  !> \author Rico Buchholz
  !>
  !> This subroutine serves three purposes (unfortunately, so far I have
  !> no idea how to solve this better).
  !> This routine initializes the module, and is responsible for reading
  !> and writing of the settings from/to a namelist file.
  !>
  !> Reading from the namelist will be done if a file unit is given and
  !> write_to_file_ is absent or false.
  !> Writing to the namelist will be done if a file unit is given and
  !> write_to_file_ is true.
  !> Initializing will be done whenever write_to_file_ is absent or
  !> false, i.e. you can not read from the namelist without the
  !> initialization prozess.
  !>
  !> \note Reading/writing of the namelist are considered optional, so
  !>   if reading/writing is not possible, the code will still continue.
  !>
  !> \note This subroutine does not call the arnoldi subroutine. Doxygen
  !>   just thinks so, because the namelist has the same name.
  !>
  !> \param[in] namelist_file_unit_number
  !>   Optional file unit, that corresponds to an open file.
  !> \param[in] write_to_file_
  !>   Optional logical, if it is present, true and also
  !>   namelist_file_unit_number is present, then the namelist will be
  !>   written to the file, instead of read.
  subroutine arnoldi_namelist(namelist_file_unit_number, write_to_file_)
    integer, intent(in), optional :: namelist_file_unit_number
    logical, intent(in), optional :: write_to_file_

    integer :: ios
    logical :: write_to_file

    namelist /arnoldi/  &
      & leigen, tol

    if (present(write_to_file_)) then
      write_to_file = write_to_file_
    end if

    if (write_to_file) then
      leigen = .false.
      ngrow = 0
      ierr = 0
      tol = 0.7
      write_to_file = .false.
    end if

    if (present(namelist_file_unit_number)) then
      if (write_to_file) then
        write(namelist_file_unit_number, nml=arnoldi, iostat=ios)
      else
        read(namelist_file_unit_number, nml=arnoldi, iostat=ios)
      end if
      if (ios .ne. 0) then
        write(*,*) "WARNING: problem reading namelist 'arnoldi'."
        write(*,*) "         Namelist is optional, i will continue."
      end if
    end if

  end subroutine arnoldi_namelist

  !> Computes m Ritz eigenvalues (approximations to extreme eigenvalues)
  !> of the iteration procedure of the vector with dimension n.
  !> Eigenvalues are computed by means of Arnoldi iterations.
  !> Optionally computes Ritz vectors (approximation to eigenvectors).
  !>
  !> Input  parameters:
  !> Formal:             n              - system dimension
  !>                     m              - number of Ritz eigenvalues
  !> (Subroutine)        next_iteration - routine computing next iteration
  !>                                      of the solution from the previous
  !> Module arnoldi_mod: leigen         - logical to compute eigenvectors
  !>                     tol            - eigenvectors are not computed for
  !>                                      eigenvalues smaller than this number
  !> Output parameters:
  !> Formal:             ritznum        - Ritz eigenvalues
  !> Module arnoldi_mod: ngrow          - number of eigenvalues larger or equal
  !>                                      to TOL
  !>                     eigvecs        - array of eigenvectors, size - (m,ngrow)
  !>                     ierr           - error code (0 - normal work, 1 - error)
  subroutine arnoldi(n,m,ritznum,next_iteration)

    use mpiprovider_module, only : mpro
    use libneo_kinds, only : real_kind, complex_kind
#ifdef PARALLEL
    use mpi
#endif

    implicit none

    procedure(interface_next_iteration) :: next_iteration
    integer, intent(in) :: n,m
    complex(kind=complex_kind), dimension(m), intent(out) :: ritznum

    integer :: k,j
    integer, parameter :: mpi_p_root = 0

    complex(kind=complex_kind), dimension(:),   allocatable :: fold,fnew,fzero
    complex(kind=complex_kind), dimension(:,:), allocatable :: qvecs,hmat,eigh

    allocate(fold(n),fnew(n),fzero(n))
    fold=(0.d0,0.d0)
    call next_iteration(n,fold,fnew)
    fzero=fnew

    if (mpro%isMaster()) then
      ierr=0
      allocate(qvecs(n,m),hmat(m,m))
      hmat=(0.d0,0.d0)
      qvecs(:,1)=fnew/sqrt(sum(conjg(fnew)*fnew))
    end if

    do k=2,m
      if (mpro%isMaster()) fold=qvecs(:,k-1)
#ifdef PARALLEL
      call MPI_BCAST(fold, n, MPI_DOUBLE_COMPLEX, mpi_p_root, MPI_COMM_WORLD, ierr)
#endif
      !print *, mype, sum(fold)
      call next_iteration(n,fold,fnew)
      !print *, mype, sum(fnew)
      if (mpro%isMaster()) then
        qvecs(:,k)=fnew-fzero
        do j=1,k-1
          hmat(j,k-1)=sum(conjg(qvecs(:,j))*qvecs(:,k))
          qvecs(:,k)=qvecs(:,k)-hmat(j,k-1)*qvecs(:,j)
        end do
        hmat(k,k-1)=sqrt(sum(conjg(qvecs(:,k))*qvecs(:,k)))
        qvecs(:,k)=qvecs(:,k)/hmat(k,k-1)
      end if
    end do

    if (leigen .and. mpro%isMaster()) then
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

    integer, intent(in) :: m
    integer, intent(out) :: ierr, ngrow
    real(kind=real_kind), intent(in) :: tol
    complex(kind=complex_kind), dimension(m), intent(out) :: ritznum
    complex(kind=complex_kind), dimension(m,m), intent(in) :: hmat
    complex(kind=complex_kind), dimension(m,m), intent(out) :: eigh

    integer :: k,j,lwork,info
    complex(kind=complex_kind)   :: tmp
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
