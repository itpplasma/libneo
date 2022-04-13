!
  implicit none
!
  integer :: ntor,mpol,nlabel
  double precision :: flabel_min,flabel_max
  double complex, dimension(:,:,:), allocatable :: apsimn,athetmn
!
  open(1,form='unformatted',file='amn.dat.old')
  read (1) ntor,mpol,nlabel,flabel_min,flabel_max
  allocate(apsimn(-mpol:mpol,ntor,nlabel))
  allocate(athetmn(-mpol:mpol,ntor,nlabel))
  read (1) apsimn,athetmn
  close(1)
!
  apsimn=apsimn*0.1d0
  athetmn=athetmn*0.1d0
!
  open(1,form='unformatted',file='amn.dat')
  write (1) ntor,mpol,nlabel,flabel_min,flabel_max
  write (1) apsimn,athetmn
  close(1)
!
  end
