!
  implicit none
!
  integer :: m0b,idummy,ns,is,imodes,nmodes
  double precision :: s,dummy
  integer,          dimension(:), allocatable :: mpol
  double precision, dimension(:), allocatable :: bmnc,bmns
!
  open(1,file='booz_test.bc')
  open(2,file='bmnc.dat')
  open(3,file='bmns.dat')
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*) m0b,idummy,ns
  nmodes = m0b +1
!
  allocate(mpol(nmodes),bmnc(nmodes),bmns(nmodes))
!
  read (1,*)
  read (1,*)
  read (1,*) s
  read (1,*)
  do imodes=1,nmodes
    read (1,*) mpol(imodes),idummy,dummy,dummy,dummy,dummy,dummy,dummy,bmnc(imodes),bmns(imodes)
  enddo
!
  write (2,*) '# s, Bmnc for m = ',mpol
  write (3,*) '# s, Bmns for m = ',mpol
  write (2,*) s,bmnc
  write (3,*) s,bmns
!
  do is=2,ns
    read (1,*)
    read (1,*)
    read (1,*) s
    read (1,*)
print *,is,s
    do imodes=1,nmodes
      read (1,*) mpol(imodes),idummy,dummy,dummy,dummy,dummy,dummy,dummy,bmnc(imodes),bmns(imodes)
    enddo
    write (2,*) s,bmnc
    write (3,*) s,bmns
  enddo
!
  close(1)
  close(2)
  close(3)
!
  end
