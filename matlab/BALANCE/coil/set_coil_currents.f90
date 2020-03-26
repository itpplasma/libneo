!
  implicit none
!
  integer, parameter :: ncoil=16
  double precision, parameter :: pi=3.14159265358979d0
  integer :: i
  double precision, dimension(ncoil) :: curr
!
  do i=1,8
!    curr(i)=sin(dfloat(i-1)*pi/2.d0)*4500.d0
    curr(i)=sin(dfloat(i-1)*pi/2.d0)*450.d0
!    curr(i+8)=curr(i)
    curr(i+8)=-curr(i)
  enddo
!
  open(1,file='cur_asd.dd')
  write(1,'(1000(1x,e17.10))') curr
  close(1)
!
  end
