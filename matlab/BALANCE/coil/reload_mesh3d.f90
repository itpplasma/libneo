!
! program kisslinger
!
  implicit none
!
  integer :: L1i,iunit1,nr,np,nz,i,j,k,iunit0
  double precision :: RT0,R0i,cbfi,BY0i,bf0
  double precision :: rmin,rmax,pmin,pmax,zmin,zmax,hrad,hphi,hzet
  double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
  double precision :: facbc
  common/cfacb/facbc
!
  iunit0=76
  iunit1=77
!
  open(iunit0,file='field.dat.old')
  open(iunit1,file='field.dat')
  read (iunit0,*) nr,np,nz,L1i
  write (iunit1,*) nr,np,nz,L1i
  read (iunit0,*) rmin,rmax
  write (iunit1,*) rmin,rmax
  read (iunit0,*) pmin,pmax
  write (iunit1,*) pmin,pmax
  read (iunit0,*) zmin,zmax
  write (iunit1,*) zmin,zmax
  do i=1,nr
    print *,i,'/',nr
    do j=1,np
      do k=1,nz
        read (iunit0,*) hcovar(1),hcovar(2),hcovar(3)
        write (iunit1,*) hcovar(1)*0.1d0,hcovar(2)*0.1d0,hcovar(3)*0.1d0
      enddo
    enddo
  enddo
  close(iunit0)
  close(iunit1)
!
  stop
  end
