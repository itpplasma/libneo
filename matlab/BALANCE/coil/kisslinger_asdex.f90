!
! program kisslinger
!
  implicit none
!
  integer :: L1i,iunit1,nr,np,nz,i,j,k
  double precision :: RT0,R0i,cbfi,BY0i,bf0
  double precision :: rmin,rmax,pmin,pmax,zmin,zmax,hrad,hphi,hzet
  double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
  double precision :: facbc
  common/cfacb/facbc
!
  facbc=1.d0
!
  iunit1=77
!
  open(iunit1,file='kisslinger.inp')
  read(iunit1,*) nr,np,nz
  read(iunit1,*) rmin, rmax
  read(iunit1,*) zmin, zmax
  close(iunit1)
!
!  call stevvo(RT0,R0i,L1i,cbfi,BY0i,bf0)
!
!  open(iunit1,file='MESH3D/stevvo_stuff.dat')
!  write(iunit1,*) RT0
!  write(iunit1,*) L1i
!  write(iunit1,*) bf0
!  write(iunit1,*) R0i
!  write(iunit1,*) cbfi
!  write(iunit1,*) BY0i
!  close(iunit1)
!
!  call set_curr
!
  L1i=1
!
  pmin=0.d0
  pmax=atan(1.d0)*8.d0/L1i
!
  hrad = (rmax - rmin)/(nr-1)  
  hphi = (pmax - pmin)/(np-1)
  hzet = (zmax - zmin)/(nz-1)
!
  open(iunit1,file='MESH3D/field.dat')
  write (iunit1,*) nr,np,nz,L1i
  write (iunit1,*) rmin,rmax
  write (iunit1,*) pmin,pmax
  write (iunit1,*) zmin,zmax
  do i=1,nr
    print *,i,'/',nr
    x(1)=rmin+hrad*(i-1)
    do j=1,np
      x(2)=pmin+hphi*(j-1)
      do k=1,nz
        x(3)=zmin+hzet*(k-1)
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        write (iunit1,*) hcovar(1)*bmod,hcovar(2)*bmod/x(1),hcovar(3)*bmod
      enddo
    enddo
  enddo
  close(iunit1)
!
  stop
  end
