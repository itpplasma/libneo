!
  implicit none
  double precision :: r,Vz,scale
!
  print *,'scaling factor?'
  read *,scale
!
  open(1,file='Vz_org.dat')
  open(2,file='Vz.dat')
  do
    read(1,*,end=1) r,Vz
    write(2,*) r,Vz*scale
  enddo
1 continue
  close(1)
  close(2)
  open(1,file='scaling_factor.dat')
  write(1,*) scale
  close(1)
  end
