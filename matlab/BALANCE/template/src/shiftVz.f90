!
  implicit none
  double precision :: r,Vz,scale,Q_spline,Q_out,rsepar = 67.d0
  logical :: firstpoint
!
  print *,'scaling factor?'
  read *,scale
!
  open(1,file='Vz_org.dat')
  open(2,file='Vz.dat')
  firstpoint = .true.
  do
    read(1,*,end=1) r,Vz
    if (firstpoint) then
        firstpoint = .false.
        Q_spline = Vz
    endif
    Q_out = Q_spline/(1.0d0 + exp((r-rsepar)/0.3d0))
    write(2,*) r,Vz+Q_out*scale
  enddo
1 continue
  close(1)
  close(2)
  open(1,file='scaling_factor.dat')
  write(1,*) scale
  close(1)
  end
