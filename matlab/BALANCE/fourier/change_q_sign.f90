!
  implicit none
!
  integer :: mpol,nrad,nstep,nsqpsi,i
  double precision :: r,q,psipol,phitor,dphidpsi,rsmall,volume
!
!
  open(1,file='fouriermodes.inp')
  read (1,*) mpol     !number of poloidal modes
  read (1,*) nstep    !number of integration steps
  read (1,*) nsqpsi   !grid size over radial variable for 1D quantities
  close(1)
!
  open(1,file='equil_r_q_psi.dat')
  open(2,file='equil_r_q_psi.dat.qneg')
  read (1,*)
  read (1,*)
  read (1,*)
  do i=1,nsqpsi
    read (1,*) r,q,psipol,phitor,dphidpsi,rsmall,volume
    write (2,*) r,-q,psipol,phitor,dphidpsi,rsmall,volume
  enddo
  close(1)
  close(2)
!
  stop
  end
