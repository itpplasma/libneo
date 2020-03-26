!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  implicit none
!
  character*1024 :: fluxdatapath
!
  integer :: ntor,mpol,nsqpsi
  integer :: m,n,k,i,nsqpsi_max,nsqpsi_rq
  double precision :: rtor,sqpsimin,sqpsimax,sqpsi,hsqpsi_rq,sqpsimin_rq,w
  double precision :: r,q
  double precision, dimension(:), allocatable :: rsmall,qsaf,psisurf
  double complex, dimension(:,:,:), allocatable :: apsimn,athetmn,br
!
  fluxdatapath='./'
!
  open(1,file='btor_rbig.dat')
  read (1,*) r,rtor
  close(1)
!
! Fourier ampitudes of the original field:
!
  open(1,form='unformatted',file=trim(fluxdatapath)//'/amn.dat')
  read (1) ntor,mpol,nsqpsi,sqpsimin,sqpsimax
  allocate(apsimn(-mpol:mpol,ntor,nsqpsi))
  allocate(athetmn(-mpol:mpol,ntor,nsqpsi))
  allocate(br(-mpol:mpol,ntor,nsqpsi))
  read (1) apsimn,athetmn
  close(1)
!
  nsqpsi_max=1000000
  allocate(rsmall(nsqpsi_max),qsaf(nsqpsi_max),psisurf(nsqpsi_max))
!
  nsqpsi_rq=0
  open(1,file='equil_r_q_psi.dat')
  read (1,*)
  read (1,*)
  read (1,*)
  do i=1,nsqpsi_max
    read (1,*,end=1) rsmall(i),qsaf(i),psisurf(i)
    nsqpsi_rq=nsqpsi_rq+1
  enddo
1 continue
  close(1)
  sqpsimin_rq=sqrt(psisurf(1))
  hsqpsi_rq=(sqrt(psisurf(nsqpsi_rq))-sqrt(psisurf(1)))/(nsqpsi_rq-1)
!
  print *,"m,n?"
  read *,m,n
  if(n.lt.0.or.n.gt.ntor.or.m.lt.-mpol.or.m.gt.mpol) then
    print *,'out of range, mpol = ',mpol,' ntor = ',ntor
  endif
!
  open(1,file="Br_jet.dat")
  do k=1,nsqpsi_rq
    w=(nsqpsi-1)*(rsmall(k)*qsaf(k)-sqpsimin)/(sqpsimax-sqpsimin)
    i=int(w)
    i=max(1,min(nsqpsi-1,i-1))
    w=w-i+1
    r=rsmall(k)
    q=qsaf(k)
    br(m,n,i)=(0.d0,1.d0)*n*(athetmn(m,n,i)*(1.d0-w)+athetmn(m,n,i+1)*w)/r/rtor
    write (1,*) r,q,real(br(m,n,i)),aimag(br(m,n,i)), &
                real(br(m,n,i))**2+aimag(br(m,n,i))**2
  enddo
  close(1)
!
  end
