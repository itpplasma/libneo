program efit_to_boozer

  implicit none

  integer :: i,k,is,it,nsurf,nt,mpol,ntor,nmodes,m,iunit,modfactor,nplotsurf
  double precision :: s,phi,hs,htheta,rho_tor,hrho,w,s_plot,s_old

  integer,          dimension(:), allocatable :: m_arr,n_arr
  double precision, dimension(:), allocatable :: Rmn_c,Rmn_s,Zmn_c,Zmn_s,almn_c,almn_s,Bmn_c,Bmn_s
  double precision, dimension(:), allocatable :: theta,Rnew,Znew,Rold,Zold

  nplotsurf=50
  hrho=1.d0/dble(nplotsurf)
!  phi=0.d0
  print *,'Enter phi in pi units:'
  read *,phi
  phi=phi*atan(1.d0)*4.d0

  iunit=71
  open(iunit,file='boozfile_to_plot.bc')
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*) mpol,ntor,nsurf

  nmodes=(mpol+1)*(2*ntor+1)
  allocate(m_arr(nmodes),n_arr(nmodes),Rmn_c(nmodes),Rmn_s(nmodes),Zmn_c(nmodes),Zmn_s(nmodes))
  modfactor=30
  nt=mpol*modfactor
  htheta=atan(1.d0)*8.d0/dble(nt)
  allocate(theta(0:nt),Rnew(0:nt),Znew(0:nt),Rold(0:nt),Zold(0:nt))

  do i=0,nt
    theta(i)=htheta*dble(i)
  enddo

  Rnew=0.d0
  Znew=0.d0
  rho_tor=hrho
  s_plot=rho_tor**2
  s=0.d0
  open(1001,file='surfplot.dat')

  do is=1,nsurf
    s_old=s
    read(iunit,*)
    read(iunit,*)
    read(iunit,*) s
    read(iunit,*)
    do k=1,nmodes
      read(iunit,*) m_arr(k),n_arr(k),Rmn_c(k),Rmn_s(k),Zmn_c(k),Zmn_s(k)
    enddo

    Rold=Rnew
    Zold=Znew
    do i=0,nt
      Rnew(i)=sum(Rmn_c*cos(m_arr*theta(i)-n_arr*phi))+sum(Rmn_s*sin(m_arr*theta(i)-n_arr*phi))
      Znew(i)=sum(Zmn_c*cos(m_arr*theta(i)-n_arr*phi))+sum(Zmn_s*sin(m_arr*theta(i)-n_arr*phi))
    enddo

    if(s.gt.s_plot) then
      w=(s_plot-s_old)/(s-s_old)
      do i=0,nt
        write(1001,*) Rnew(i)*w+Rold(i)*(1.d0-w),Znew(i)*w+Zold(i)*(1.d0-w)
      enddo
      write(1001,*) ' '
      rho_tor=rho_tor+hrho
      s_plot=rho_tor**2
    endif

  enddo

  close(iunit)
  close(1001)

  deallocate(m_arr,n_arr,Rmn_c,Rmn_s,Zmn_c,Zmn_s,theta,Rnew,Znew,Rold,Zold)

end program efit_to_boozer
