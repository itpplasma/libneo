!
  subroutine gengrid(npoimin)
!
! Generates the grid for two types of boundary conditions at the outer
! boundary: iboutype=1 - fixed parameters, iboutype=2 - fixed fluxes
! At the inner boundary fixed fluxes = 0 are always assumed
!
  use grid_mod
!
  implicit none
!
  integer :: npoimin,ipoib,nder,ipb,ipe
  double precision :: hrmax,r,rnext,recnsp,rscale
  double precision, dimension(:),   allocatable :: x
  double precision, dimension(:,:), allocatable :: coef
!
  nbaleqs=4
!
  nder=1
  npoi_der=4
  allocate(x(npoi_der),coef(0:nder,npoi_der))
!
  hrmax=(rmax-rmin)/(npoimin+1)
!
  npoib=1
  r=rmin
!
  do while(r.lt.rmax)
    call recnsplit(r,recnsp)
    rnext=r+hrmax/recnsp
    call recnsplit(rnext,recnsp)
    r=0.5d0*(rnext+r+hrmax/recnsp)
    npoib=npoib+1
  enddo
  npoic=npoib-1
!
  allocate(rb(npoib),rc(npoic))
  allocate(Sb(npoib),Sc(npoic))
!
  r=rmin
  rb(1)=r
!
  do ipoib=2,npoib
    call recnsplit(r,recnsp)
    rnext=r+hrmax/recnsp
    call recnsplit(rnext,recnsp)
    r=0.5d0*(rnext+r+hrmax/recnsp)
    rb(ipoib)=r
    rc(ipoib-1)=0.5*(rb(ipoib-1)+rb(ipoib))
  enddo
!
  if(iboutype.eq.1) then
    rscale=(rmax-rmin)/(rc(npoic)-rmin)
  else
    rscale=(rmax-rmin)/(rb(npoib)-rmin)
  endif
  !rb=rmin+rscale*(rb-rmin)
  !rc=rmin+rscale*(rc-rmin)
!
  if(npoi_der.gt.npoic) then
    print *,'gengrid : not enough grid points for derivatives'
    stop
  endif
!
  allocate(deriv_coef(npoi_der,npoib),ipbeg(npoib),ipend(npoib))
  allocate(reint_coef(npoi_der,npoib))
!
  do ipoib=1,npoib
    ipb=ipoib-npoi_der/2
    ipe=ipb+npoi_der-1
    if(ipb.lt.1) then
      ipb=1
      ipe=ipb+npoi_der-1
    elseif(ipe.gt.npoic) then
      ipe=npoic
      ipb=ipe-npoi_der+1
    endif
    ipbeg(ipoib)=ipb
    ipend(ipoib)=ipe
    call plag_coeff(npoi_der,nder,rb(ipoib),rc(ipb:ipe),coef)
    deriv_coef(:,ipoib)=coef(1,:)
    reint_coef(:,ipoib)=coef(0,:)
  enddo
!
  deallocate(coef)
!
  allocate(params(nbaleqs,npoic),dot_params(nbaleqs,npoic))
  allocate(params_b(nbaleqs,npoib),ddr_params(nbaleqs,npoib),ddr_params_nl(nbaleqs,npoib))
  allocate(params_lin(nbaleqs,npoic),params_b_lin(nbaleqs,npoib))
  allocate(fluxes_dif(nbaleqs,npoib),fluxes_con(nbaleqs,npoib),fluxes_con_nl(nbaleqs,npoib))
!
  if(iboutype.eq.1) then
    neqset=nbaleqs*(npoic-1)
  else
    neqset=nbaleqs*npoic
  endif
  allocate(y(neqset),dery(neqset),dery_equisource(neqset))
  allocate(alpha(neqset*neqset))
  allocate(source_term(neqset))
!
  allocate(dae11(npoib),dae12(npoib),dae22(npoib))
  allocate(dai11(npoib),dai12(npoib),dai22(npoib))
  allocate(dni22(npoib),visca(npoib))
  allocate(dqle11(npoib),dqle12(npoib),dqle21(npoib),dqle22(npoib))
  allocate(dqli11(npoib),dqli12(npoib),dqli21(npoib),dqli22(npoib))
  allocate(de11(npoib),de12(npoib),de21(npoib),de22(npoib))
  allocate(di11(npoib),di12(npoib),di21(npoib),di22(npoib))
  allocate(polforce(npoib),qlheat_e(npoib),qlheat_i(npoib))
!
  dni22=0.d0
!
  allocate(cneo(npoib),gpp_av(npoib))
  allocate(qsafb(npoib),qsaf(npoic))
  allocate(sqg_bthet_overc(npoib),Ercov(npoib))
  allocate(Ercov_lin(npoib))
!
  return
  end subroutine gengrid
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine geomparprof
!
  use grid_mod
  use baseparam_mod
!
  implicit none
!
  integer :: ipoi
  double precision :: cneo_0,coullog,om_ci
!
  Sb=rb
  Sc=rc
  gpp_av=rtor**2
!
  coullog=15.d0
  om_ci=Z_i*e_charge*btor/(am*p_mass*c)
!  print *,'om_ci = ',om_ci
  cneo_0=1.32*4.d0*sqrt(pi)*Z_i**3*e_charge**4*coullog &
        /(3.d0*(am*p_mass)**1.5d0*om_ci**2)
  do ipoi=1,npoib
    qsafb(ipoi)=sum(qsaf(ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
    cneo(ipoi)=(rtor/rb(ipoi))**1.5d0*qsafb(ipoi)**2*cneo_0
    sqg_bthet_overc(ipoi)=btor*rb(ipoi)/qsafb(ipoi)/c
  enddo
!
  end subroutine geomparprof
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine recnsplit(r,recnsp)

  use resonances_mod
!  use grid_mod, only: gg_width, gg_factor, gg_r_res;
!
  implicit none;
!
  logical :: prop=.true.
  integer :: k
  double precision :: r, recnsp;
!
  if(prop) then
    prop=.false.
    call prepare_resonances
  endif
!
!  recnsp = 1.d0 + gg_factor*exp(-((r-gg_r_res)/gg_width)**2);
  recnsp = 1.d0
  do k=1,numres
    recnsp = recnsp + ampl_res(k)*exp(-((r-r_res(k))/width_res(k))**2)
  enddo
!
  return
  end subroutine recnsplit
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine prepare_resonances
!
  use resonances_mod
  use grid_mod, only: gg_width, gg_factor,r_resonant
  use mpi
!
  implicit none
!
  integer :: m,n,i,j,k,nr,jj,numres_orig,irank
  double precision :: qres,qmin,qmax
  complex :: a
  integer, dimension(:), allocatable :: m_a,n_a
  integer, dimension(:), allocatable :: m_aa,n_aa
  double precision, dimension(:), allocatable :: r,q
!
  iunit_res=157
!
  nr=0
  open(iunit_res,file='profiles/q.dat')
  do
    read(iunit_res,*,end=1)
    nr=nr+1
  enddo
1 continue
  close(iunit_res)
  allocate(r(nr),q(nr))
!
  open(iunit_res,file='profiles/q.dat')
  do i=1,nr
    read(iunit_res,*) r(i),q(i)
  enddo
  close(iunit_res)
  q=abs(q)
  qmin=minval(q)
  qmax=maxval(q)
!
  open(iunit_res,file='flre/antenna.in')
  read(iunit_res,*)
  read(iunit_res,*)
  read(iunit_res,*)
  read(iunit_res,*)
  read(iunit_res,*)
  read(iunit_res,*) numres
  close(iunit_res)
!
  allocate(r_res(numres),width_res(numres),ampl_res(numres))
  allocate(r_resonant(numres))
  width_res=gg_width
  ampl_res=gg_factor
  allocate(m_a(numres),n_a(numres))
  allocate(m_aa(numres),n_aa(numres))
  numres_orig=numres
!
  open(iunit_res,file='flre/modes.in')
  k=1
  jj=1
  read(iunit_res,*) a
  m=-abs(nint(real(a)))
  n=abs(nint(imag(a)))
  m_a(k)=m
  n_a(k)=n
  m_aa(jj)=m
  n_aa(jj)=n
  r_res(k)=abs(dfloat(m)/dfloat(n))
  outer: do i=2,numres
    read(iunit_res,*) a
    m=-abs(nint(real(a)))
    n=abs(nint(imag(a)))
    qres=abs(dfloat(m)/dfloat(n))
!check for existence of resonant point:
    if(qres.lt.qmin.or.qres.gt.qmax) cycle
!check for repeated resonance radii:
    jj=jj+1
    m_aa(jj)=m
    n_aa(jj)=n
    do j=1,k
      if(m*n_a(j)-n*m_a(j).eq.0) cycle outer
    enddo
    k=k+1
    m_a(k)=m
    n_a(k)=n
    r_res(k)=qres
  enddo outer
  close(iunit_res)
  numres=k
!
  do i=1,numres
    qres=r_res(i)
    do j=2,nr
      if(qres.gt.q(j-1).and.qres.lt.q(j)) then
        r_res(i)=(r(j-1)*(q(j)-qres)+r(j)*(qres-q(j-1)))/(q(j)-q(j-1))
        exit
      endif
    enddo
  enddo
  if (irank .eq. 0 ) then
    print *,'gengrid: number of resonance points = ',numres
    print *,'resonant radii: ',r_res(1:numres)
  endif
!
  do i=1,numres_orig
    do j=1,numres
      if(m_aa(i)*n_a(j)-n_aa(i)*m_a(j).eq.0) r_resonant(i)=r_res(j)
    enddo
  enddo
  if (irank .eq. 0 ) then
    print *,'resonant radii for all modes: ',r_resonant(1:numres_orig)
  endif
  
!
  deallocate(m_a,n_a,r,q)
!
  end subroutine prepare_resonances
