!
  module bdivfree_mod
    integer :: nr,nz,ntor,icp
    integer, dimension(:,:), allocatable :: ipoint
    double precision :: rmin,zmin,hr,hz,pmin,pfac
    double precision, dimension(:),       allocatable :: rpoi,zpoi
    double precision, dimension(:,:,:),   allocatable :: apav,rbpav_coef
    double precision, dimension(:,:,:,:), allocatable :: aznre,aznim,arnre,arnim
  end module bdivfree_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module theta_rz_mod
    integer :: icall=0
    integer :: nsqp,nlab,nthe,icp_pt,n_subgrid_fl
    integer,      dimension(:),     allocatable :: ind_subgrid_fl
    integer,      dimension(:,:),   allocatable :: ipoint_pt
    double precision :: h_subgrid_fl
    real(kind=8) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
    real(kind=8), dimension(:,:),   allocatable :: spllabel
    real(kind=8), dimension(:,:,:), allocatable :: splthet
    real(kind=8), dimension(:),     allocatable :: sqpsi,flab,theqt
    real(kind=8), dimension(:),     allocatable :: sqpsi_lab
  end module theta_rz_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module new_amn_mod
! Fourier ampitudes of the original field:
    logical :: prop=.true.
    integer, parameter :: nhalfband=2, ncoefs=2*nhalfband+2
    integer :: mpol,ntor_act,nlabel_tot,n_subgrid,isw_gauge
    double precision :: h_subgrid
    integer,          dimension(:),     allocatable :: ntor_arr,ind_subgrid
    integer,          dimension(:),     allocatable :: ibeg_arr,iend_arr
    double precision, dimension(:),     allocatable :: flabel_arr
    double precision, dimension(:,:),   allocatable :: splsqggpp
    double complex,   dimension(:,:,:), allocatable :: splaalp,splatet
    double complex,   dimension(:,:,:), allocatable :: splaeta,splaphi
    double complex,   dimension(:), allocatable :: amalp,   amtet,     &
                                                   amalp_s, amtet_s,   &
                                                   amalp_ss,amtet_ss
    double complex,   dimension(:), allocatable :: ameta,   amphi,     &
                                                   ameta_s, amphi_s,   &
                                                   ameta_ss,amphi_ss
    double complex,   dimension(:), allocatable :: expthe,expphi
    double complex,   dimension(:), allocatable :: dthe,ddthe
  end module new_amn_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module extract_fluxcoord_mod
    integer :: load_extract_fluxcoord=1
    integer :: nphinorm
    double precision :: psif_extract,theta_extract,psifmin,hpsif
    double precision :: psifmax,phifmax,sigcos
    double precision, dimension(:), allocatable :: phinorm_arr
  end module extract_fluxcoord_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module metrics_spectra_mod
!
    logical :: load_metrics=.true.
!
    integer :: n_fix=0
    integer :: nlabel_gik,mpol_amn,ntor_amn
!
    double precision :: flabel,sqrtg_over_gphiphi,sqrtpsi_lab,dsqrtpsi_lab
!
    double precision, dimension(:),     allocatable :: sqpsi_lab,gik_lab
    double precision, dimension(:),     allocatable :: sqg_over_gpp2
!
    double precision, dimension(:,:),   allocatable :: ai_sqg_over_gpp2
    double precision, dimension(:,:),   allocatable :: ai_sqpsi_lab
!
    double complex,   dimension(:),     allocatable :: one_over_gpsipsi
    double complex,   dimension(:),     allocatable :: one_over_sqrtggpsipsi
    double complex,   dimension(:),     allocatable :: gpsitheta_over_gpsipsi
    double complex,   dimension(:),     allocatable :: athetm_vac
    double complex,   dimension(:),     allocatable :: apsim_vac
!
    double complex,   dimension(:,:),   allocatable :: one_over_gpp
    double complex,   dimension(:,:),   allocatable :: gpt_over_gpp
    double complex,   dimension(:,:),   allocatable :: one_over_sqggpp
!
    double complex,   dimension(:,:,:), allocatable :: ai_one_over_gpp
    double complex,   dimension(:,:,:), allocatable :: ai_gpt_over_gpp
    double complex,   dimension(:,:,:), allocatable :: ai_one_over_sqggpp
    double complex,   dimension(:,:,:), allocatable :: athetmn_vac
    double complex,   dimension(:,:,:), allocatable :: apsimn_vac
    double complex,   dimension(:,:,:), allocatable :: ai_athetm
    double complex,   dimension(:,:,:), allocatable :: ai_apsim
!
  end module metrics_spectra_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module rhs_coulgauge_mod
    integer :: neq_amn,nmodes,ncomps,ntor_one
    double complex,   dimension(:),   allocatable :: athetm,alpham
    double complex,   dimension(:),   allocatable :: dathetm,dalpham
    double complex,   dimension(:),   allocatable :: etam,chim
    double complex,   dimension(:),   allocatable :: detam,dchim
  end module rhs_coulgauge_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module start_end_coulgauge_mod
    double complex, dimension(:,:), allocatable :: athetmn_beg
    double complex, dimension(:,:), allocatable :: athetmn_end
    double complex, dimension(:,:), allocatable :: apsimn_beg
    double complex, dimension(:,:), allocatable :: apsimn_end
  end module start_end_coulgauge_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine vector_potentials(nr_in,np_in,nz_in,ntor_in,      &
             rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in,  &
             br,bp,bz)
!
  use bdivfree_mod
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nr_in,np_in,nz_in,ntor_in,ip,np,n,ir,iz
  integer, dimension(:), allocatable :: imi,ima,jmi,jma
!
  integer :: nashli_rukami
  integer :: irmin, irmax, i,j
  double precision, dimension(4), parameter :: weight=(/-1., 13., 13., -1./)/24.
!
  double precision :: rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in
  double precision :: hp,r,rm,zm,sumbz,hrm1,hzm1
  double precision, dimension(nr_in,np_in,nz_in)  :: br,bp,bz
  double precision, dimension(:),     allocatable :: dummy
  double precision, dimension(:,:),   allocatable :: a_re, a_im, rbpav_dummy
  double precision, dimension(:,:),   allocatable :: brm,bpm,bzm
!
  double complex :: four_ampl
  double complex, dimension(:,:), allocatable :: expon
!
  integer, parameter :: mp=4 ! power of Lagrange polynomial =3
  integer,          dimension(mp)    :: indx,indy
  double precision, dimension(mp)    :: xp,yp
  double precision, dimension(mp,mp) :: fp
!
  nr=nr_in
  nz=nz_in
  np=np_in-1
  ntor=ntor_in
  nashli_rukami=(nr_in+1)/2
!
  rmin=rmin_in
  zmin=zmin_in
  hr=(rmax_in-rmin_in)/(nr-1)
  hz=(zmax_in-zmin_in)/(nz-1)
  hp=2.d0*pi/np
  pmin=pmin_in
  pfac = dble(nint(2.d0*pi/(pmax_in-pmin_in)))
!
  allocate(expon(np,ntor),a_re(nr,nz),a_im(nr,nz),rbpav_dummy(nr,nz))
  allocate(imi(nz),ima(nz),jmi(nr),jma(nr), dummy(nr))
  allocate(rpoi(nr),zpoi(nz))
  allocate(brm(nr,nz),bpm(nr,nz),bzm(nr,nz))
!
  imi=1
  ima=nr
  jmi=1
  jma=nz
  do ir=1,nr
    rpoi(ir)=rmin+hr*(ir-1)
  enddo
  do iz=1,nz
    zpoi(iz)=zmin+hz*(iz-1)
  enddo
!
! Truncation of data outside the limiting convex:
!
  hrm1=1.d0/hr
  hzm1=1.d0/hz
  do ip=1,np
    do ir=1,nr
      do iz=1,nz
        call stretch_coords(rpoi(ir),zpoi(iz),rm,zm)
        call indef_bdf(rm,rmin,hrm1,nr,indx)
        call indef_bdf(zm,zmin,hzm1,nz,indy)
!
        do i=1,mp
          xp(i) = rpoi(indx(i))
          yp(i) = zpoi(indy(i))
        enddo
!
        do j=1,mp
          do i=1,mp
            fp(i,j) = Br(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,brm(ir,iz))
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bp(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bpm(ir,iz))
        bpm(ir,iz)=bpm(ir,iz)*rm/rpoi(ir)
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bz(indx(i),ip,indy(j))
          enddo
        enddo
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bzm(ir,iz))
!
      enddo
    enddo
    Br(:,ip,:)=brm
    Bp(:,ip,:)=bpm
    Bz(:,ip,:)=bzm
  enddo
!
! End of data truncation
!
  allocate(ipoint(nr,nz))
  icp=nr*nz
  allocate(aznre(6,6,icp,ntor),aznim(6,6,icp,ntor))
  allocate(arnre(6,6,icp,ntor),arnim(6,6,icp,ntor))
  allocate(apav(6,6,icp),rbpav_coef(6,6,icp))
!
  do n=1,ntor
    do ip=1,np
      expon(ip,n)=exp(dcmplx(0.d0,-n*(ip-1)*hp))/np
    enddo
  enddo
!
  do n=1,ntor
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(br(ir,1:np,iz)*expon(:,n))*dcmplx(0.d0,-r/(n*pfac))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      enddo
    enddo
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,aznre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,aznim(:,:,:,n),ipoint)
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(bz(ir,1:np,iz)*expon(:,n))*dcmplx(0.d0,r/(n*pfac))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      enddo
    enddo
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,arnre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,arnim(:,:,:,n),ipoint)
  enddo
!
  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        dummy(ir) = sum(bz(ir,1:np,iz))*hr*r/np
     enddo
     a_re(nashli_rukami,iz) = 0.
     sumbz=0.d0
     do ir=nashli_rukami+1,nr
        irmax = min(ir+1,nr)
        irmin = irmax - 3
        sumbz = sumbz + sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     enddo
     sumbz=0.d0
     do ir=nashli_rukami-1,1,-1
        irmin = max(ir-1,1)
        irmax = irmin + 3
        sumbz = sumbz - sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     enddo
  enddo
!
  call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,apav,ipoint)
!
  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        rbpav_dummy(ir,iz) = r*sum(bp(ir,1:np,iz))/np
     enddo
  enddo
!
  call s2dcut(nr,nz,hr,hz,rbpav_dummy,imi,ima,jmi,jma,icp,rbpav_coef,ipoint)
!
  deallocate(expon,a_re,a_im,rbpav_dummy,imi,ima,jmi,jma,dummy,brm,bpm,bzm)
!
  return
  end subroutine vector_potentials
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
  use bdivfree_mod
!
  implicit none
!
  integer :: n,ierr
  double precision :: r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: f,fr,fz,frr,frz,fzz
  double precision :: g,gr,gz,grr,grz,gzz
  double precision :: delbr,delbz,delbp
  double precision :: deldBrdR,deldBrdp,deldBrdZ
  double precision :: deldBpdR,deldBpdp,deldBpdZ
  double precision :: deldBzdR,deldBzdp,deldBzdZ
  double complex :: expon,anr,anz,anr_r,anr_z,anz_r,anz_z
  double complex :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz
!
!
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,rbpav_coef,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)
  Bp = f/r
  dBpdR = fr/r - Bp/r
  dBpdZ = fz/r
  dBpdp = 0.d0
!
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,apav,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)
!
  Br=-fz/r
  Bz=fr/r
  dBrdR=fz/r**2-frz/r
  dBrdZ=-fzz/r
  dBzdR=-fr/r**2+frr/r
  dBzdZ=frz/r
  dBrdp=0.d0
  dBzdp=0.d0
!
  do n=1,ntor
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anr=dcmplx(f,g)
    anr_r=dcmplx(fr,gr)
    anr_z=dcmplx(fz,gz)
    anr_rr=dcmplx(frr,grr)
    anr_rz=dcmplx(frz,grz)
    anr_zz=dcmplx(fzz,gzz)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anz=dcmplx(f,g)
    anz_r=dcmplx(fr,gr)
    anz_z=dcmplx(fz,gz)
    anz_rr=dcmplx(frr,grr)
    anz_rz=dcmplx(frz,grz)
    anz_zz=dcmplx(fzz,gzz)
!
    expon=exp(dcmplx(0.d0,n*pfac*phi))
    delbr=2.d0*dble(dcmplx(0.d0,pfac)*n*anz*expon/r)
    delbz=-2.d0*dble(dcmplx(0.d0,pfac)*n*anr*expon/r)
    delbp=2.d0*dble((anr_z-anz_r)*expon)
    deldBrdR=-delbr/r+2.d0*dble(dcmplx(0.d0,pfac)*n*anz_r*expon/r)
    deldBrdZ=2.d0*dble(dcmplx(0.d0,pfac)*n*anz_z*expon/r)
    deldBrdp=-2.d0*dble((pfac*n)**2*anz*expon/r)
    deldBzdR=-delbz/r-2.d0*dble(dcmplx(0.d0,pfac)*n*anr_r*expon/r)
    deldBzdZ=-2.d0*dble(dcmplx(0.d0,pfac)*n*anr_z*expon/r)
    deldBzdp=2.d0*dble((pfac*n)**2*anr*expon/r)
    deldBpdR=2.d0*dble((anr_rz-anz_rr)*expon)
    deldBpdZ=2.d0*dble((anr_zz-anz_rz)*expon)
    deldBpdp=2.d0*dble(dcmplx(0.d0,pfac)*n*(anr_z-anz_r)*expon)
!
    br=br+delbr
    bz=bz+delbz
    bp=bp+delbp
    dBrdR=dBrdR+deldBrdR
    dBrdZ=dBrdZ+deldBrdZ
    dBrdp=dBrdp+deldBrdp
    dBzdR=dBzdR+deldBzdR
    dBzdZ=dBzdZ+deldBzdZ
    dBzdp=dBzdp+deldBzdp
    dBpdR=dBpdR+deldBpdR
    dBpdZ=dBpdZ+deldBpdZ
    dBpdp=dBpdp+deldBpdp
!
  enddo
!
  return
  end subroutine field_divfree
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine indef_bdf(u,umin,dum1,nup,indu)
! defines interval for 1D interpolation on uniform mesh, normally
! looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    u - coordinate of a point (to be interpolated)
!    umin - minimal value of u
!    dum1 = 1./h reciprocal to length of mesh interval
!    nup - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points
!
! the power 3 of polinomial is fixed strictly:
!
      implicit double precision (a-h,o-z)
!
      parameter(mp=4)
      integer indu(mp)

      indu(1) = int((u-umin)*dum1)
      if( indu(1) .le. 0 ) indu(1) = 1
      indu(mp) = indu(1) + mp - 1
      if( indu(mp) .gt. nup ) then
         indu(mp) = nup
         indu(1) = indu(mp) - mp + 1
      endif
      do i=2,mp-1
         indu(i) = indu(i-1) + 1
      enddo

      return
      end
!---------------------------------------------------------------------
      subroutine indsmp_bdf(index,nup,indu)
! defines interval for 1D interpolation on uniform mesh
! by known index.
! Normally looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    index - number of a cell on the mesh
!    nup - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points

! the power 3 of polinomial is fixed strictly:
      parameter(mp=4)
      integer indu(mp)

      indu(1) = index - 1
      if( indu(1) .le. 0 ) indu(1) = 1
      indu(mp) = indu(1) + mp - 1
      if( indu(mp) .gt. nup ) then
         indu(mp) = nup
         indu(1) = indu(mp) - mp + 1
      endif
      do i=2,mp-1
         indu(i) = indu(i-1) + 1
      enddo

      return
      end
!---------------------------------------------------------------------
      subroutine plag2d_bdf(x,y,fp,dxm1,dym1,xp,yp,polyl2d)
!
      implicit double precision (a-h,o-z)
!
! 2D interpolation by means of Lagrange polynomial
! the power 3 is fixed strictly:
      parameter(mp=4)
! uniform mesh (increasingly ordered) in all dimensions is implied
!
! Input parameters:
!   x,y - coordinates of the point for interpolation
!   dxm1,dym1 - inverse steps in each direction
!   xp,yp - vertices of stencil
!
! Output parameters:
! polyl2d - polynomial itself
      dimension cx(mp),cy(mp),fp(mp,mp),xp(mp),yp(mp)
!
      call coefs_bdf(x,xp,dxm1,cx)
      call coefs_bdf(y,yp,dym1,cy)
!
      polyl2d = 0.d0
      do j=1,mp
        do i=1,mp
          polyl2d = polyl2d + fp(i,j)*cx(i)*cy(j)
        enddo
      enddo
!
      return
      end
!---------------------------------------------------------------------
      subroutine coefs_bdf(u,up,dum1,cu)
!
      implicit double precision (a-h,o-z)
!
      parameter(mp=4)
      dimension up(mp),cu(mp)
      data one6/0.16666666666667d0/
      du3 = dum1**3
      cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-one6*du3)
      cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5d0*du3)
      cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5d0*du3)
      cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (one6*du3)
      return
      end
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine s2dring(nx,ny,hx,hy,f,icount,spl,ipoint)
!
! Calculates coefficients of a 2D spline for a ring domain
! (periodic over y variable)
! equidistant mesh, but hx must not be = hy
!
!  Input parameters:
!                    nx           - horizontal size of the mesh (over x)
!                    ny           - vertical size of the mesh (over y)
!                    hx           - step of the mesh over x
!                    hy           - step of the mesh over y
!                    f(i,j)       - array of values to be interpolated
!                                   (i = 1, ..., nx; j = 1, ..., ny).
!
!                                   For the case of non-rectangular domain:
!                                   numbers of the mesh points which
!                                   correspond to the boundaries of the
!                                   interpolation region:
!                    icount       - maximum number of entries in spl
! Output parameters:
!                    spl(l,m,k)   - spline coefficients (i,j = 1, ... , n;
!                    ipoint(i,j)    l,m = 1, ..., 4; i,j - numbers of the
!                                   mesh point in horizontal and vertical
!                                   direction (over x and over y), l,m -
!                                   the numbers of expansion power over x
!                                   and y (~ dx**(l-1)*dy**(m-1) ))
!                                   ipoint(i,j) contains the pointer to k
!
  use spl_three_to_five_mod, only: spl_five_per
  implicit double precision (a-h,o-z)
!
  dimension f(nx,ny),spl(6,6,icount),ipoint(nx,ny)
!
  integer,          dimension(:), allocatable :: imi,ima,jmi,jma
  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi
!
  nmax=max(nx,ny)
!
  allocate( ai(nmax),bi(nmax),ci(nmax),di(nmax),ei(nmax),fi(nmax) )
  allocate(imi(ny),ima(ny),jmi(nx),jma(nx))
!
  imi=1
  ima=nx
  jmi=1
  jma=ny
!
  spl=0.d0
  ipoint=-1
!
!  spline along Y-axis
!
  ic = 0
  do i=1,nx
    if(jmi(i).gt.0) then
      nsi=jma(i)-jmi(i)+1
      do j=jmi(i),jma(i)
        ai(j-jmi(i)+1)=f(i,j)
      enddo
      call spl_five_per(nsi,hy,ai,bi,ci,di,ei,fi)
      do j=jmi(i),jma(i)
        jj=j-jmi(i)+1
        ic = ic+1
        ipoint(i,j)=ic
        spl(1,1,ic)=ai(jj)
        spl(1,2,ic)=bi(jj)
        spl(1,3,ic)=ci(jj)
        spl(1,4,ic)=di(jj)
        spl(1,5,ic)=ei(jj)
        spl(1,6,ic)=fi(jj)
      enddo
    endif
  enddo
!
  if (ic .ne. icount) then
    write (6,*) 'Warning, ic, icount:  ',ic,icount
  endif
!
!  spline along X-axis
!
  do j=1,ny
    if(imi(j).gt.0) then
      nsi=ima(j)-imi(j)+1
      do l=1,6
        do i=imi(j),ima(j)
          ai(i-imi(j)+1)=spl(1,l,ipoint(i,j))
        enddo
        call spl_five_reg(nsi,hx,ai,bi,ci,di,ei,fi)
        do i=imi(j),ima(j)
          ii=i-imi(j)+1
          spl(2,l,ipoint(i,j))=bi(ii)
          spl(3,l,ipoint(i,j))=ci(ii)
          spl(4,l,ipoint(i,j))=di(ii)
          spl(5,l,ipoint(i,j))=ei(ii)
          spl(6,l,ipoint(i,j))=fi(ii)
        enddo
      enddo
    endif
  enddo
!
  deallocate( ai,bi,ci,di,ei,fi )
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
MODULE polleg_mod

  IMPLICIT NONE

  INTEGER :: dummy

CONTAINS
!
  SUBROUTINE polleg(n,coefleg)
!
! Computes coefficients of Legendre polynomials of orders from 0 to n

!
! Input parameters:
!           Formal: n            - maximum order of Legendre polynomials
! Output parameters:
!           Formal: coefleg(i,j) - j-th coefficient of Legendre polynomial
!                                  of the order i
!
  INTEGER :: n,i,j
!
  DOUBLE PRECISION :: frontfac,rearfac
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coefleg
!
  IF(ALLOCATED(coefleg)) DEALLOCATE(coefleg)
  ALLOCATE(coefleg(0:n,0:n))
!
  coefleg=0.d0
  coefleg(0,0)=1.d0
  coefleg(1,1)=1.d0
  frontfac=1.d0
!
  DO i=2,n
    frontfac=frontfac*(2.d0-1.d0/DBLE(i))
    rearfac=frontfac
    coefleg(i,i)=rearfac
    DO j=i-2,0,-2
      rearfac=-rearfac*DBLE(j+1)*DBLE(j+2)/DBLE(i-j)/DBLE(i+j+1)
      coefleg(i,j)=rearfac
    ENDDO
  ENDDO
!
  RETURN
END SUBROUTINE polleg
!
  SUBROUTINE binomial(n,coefbin)
!
! Computes binomial coefficients of orders from 0 to n
!
! Input parameters:
!           Formal: n            - maximum power of the binom
! Output parameters:
!           Formal: coefbin(i,j) - j-th coefficient of the binom
!
  INTEGER :: n,i,j
!
  DOUBLE PRECISION :: frontfac,factforw,factbackw
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coefbin
!
  IF(ALLOCATED(coefbin)) DEALLOCATE(coefbin)
  ALLOCATE(coefbin(0:n,0:n))
!
  coefbin=0.d0
  coefbin(0,0)=1.d0
  frontfac=1.d0
!
  DO i=1,n
    frontfac=frontfac*DBLE(i)
    factforw=1.d0
    factbackw=frontfac*DBLE(i+1)
    DO j=0,i
      IF(j.GT.0) factforw=factforw*DBLE(j)
      factbackw=factbackw/DBLE(i-j+1)
      coefbin(i,j)=frontfac/factforw/factbackw
    ENDDO
  ENDDO
!
  RETURN
END SUBROUTINE binomial
END MODULE
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module binsrc_mod
contains
!
  subroutine binsrc(p,nmin,nmax,xi,i)
!
! Finds the index  i  of the array of increasing numbers   p  with dimension  n
! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.
!
  implicit none
!
  integer                                :: n,nmin,nmax,i,imin,imax,k
  double precision                       :: xi
  double precision, dimension(nmin:nmax) :: p
!
  imin=nmin
  imax=nmax
  n=nmax-nmin
!
  do k=1,n
    i=(imax-imin)/2+imin
    if(p(i).gt.xi) then
      imax=i
    else
      imin=i
    endif
    if(imax.eq.imin+1) exit
  enddo
!
  i=imax
!
  return
  end subroutine binsrc
end module binsrc_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
MODULE oddorderspline_mod
CONTAINS
!
  SUBROUTINE oddorderspline(nshift,nbands,nfun,eta,f,ai,ierr)
!
! Computes natural spline coefficients for the arbitrary odd order spline
! Spline order is npower=2*nshift+1
! nbands                      - number of points
! nfun                        - number of splined functions
! eta(1:nbands)               - argument
! f(1:nbands,nfun)            - splined functions
! ai(0:npowers,1:nbands,nfun) - spline coefficients
!
  USE polleg_mod
!
  IMPLICIT NONE
!
  INTEGER :: nshift,nbands,nfun,npowers
  INTEGER :: i,k,ierr,kk,ndim,ind1,ind2,info
  INTEGER :: kl,ku,ldab,indbeg
  INTEGER,          DIMENSION(:), ALLOCATABLE     :: ipiv
  DOUBLE PRECISION                    :: deleta
  DOUBLE PRECISION, DIMENSION(1:nbands)           :: eta
  DOUBLE PRECISION, DIMENSION(1:nbands,nfun)      :: f
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: powdeleta
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: amat,bvec,amatch,coefbin
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: ails
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: ai
!
  npowers=2*nshift+1
  ndim=(npowers+1)*(nbands-1)+2*nshift+1
  kl=3*nshift+1
  ku=nshift+1
  ldab=2*kl+ku+1
  indbeg=kl+ku+1
!
  IF(npowers+1.GT.nbands) THEN
    ierr=1
    PRINT *,                                                                  &
         'oddorderspline: number of points is not enough for this spline order'
    RETURN
  ENDIF
!
  IF(ALLOCATED(ai)) DEALLOCATE(ai)
  ALLOCATE(ai(0:npowers,1:nbands,nfun),ails(npowers+1,nfun))
!
  ALLOCATE(amat(ldab,ndim),bvec(ndim,nfun),ipiv(ndim))
  ALLOCATE(powdeleta(0:npowers),amatch(0:npowers,0:npowers), &
           coefbin(0:npowers,0:npowers))
!
  amatch=0.d0
!
  CALL binomial(npowers,coefbin)
!
  amat=0.d0
  bvec=0.d0
!
  DO i=1,nbands-1
!
    ind1=(npowers+1)*(i-1)+1+nshift
    bvec(ind1,:)=f(i,:)
    ind2=(npowers+1)*(i-1)+1
!    amat(ind1,ind2)=1.d0
    amat(indbeg+ind1-ind2,ind2)=1.d0
!
! Spline continuity conditions:
!
    deleta=eta(i+1)-eta(i)
    powdeleta(0)=1.d0
    DO k=1,npowers
      powdeleta(k)=powdeleta(k-1)*deleta
    ENDDO
!
    DO k=0,npowers
      amatch(k,k:npowers)=powdeleta(0:npowers-k)*coefbin(k:npowers,k)
    ENDDO
!
    DO k=0,npowers-1
      ind1=(npowers+1)*(i-1)+k+2+nshift
      ind2=(npowers+1)*i+k+1
!      amat(ind1,ind2)=-1
      amat(indbeg+ind1-ind2,ind2)=-1.d0
      DO kk=0,npowers
        ind2=(npowers+1)*(i-1)+kk+1
!        amat(ind1,ind2)=amatch(k,kk)
        amat(indbeg+ind1-ind2,ind2)=amatch(k,kk)
      ENDDO
    ENDDO
!
! End spline continuity conditions
!
  ENDDO
!
  CALL locspline_mf(npowers+1,nfun,eta(1:npowers+1),f(1:npowers+1,:), &
                 eta(1),ails)
  DO i=1,nshift
    ind1=i
    bvec(ind1,:)=ails(i+1,:)
    ind2=i+1
!    amat(ind1,ind2)=1.d0
    amat(indbeg+ind1-ind2,ind2)=1.d0
  ENDDO
!
  CALL locspline_mf(npowers+1,nfun,eta(nbands-npowers:nbands), &
                    f(nbands-npowers:nbands,:),eta(nbands),ails)
!
  ind1=(npowers+1)*(nbands-1)+nshift+1
  bvec(ind1,:)=f(nbands,:)
  ind2=(npowers+1)*(nbands-1)+1
!  amat(ind1,ind2)=1.d0
  amat(indbeg+ind1-ind2,ind2)=1.d0
  DO i=1,nshift
    ind1=(npowers+1)*(nbands-1)+nshift+1+i
    bvec(ind1,:)=ails(i+1,:)
    ind2=(npowers+1)*(nbands-1)+1+i
!    amat(ind1,ind2)=1.d0
    amat(indbeg+ind1-ind2,ind2)=1.d0
  ENDDO
!
  CALL dgbsv(ndim,kl,ku,nfun,amat,ldab,ipiv,bvec,ndim,info)
!
  IF(info.NE.0) THEN
    PRINT *,'oddorderspline: dgbsv error',info
    ierr=2
    DEALLOCATE(amat,bvec,ipiv,powdeleta,amatch,coefbin)
    RETURN
  ENDIF
!
  DO i=1,nbands-1
    ai(0:npowers,i,:)=bvec((npowers+1)*(i-1)+1:(npowers+1)*i,:)
  ENDDO
!  ai(0,nbands,:)=f(nbands,:)
  ai(:,nbands,:)=ails(:,:)
!
  DEALLOCATE(amat,bvec,ipiv,powdeleta,amatch,coefbin)
  ierr=0
!
  END SUBROUTINE oddorderspline
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE locspline_mf(npoi,nfun,xp,fp,x0,ai)
!
! Given "npoi" argument values, "xp", and "npoi" function values "fp"
! computes interpolation polynomial coefficients, "ai" at the point "x0".
! Order of polynomial is "npoi-1". Number of interpolated functions is "nfun"
!
  IMPLICIT NONE
!
  INTEGER :: npoi,nfun,info,i,ndim,j,ldab,kl,ku,indbeg
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipiv
  DOUBLE PRECISION                              :: x0
  DOUBLE PRECISION, DIMENSION(npoi)             :: xp
  DOUBLE PRECISION, DIMENSION(npoi,nfun)        :: fp,ai
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,ab,bvec
!
  ndim=npoi
  kl=ndim-1
  ku=ndim-1
  ldab=2*kl+ku+1
  indbeg=kl+ku+1
  ALLOCATE(amat(ndim,ndim),ab(ldab,ndim),bvec(ndim,nfun),ipiv(ndim))
  amat(:,1)=1.d0
  amat(:,2)=xp-x0
  DO i=3,npoi
    amat(:,i)=amat(:,i-1)*amat(:,2)
  ENDDO
!
  do i=1,ndim
    do j=1,ndim
      ab(indbeg+i-j,j)=amat(i,j)
    enddo
  enddo
!
  bvec=fp
!
  CALL dgbsv(ndim,kl,ku,nfun,ab,ldab,ipiv,bvec,ndim,info)
!
  ai=bvec
!
  DEALLOCATE(amat,ab,bvec,ipiv)
!
  END SUBROUTINE locspline_mf
!
END MODULE oddorderspline_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine invert_mono_reg(nx,arry,xmin,xmax,ny,arrx,ymin,ymax)
!
! Inverts the monotonous function y(x) given on the equidistant grid
! of x values on the interval [xmin,xmax] by the array y_i=arry(i).
! The result, function x(y), is given on the equidistant grid of y values
! at the interval [ymin,ymax] by the array x_i=arrx(i).
!
  USE oddorderspline_mod
!
  implicit none
!
  integer, parameter :: nhalfband=2, npowers=2*nhalfband+1
!
  integer :: ny,nx,iy,ix,npoints,k,m,ierr
!
  double precision :: xmin,xmax,ymin,ymax,hy,hx,yeq,dely
  double precision, dimension(0:nx) :: arry
  double precision, dimension(0:ny) :: arrx
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: y
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: f
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: ai
!
  npoints=nx+1
  allocate(y(1:npoints),f(1:npoints,1))
!
  ymin=arry(0)
  ymax=arry(nx)
!
  hy=(ymax-ymin)/ny
  hx=(xmax-xmin)/nx
!
  if(ny.ge.nx) then
    y=arry
    do ix=0,nx
      f(ix+1,1)=xmin+hx*dfloat(ix)
    enddo
  else
    npoints=1
    y(1)=ymin
    f(1,1)=xmin
    do ix=1,nx
      if((arry(ix)-ymin)/hy.ge.dfloat(npoints).or.ix.eq.nx) then
        npoints=npoints+1
        y(npoints)=arry(ix)
        f(npoints,1)=xmin+hx*dfloat(ix)
      endif
    enddo
  endif
!
  CALL oddorderspline(nhalfband,npoints,1,y(1:npoints),f(1:npoints,1),ai,ierr)
!
  arrx(0)=xmin
  arrx(ny)=xmax
!
  k=0
  do iy=1,ny-1
    yeq=ymin+hy*dfloat(iy)
    do while(y(k+1).lt.yeq)
      k=k+1
    enddo
    dely=yeq-y(k)
    arrx(iy)=0.d0
    DO m=npowers,0,-1
      arrx(iy)=arrx(iy)*dely+ai(m,k,1)
    ENDDO
  enddo
!
  deallocate(y,f,ai)
!
  return
  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine invert_mono_per(nx,arry_in,xmin,xmax,ny,arrx,ymin,ymax)
!
! Inverts the monotonous function y(x) given on the equidistant grid
! of x values on the interval [xmin,xmax] by the array y_i=arry(i).
! Special case: y(x) is a sum of linear and periodic functions.
! The result, function x(y), is given on the equdistant grid of y values
! at the interval [ymin,ymax] by the array x_i=arrx(i).
!
  implicit none
!
  integer :: ny,nx,iy,ix,ixfix,ix1,ix2,ix3,ix4
!
  double precision :: xmin,xmax,ymin,ymax,hy,y,hx,x1,x2,x3,x4,y1,y2,y3,y4
  double precision, dimension(0:nx) :: arry_in
  double precision, dimension(0:ny) :: arrx
  double precision, dimension(:), allocatable :: arry
!
  allocate(arry(-1:nx+1))
  arry(0:nx)=arry_in
!
  ymin=arry(0)
  ymax=arry(nx)
  arry(-1)=arry(nx-1)-ymax+ymin
  arry(nx+1)=arry(1)+ymax-ymin
!
  hy=(ymax-ymin)/ny
  hx=(xmax-xmin)/nx
!
  arrx(0)=xmin
  arrx(ny)=xmax
!
  ixfix=1
  do iy=1,ny-1
    y=ymin+iy*hy
    do ix=0,nx
      if(arry(ix).gt.y) then
        ixfix=ix-3
        exit
      endif
    enddo
!    ixfix=max(ixfix,-1)
!    ixfix=min(ixfix,nx-4)
    ix1=ixfix+1
    ix2=ixfix+2
    ix3=ixfix+3
    ix4=ixfix+4
    x1=xmin+ix1*hx
    x2=xmin+ix2*hx
    x3=xmin+ix3*hx
    x4=xmin+ix4*hx
    y1=arry(ix1)
    y2=arry(ix2)
    y3=arry(ix3)
    y4=arry(ix4)
    arrx(iy) = x1*(y-y2)/(y1-y2)*(y-y3)/(y1-y3)*(y-y4)/(y1-y4)    &
             + x2*(y-y3)/(y2-y3)*(y-y4)/(y2-y4)*(y-y1)/(y2-y1)    &
             + x3*(y-y4)/(y3-y4)*(y-y1)/(y3-y1)*(y-y2)/(y3-y2)    &
             + x4*(y-y1)/(y4-y1)*(y-y2)/(y4-y2)*(y-y3)/(y4-y3)
  enddo
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_theta
!
  use theta_rz_mod
  use input_files, only : iunit,fluxdatapath
  use oddorderspline_mod
!
  implicit none
!
  integer, parameter :: nhalfband=2, ncoefs=2*nhalfband+2, nfun=1
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: nsqpsi,nlabel,ntheta,i,ierr,k
  real(kind=8) :: sqpsimin,sqpsimax
  real(kind=8) :: flabel_min,flabel_max,xarrfine
  double precision, dimension(:,:),   allocatable :: fa
  double precision, dimension(:,:,:), allocatable :: ai
  real(kind=8), dimension(:,:), allocatable :: theta_of_theta_qt
!
  open(iunit,form='unformatted',                                 &
       file=trim(fluxdatapath)//'/theta_of_theta_qt_flabel_new.dat')
  read (iunit) nsqpsi,nlabel,ntheta,sqpsimin,sqpsimax,flabel_min,flabel_max &
              ,raxis,zaxis,psiaxis,sigma_qt
  allocate(theta_of_theta_qt(nlabel,0:ntheta),sqpsi_lab(0:nlabel))
  read (iunit) theta_of_theta_qt
  read (iunit) sqpsi_lab
  close(iunit)
!
  nsqp=nsqpsi
  nlab=nlabel
  nthe=ntheta+1
!
  hsqpsi=(sqpsimax-sqpsimin)/(nsqp-1)
  hlabel=(flabel_max-flabel_min)/(nlab-1)
  htheqt=2.d0*pi/ntheta
!
  allocate(sqpsi(nsqp),fa(0:nlab,nfun),theqt(nthe),flab(nlab))
!
  do i=1,nsqp
    sqpsi(i)=sqpsimin+hsqpsi*(i-1)
  enddo
!
  do i=0,nlab
    fa(i,1)=hlabel*dfloat(i)
  enddo
!
  flab=fa(1:nlab,1)
!
  do i=1,nthe
    theqt(i)=htheqt*(i-1)
  enddo
!
  icp_pt=nthe*nlab
  allocate( splthet(6,6,icp_pt), ipoint_pt(nlab,nthe) )
!
  call s2dring(nlab,nthe,hlabel,htheqt,theta_of_theta_qt(:,0:ntheta), &
               icp_pt,splthet,ipoint_pt)
!
  allocate(spllabel(ncoefs,0:nlab))
!
  call oddorderspline(nhalfband,nlab+1,nfun,sqpsi_lab,fa,ai,ierr)
!
  spllabel=ai(:,:,1)
  deallocate(fa,ai)
  h_subgrid_fl=minval(sqpsi_lab(1:nlab)-sqpsi_lab(0:nlab-1))
  n_subgrid_fl=ceiling((sqpsi_lab(nlab)-sqpsi_lab(0))/h_subgrid_fl)
  allocate(ind_subgrid_fl(0:n_subgrid_fl))
!
  k=0
  do i=0,n_subgrid_fl
    xarrfine=sqpsi_lab(0)+h_subgrid_fl*dfloat(i-1)
    if(xarrfine.ge.sqpsi_lab(k+1)) k=k+1
    k=min(k,nlab)
    ind_subgrid_fl(i)=k
  enddo
!
  return
  end
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine psithet_rz(rrr,zzz,                                          &
                        theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                        flabel,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)
!
  use theta_rz_mod
  use field_eq_mod, only : nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint  &
                         , psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use extract_fluxcoord_mod, only : psif_extract,theta_extract
!
  implicit none
!
  real(kind=8), parameter :: pi=3.14159265358979d0
!
  integer :: i,ierr,k
  real(kind=8) :: rrr,zzz,theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz
  real(kind=8) :: theta_s,theta_t,theta_ss,theta_st,theta_tt
  real(kind=8) :: sqpsi_qt,s_r,s_z,s_rr,s_rz,s_zz
  real(kind=8) :: theta_qt,t_r,t_z,t_rr,t_rz,t_zz
  real(kind=8) :: rho2,rho4,dr,dz,flabel,dflabel,ddflabel,dx,dfl_dpsi,ddfl_dpsi
  real(kind=8) :: s0,ds0ds,dds0ds
!
  if(icall.eq.0) then
    icall=1
    call load_theta
  endif
!
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
              psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
!
  sqpsi_qt=sqrt(abs(psif-psiaxis))
!
  i=int(sqpsi_qt/h_subgrid_fl)
  i=min(n_subgrid_fl,i)
  k=ind_subgrid_fl(i)
  dx=sqpsi_qt-sqpsi_lab(k)
  flabel=spllabel(1,k)+dx*(spllabel(2,k)+dx*(spllabel(3,k)           &
        +dx*(spllabel(4,k)+dx*(spllabel(5,k)+dx*spllabel(6,k)))))
  dflabel=spllabel(2,k)+dx*(2.d0*spllabel(3,k)                       &
          +dx*(3.d0*spllabel(4,k)+dx*(4.d0*spllabel(5,k)             &
          +dx*5.d0*spllabel(6,k))))
  ddflabel=2.d0*spllabel(3,k)+dx*(6.d0*spllabel(4,k)                 &
           +dx*(12.d0*spllabel(5,k)+dx*20.d0*spllabel(6,k)))
!
  if(abs(sqpsi_qt).gt.100.d0*hsqpsi*epsilon(1.d0)) then
    dfl_dpsi=0.5d0*dflabel/sqpsi_qt
    ddfl_dpsi=0.25d0*(ddflabel-dflabel/sqpsi_qt)/abs(psif-psiaxis)
  else
    dfl_dpsi=0.5d0*dflabel/(100.d0*hsqpsi*epsilon(1.d0))
    ddfl_dpsi=0.25d0*(ddflabel-dflabel/(100.d0*hsqpsi*epsilon(1.d0)))&
             /(100.d0*hsqpsi*epsilon(1.d0))
  endif
  s_r=dpsidr*dfl_dpsi
  s_z=dpsidz*dfl_dpsi
  s_rr=d2psidr2*dfl_dpsi+dpsidr**2*ddfl_dpsi
  s_rz=d2psidrdz*dfl_dpsi+dpsidr*dpsidz*ddfl_dpsi
  s_zz=d2psidz2*dfl_dpsi+dpsidz**2*ddfl_dpsi
!
  s0=sqpsi_qt
  ds0ds=1.d0/dflabel
  dds0ds=-ds0ds**3*ddflabel
!
  dr=rrr-raxis
  if(abs(dr).lt.raxis*epsilon(1.d0)) dr=raxis*epsilon(1.d0)
  dz=zzz-zaxis
  rho2=dr**2+dz**2
  rho4=rho2**2
  theta_qt=mod(sigma_qt*atan2(dz,dr)+2*pi,2*pi)
  t_r=-sigma_qt*dz/rho2
  t_z=sigma_qt*dr/rho2
  t_rr=2.d0*sigma_qt*dr*dz/rho4
  t_zz=-t_rr
  t_rz=sigma_qt*(dz**2-dr**2)/rho4
!
  call spline(nlab,nthe,flab,theqt,hlabel,htheqt,icp_pt,splthet,ipoint_pt, &
              flabel,theta_qt,                                             &
              theta,theta_s,theta_t,theta_ss,theta_st,theta_tt,ierr)
!
  theta=theta+theta_qt
  theta_r=theta_s*s_r+(theta_t+1.d0)*t_r
  theta_z=theta_s*s_z+(theta_t+1.d0)*t_z
  theta_rr=theta_ss*s_r**2+2.d0*theta_st*s_r*t_r+theta_tt*t_r**2 &
          +theta_s*s_rr+(theta_t+1.d0)*t_rr
  theta_rz=theta_ss*s_r*s_z+theta_st*(s_r*t_z+s_z*t_r)+theta_tt*t_r*t_z &
          +theta_s*s_rz+(theta_t+1.d0)*t_rz
  theta_zz=theta_ss*s_z**2+2.d0*theta_st*s_z*t_z+theta_tt*t_z**2 &
          +theta_s*s_zz+(theta_t+1.d0)*t_zz
!
  psif_extract=psif
  theta_extract=theta
!
  return
  end subroutine psithet_rz
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine cspl_five_reg(n,h,a,b,c,d,e,f)
!
  implicit none
!
  integer :: n,i,ip1
  double precision :: h,rhop,rhom,fac
  double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,det
  double complex :: abeg,bbeg,cbeg,dbeg,ebeg,fbeg
  double complex :: aend,bend,cend,dend,eend,fend
  double complex :: b1,b2,b3
  double complex, dimension(n) :: a,b,c,d,e,f
  double complex, dimension(:), allocatable :: alp,bet,gam
!
  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)
!
  a11=1.d0
  a12=1.d0/4.d0
  a13=1.d0/16.d0
  a21=3.d0
  a22=27.d0/4.d0
  a23=9.d0*27.d0/16.d0
  a31=5.d0
  a32=125.d0/4.d0
  a33=5.d0**5/16.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)-a(3)
  b2=a(5)-a(2)
  b3=a(6)-a(1)
  bbeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bbeg=bbeg/det
  dbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dbeg=dbeg/det
  fbeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fbeg=fbeg/det
  b1=a(n-2)-a(n-3)
  b2=a(n-1)-a(n-4)
  b3=a(n)-a(n-5)
  bend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  bend=bend/det
  dend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  dend=dend/det
  fend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  fend=fend/det
  a11=2.d0
  a12=1.d0/2.d0
  a13=1.d0/8.d0
  a21=2.d0
  a22=9.d0/2.d0
  a23=81.d0/8.d0
  a31=2.d0
  a32=25.d0/2.d0
  a33=625.d0/8.d0
  det=a11*a22*a33+a12*a23*a31+a13*a21*a32-a12*a21*a33-a13*a22*a31-a11*a23*a32
  b1=a(4)+a(3)
  b2=a(5)+a(2)
  b3=a(6)+a(1)
  abeg=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  abeg=abeg/det
  cbeg=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cbeg=cbeg/det
  ebeg=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  ebeg=ebeg/det
  b1=a(n-2)+a(n-3)
  b2=a(n-1)+a(n-4)
  b3=a(n)+a(n-5)
  aend=b1*a22*a33+a12*a23*b3+a13*b2*a32-a12*b2*a33-a13*a22*b3-b1*a23*a32
  aend=aend/det
  cend=a11*b2*a33+b1*a23*a31+a13*a21*b3-b1*a21*a33-a13*b2*a31-a11*a23*b3
  cend=cend/det
  eend=a11*a22*b3+a12*b2*a31+b1*a21*a32-a12*a21*b3-b1*a22*a31-a11*b2*a32
  eend=eend/det
!
  allocate(alp(n),bet(n),gam(n))
!
  alp(1)=0.0d0
  bet(1)=ebeg*(2.d0+rhom)-5.d0*fbeg*(3.d0+1.5d0*rhom) !gamma1
!
  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  enddo
!
  gam(n-2)=eend*(2.d0+rhom)+5.d0*fend*(3.d0+1.5d0*rhom) !gamma
  do i=n-3,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  enddo
!
  alp(1)=0.0d0
  bet(1)=ebeg-2.5d0*5.d0*fbeg !e1
!
  do i=1,n-2
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  enddo
!
  e(n)=eend+2.5d0*5.d0*fend
  e(n-1)=e(n)*alp(n-1)+bet(n-1)
  f(n-1)=(e(n)-e(n-1))/5.d0
  e(n-2)=e(n-1)*alp(n-2)+bet(n-2)
  f(n-2)=(e(n-1)-e(n-2))/5.d0
  d(n-2)=dend+1.5d0*4.d0*eend+1.5d0**2*10.d0*fend
!
  do i=n-3,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
    f(i)=(e(i+1)-e(i))/5.d0
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  enddo
!
  do i=n-3,n
    b(i)=b(i-1)+2.d0*c(i-1)+3.d0*d(i-1)+4.d0*e(i-1)+5.d0*f(i-1)
    c(i)=c(i-1)+3.d0*d(i-1)+6.d0*e(i-1)+10.d0*f(i-1)
    d(i)=d(i-1)+4.d0*e(i-1)+10.d0*f(i-1)
    if(i.ne.n) f(i)= a(i+1)-a(i)-b(i)-c(i)-d(i)-e(i)
  enddo
  f(n)=f(n-1)
!
  fac=1.d0/h
  b=b*fac
  fac=fac/h
  c=c*fac
  fac=fac/h
  d=d*fac
  fac=fac/h
  e=e*fac
  fac=fac/h
  f=f*fac
!
  deallocate(alp,bet,gam)
!
  return
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_fourier(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ              &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Caution: derivatives are not computed, for derivatives call
! a driver routine "field_fourier_derivs"
!
  use new_amn_mod
  use input_files,           only : iunit
  use inthecore_mod, only : incore,plaf
  use oddorderspline_mod
  use binsrc_mod
!
  implicit none
!
  integer :: m,n,i,k,ierr,nfun,npoi_amn,ibeg,iend
  double precision :: rrinv
  double precision :: r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ                &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  double precision :: flabel,dx,g11,g12,g11_r,g11_z,g12_r,g12_z
  double precision :: theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                      s_r,s_z,s_rr,s_rz,s_zz
  double precision :: apsi,apsi_s,apsi_t,apsi_p
  double precision :: athe,athe_s,athe_t,athe_p
  double precision :: aphi,aphi_s,aphi_t,aphi_p
  double precision :: delbr,delbz,delbp,delar,delaz
  double precision :: fcjac,g11_t,g12_t,s0,ds0ds,dds0ds
  double precision :: sqg_gpp,sqg_gpp_s,sqg_gpp_ss
  double precision :: gpp_sqg,gpp_sqg_s,gpp_sqg_ss
!
  double complex :: expon,dexpon,ddexpon,cdummy
  double complex :: cas,cas_s,cas_t,cas_p
  double complex :: cas_ss,cas_tt,cas_pp,cas_st,cas_sp,cas_tp
  double complex :: cat,cat_s,cat_t,cat_p
  double complex :: cat_ss,cat_tt,cat_pp,cat_st,cat_sp,cat_tp
  double complex :: casrr,casrr_s,casrr_t,casrr_p
  double complex :: casrr_ss,casrr_tt,casrr_pp,casrr_st,casrr_sp,casrr_tp
  double complex :: cap,cap_s,cap_t,cap_p
  double complex :: cap_ss,cap_tt,cap_pp,cap_st,cap_sp,cap_tp
  double precision, dimension(:),     allocatable :: sqg_ov_gphph
  double precision, dimension(:,:),   allocatable :: fa
  double precision, dimension(:,:,:), allocatable :: ai
!
  double complex,   dimension(:,:),   allocatable :: athetm_plas,alpham_plas
  double complex,   dimension(:,:),   allocatable :: etam_plas,aphim_plas
!
! Initialization ------------------------------------------------------------
!
  if(prop) then
    prop=.false.
!
! Toroidally symetric part of the vacuum perturbation field - comes now
! from the cylindrical vacuum field routine
!
    isw_gauge=1
!
! Fourier ampitudes of the non-axisymmetric vacuum perturbation field:
!
    if(isw_gauge.eq.1) then
!
! Gauge $A_\varphi=0$:
!
!      open(iunit,form='unformatted',                                         &
!           file=trim(fluxdatapath)//'/amn_plas_Aphi_zero.dat')
      open(iunit,form='unformatted',file='amn_plas_Aphi_zero.dat')
      read (iunit) ntor_act,mpol,nlabel_tot
      allocate(ntor_arr(ntor_act),ibeg_arr(ntor_act),iend_arr(ntor_act))
      allocate(flabel_arr(nlabel_tot),sqg_ov_gphph(nlabel_tot))
      allocate(athetm_plas(nlabel_tot,-mpol:mpol))
      allocate(alpham_plas(nlabel_tot,-mpol:mpol))
      read (iunit) ntor_arr,ibeg_arr,iend_arr,flabel_arr,sqg_ov_gphph
      read (iunit) athetm_plas,alpham_plas
      close(iunit)
!
      allocate(splatet(-mpol:mpol,ncoefs,1:nlabel_tot))
      allocate(splaalp(-mpol:mpol,ncoefs,1:nlabel_tot))
      nfun=2*mpol+1
!
      do n=1,ntor_act
        npoi_amn=iend_arr(n)-ibeg_arr(n)+1
        allocate(fa(1:npoi_amn,nfun))
        ibeg=ibeg_arr(n)
        iend=iend_arr(n)
!
! spline $A_\vartheta$:
        fa=real(athetm_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splatet(:,i,ibeg:iend)=transpose(dcmplx(ai(i-1,:,:),0.d0))
        enddo
        fa=dimag(athetm_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splatet(:,i,ibeg:iend)=splatet(:,i,ibeg:iend)                        &
                                +transpose(dcmplx(0.d0,ai(i-1,:,:)))
        enddo
!
! spline $\alpha$:
        fa=real(alpham_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaalp(:,i,ibeg:iend)=transpose(dcmplx(ai(i-1,:,:),0.d0))
        enddo
        fa=dimag(alpham_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaalp(:,i,ibeg:iend)=splaalp(:,i,ibeg:iend)                        &
                                +transpose(dcmplx(0.d0,ai(i-1,:,:)))
        enddo
!
        deallocate(fa)
      enddo
!
      deallocate(ai,athetm_plas,alpham_plas)
!
      allocate(amalp(-mpol:mpol),amtet(-mpol:mpol))
      allocate(amalp_s(-mpol:mpol),amtet_s(-mpol:mpol))
      allocate(amalp_ss(-mpol:mpol),amtet_ss(-mpol:mpol))
!
    else
!
! Coulomb gauge:
!
!      open(iunit,form='unformatted',                                         &
!           file=trim(fluxdatapath)//'/amn_plas_Coul.dat')
      open(iunit,form='unformatted',file='amn_plas_Coul.dat')
      read (iunit) ntor_act,mpol,nlabel_tot
      allocate(ntor_arr(ntor_act),ibeg_arr(ntor_act),iend_arr(ntor_act))
      allocate(flabel_arr(nlabel_tot),sqg_ov_gphph(nlabel_tot))
      allocate(etam_plas(nlabel_tot,-mpol:mpol))
      allocate(athetm_plas(nlabel_tot,-mpol:mpol))
      allocate(aphim_plas(nlabel_tot,-mpol:mpol))
      read (iunit) ntor_arr,ibeg_arr,iend_arr,flabel_arr,sqg_ov_gphph
      read (iunit) etam_plas,athetm_plas,aphim_plas
      close(iunit)
!
      allocate(splaeta(-mpol:mpol,ncoefs,1:nlabel_tot))
      allocate(splatet(-mpol:mpol,ncoefs,1:nlabel_tot))
      allocate(splaphi(-mpol:mpol,ncoefs,1:nlabel_tot))
      nfun=2*mpol+1
!
      do n=1,ntor_act
        npoi_amn=iend_arr(n)-ibeg_arr(n)+1
        allocate(fa(1:npoi_amn,nfun))
        ibeg=ibeg_arr(n)
        iend=iend_arr(n)
!
! spline $\eta$:
        fa=real(etam_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaeta(:,i,ibeg:iend)=transpose(dcmplx(ai(i-1,:,:),0.d0))
        enddo
        fa=dimag(etam_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaeta(:,i,ibeg:iend)=splaeta(:,i,ibeg:iend)                        &
                                +transpose(dcmplx(0.d0,ai(i-1,:,:)))
        enddo
!
! spline $A_\vartheta$:
        fa=real(athetm_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splatet(:,i,ibeg:iend)=transpose(dcmplx(ai(i-1,:,:),0.d0))
        enddo
        fa=dimag(athetm_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splatet(:,i,ibeg:iend)=splatet(:,i,ibeg:iend)                        &
                                +transpose(dcmplx(0.d0,ai(i-1,:,:)))
        enddo
!
! spline $A_\varphi$:
        fa=real(aphim_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaphi(:,i,ibeg:iend)=transpose(dcmplx(ai(i-1,:,:),0.d0))
        enddo
        fa=dimag(aphim_plas(ibeg:iend,:))
!
        call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                            ai,ierr)
!
        do i=1,ncoefs
          splaphi(:,i,ibeg:iend)=splaphi(:,i,ibeg:iend)                        &
                                +transpose(dcmplx(0.d0,ai(i-1,:,:)))
        enddo
!
        deallocate(fa)
      enddo
!
      deallocate(ai,athetm_plas,etam_plas,aphim_plas)
!
      allocate(ameta(-mpol:mpol),amtet(-mpol:mpol),amphi(-mpol:mpol))
      allocate(ameta_s(-mpol:mpol),amtet_s(-mpol:mpol),amphi_s(-mpol:mpol))
      allocate(ameta_ss(-mpol:mpol),amtet_ss(-mpol:mpol),amphi_ss(-mpol:mpol))
!
    endif
!
    nfun=1
    allocate(splsqggpp(ncoefs,1:nlabel_tot))
!
    do n=1,ntor_act
      npoi_amn=iend_arr(n)-ibeg_arr(n)+1
      allocate(fa(1:npoi_amn,nfun))
      ibeg=ibeg_arr(n)
      iend=iend_arr(n)
      fa(:,1)=sqg_ov_gphph(ibeg:iend)
!
      call oddorderspline(nhalfband,npoi_amn,nfun,flabel_arr(ibeg:iend),fa,  &
                          ai,ierr)
!
      do i=1,ncoefs
        splsqggpp(i,ibeg:iend)=ai(i-1,:,1)
      enddo
!
      deallocate(fa)
    enddo
!
    deallocate(ai,sqg_ov_gphph)
!
    allocate(expthe(-mpol:mpol),expphi(ntor_arr(ntor_act)))
    allocate(dthe(-mpol:mpol),ddthe(-mpol:mpol))
!
    do m=-mpol,mpol
      dthe(m)=dcmplx(0.d0,dfloat(m))
      ddthe(m)=dcmplx(-dfloat(m)**2,0.d0)
    enddo
!
  endif
  sqg_gpp=0.d0
!
! End of initialization ------------------------------------------------------
!
! Toroidally symmetric part of the perturbation field - not computed, comes
! from the vacuum routine
!
  Br=0.d0
  Bp=0.d0
  Bz=0.d0
  dBrdR=0.d0
  dBrdZ=0.d0
  dBrdp=0.d0
  dBpdR=0.d0
  dBpdZ=0.d0
  dBpdp=0.d0
  dBzdR=0.d0
  dBzdZ=0.d0
  dBzdp=0.d0
!
! Asymmetric part of the perturbation field:
!
  call psithet_rz(r,z,                                              &
                  theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                  flabel,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)
!
  g11=s_r**2+s_z**2
  g12=s_r*theta_r+s_z*theta_z
  g11_r=2.d0*s_r*s_rr+2.d0*s_z*s_rz
  g11_z=2.d0*s_r*s_rz+2.d0*s_z*s_zz
  g12_r=s_rr*theta_r+s_r*theta_rr+s_rz*theta_z+s_z*theta_rz
  g12_z=s_rz*theta_r+s_r*theta_rz+s_zz*theta_z+s_z*theta_zz
  fcjac=s_r*theta_z-s_z*theta_r
  g11_t=(s_r*g11_z-s_z*g11_r)/fcjac
  g12_t=(s_r*g12_z-s_z*g12_r)/fcjac
!
  expthe(0)=(1.d0,0.d0)
  expthe(1)=exp(dcmplx(0.d0,theta))
  expthe(-1)=conjg(expthe(1))
  do m=2,mpol
    expthe(m)=expthe(m-1)*expthe(1)
    expthe(-m)=conjg(expthe(m))
  enddo
  expphi(1)=exp(dcmplx(0.d0,phi))
  do n=2,ntor_arr(ntor_act)
    expphi(n)=expphi(n-1)*expphi(1)
  enddo
!
  if(isw_gauge.eq.1) then
!
! Gauge $A_\varphi=0$:
!
    cas=(0.d0,0.d0)
    cas_s=(0.d0,0.d0)
    cas_t=(0.d0,0.d0)
    cas_p=(0.d0,0.d0)
    cas_ss=(0.d0,0.d0)
    cas_tt=(0.d0,0.d0)
    cas_pp=(0.d0,0.d0)
    cas_sp=(0.d0,0.d0)
    cas_st=(0.d0,0.d0)
    cas_tp=(0.d0,0.d0)
    cat=(0.d0,0.d0)
    cat_s=(0.d0,0.d0)
    cat_t=(0.d0,0.d0)
    cat_p=(0.d0,0.d0)
    cat_ss=(0.d0,0.d0)
    cat_tt=(0.d0,0.d0)
    cat_pp=(0.d0,0.d0)
    cat_sp=(0.d0,0.d0)
    cat_st=(0.d0,0.d0)
    cat_tp=(0.d0,0.d0)
!
    do n=1,ntor_act
      ibeg=ibeg_arr(n)
      iend=iend_arr(n)
      expon=expphi(ntor_arr(n))
      dexpon=expon*dcmplx(0.d0,ntor_arr(n))
      ddexpon=dexpon*dcmplx(0.d0,ntor_arr(n))
!
      call binsrc(flabel_arr(ibeg:iend),ibeg,iend,flabel,k)
!
      k=max(ibeg,k-1)
      dx=flabel-flabel_arr(k)
      amalp=(0.d0,0.d0)
      amtet=(0.d0,0.d0)
      sqg_gpp=0.d0
      do i=ncoefs,1,-1
        amalp=amalp*dx+splaalp(:,i,k)
        amtet=amtet*dx+splatet(:,i,k)
        sqg_gpp=sqg_gpp*dx+splsqggpp(i,k)
      enddo
      cdummy=sum(amalp*expthe)
      cas=cas+cdummy*expon
      cas_p=cas_p+cdummy*dexpon
      cas_pp=cas_pp+cdummy*ddexpon
      cdummy=sum(amalp*expthe*dthe)
      cas_t=cas_t+cdummy*expon
      cas_tp=cas_tp+cdummy*dexpon
      cdummy=sum(amalp*expthe*ddthe)
      cas_tt=cas_tt+cdummy*expon
      cdummy=sum(amtet*expthe)
      cat=cat+cdummy*expon
      cat_p=cat_p+cdummy*dexpon
      cat_pp=cat_pp+cdummy*ddexpon
      cdummy=sum(amtet*expthe*dthe)
      cat_t=cat_t+cdummy*expon
      cat_tp=cat_tp+cdummy*dexpon
      cdummy=sum(amtet*expthe*ddthe)
      cat_tt=cat_tt+cdummy*expon
!
      amalp_s=(0.d0,0.d0)
      amtet_s=(0.d0,0.d0)
      sqg_gpp_s=0.d0
      do i=ncoefs,2,-1
        amalp_s=amalp_s*dx+splaalp(:,i,k)*dfloat(i-1)
        amtet_s=amtet_s*dx+splatet(:,i,k)*dfloat(i-1)
        sqg_gpp_s=sqg_gpp_s*dx+splsqggpp(i,k)*dfloat(i-1)
      enddo
      cdummy=sum(amalp_s*expthe)
      cas_s=cas_s+cdummy*expon
      cas_sp=cas_sp+cdummy*dexpon
      cdummy=sum(amalp_s*expthe*dthe)
      cas_st=cas_st+cdummy*expon
      cdummy=sum(amtet_s*expthe)
      cat_s=cat_s+cdummy*expon
      cat_sp=cat_sp+cdummy*dexpon
      cdummy=sum(amtet_s*expthe*dthe)
      cat_st=cat_st+cdummy*expon
!
      amalp_ss=(0.d0,0.d0)
      amtet_ss=(0.d0,0.d0)
      sqg_gpp_ss=0.d0
      do i=ncoefs,3,-1
        amalp_ss=amalp_ss*dx+splaalp(:,i,k)*dfloat(i-1)*dfloat(i-2)
        amtet_ss=amtet_ss*dx+splatet(:,i,k)*dfloat(i-1)*dfloat(i-2)
        sqg_gpp_ss=sqg_gpp_ss*dx+splsqggpp(i,k)*dfloat(i-1)*dfloat(i-2)
      enddo
      cdummy=sum(amalp_ss*expthe)
      cas_ss=cas_ss+cdummy*expon
      cdummy=sum(amtet_ss*expthe)
      cat_ss=cat_ss+cdummy*expon
    enddo
!
    gpp_sqg=1.d0/sqg_gpp
    gpp_sqg_s=-sqg_gpp_s*gpp_sqg**2
    gpp_sqg_ss=2.d0*gpp_sqg_s**2*sqg_gpp-sqg_gpp_ss*gpp_sqg**2
!
    cas_ss=cas_ss*gpp_sqg+2.d0*cas_s*gpp_sqg_s+cas*gpp_sqg_ss
    cas_st=cas_st*gpp_sqg+cas_t*gpp_sqg_s
    cas_sp=cas_sp*gpp_sqg+cas_p*gpp_sqg_s
    cas_tt=cas_tt*gpp_sqg
    cas_tp=cas_tp*gpp_sqg
    cas_pp=cas_pp*gpp_sqg
    cas_s=cas_s*gpp_sqg+cas*gpp_sqg_s
    cas_t=cas_t*gpp_sqg
    cas_p=cas_p*gpp_sqg
    cas=cas*gpp_sqg
!
    apsi=2.d0*real(cas)
    apsi_s=2.d0*real(cas_s)
    apsi_t=2.d0*real(cas_t)
    apsi_p=2.d0*real(cas_p)
    athe=2.d0*real(cat)
    athe_s=2.d0*real(cat_s)
    athe_t=2.d0*real(cat_t)
    athe_p=2.d0*real(cat_p)
!
    delar=(apsi-g12*athe)/g11*s_r+athe*theta_r
    delaz=(apsi-g12*athe)/g11*s_z+athe*theta_z
    delbr=((apsi_p-g12*athe_p)/g11*s_z+athe_p*theta_z)/r
    delbz=-((apsi_p-g12*athe_p)/g11*s_r+athe_p*theta_r)/r
    delbp=fcjac*( (apsi_t-g12*athe_t-g12_t*athe)/g11                           &
         -(apsi-g12*athe)*g11_t/g11**2 )                                       &
         +athe_s*(theta_r*s_z-theta_z*s_r)
  else
!
! Coulomb gauge:
!
    casrr=(0.d0,0.d0)
    casrr_s=(0.d0,0.d0)
    casrr_t=(0.d0,0.d0)
    casrr_p=(0.d0,0.d0)
    casrr_ss=(0.d0,0.d0)
    casrr_tt=(0.d0,0.d0)
    casrr_pp=(0.d0,0.d0)
    casrr_sp=(0.d0,0.d0)
    casrr_st=(0.d0,0.d0)
    casrr_tp=(0.d0,0.d0)
    cat=(0.d0,0.d0)
    cat_s=(0.d0,0.d0)
    cat_t=(0.d0,0.d0)
    cat_p=(0.d0,0.d0)
    cat_ss=(0.d0,0.d0)
    cat_tt=(0.d0,0.d0)
    cat_pp=(0.d0,0.d0)
    cat_sp=(0.d0,0.d0)
    cat_st=(0.d0,0.d0)
    cat_tp=(0.d0,0.d0)
    cap=(0.d0,0.d0)
    cap_s=(0.d0,0.d0)
    cap_t=(0.d0,0.d0)
    cap_p=(0.d0,0.d0)
    cap_ss=(0.d0,0.d0)
    cap_tt=(0.d0,0.d0)
    cap_pp=(0.d0,0.d0)
    cap_sp=(0.d0,0.d0)
    cap_st=(0.d0,0.d0)
    cap_tp=(0.d0,0.d0)
!
    do n=1,ntor_act
      ibeg=ibeg_arr(n)
      iend=iend_arr(n)
      expon=expphi(ntor_arr(n))
      dexpon=expon*dcmplx(0.d0,ntor_arr(n))
      ddexpon=dexpon*dcmplx(0.d0,ntor_arr(n))
!
      call binsrc(flabel_arr(ibeg:iend),ibeg,iend,flabel,k)
!
      k=max(ibeg,k-1)
      dx=flabel-flabel_arr(k)
      ameta=(0.d0,0.d0)
      amtet=(0.d0,0.d0)
      amphi=(0.d0,0.d0)
      sqg_gpp=0.d0
      do i=ncoefs,1,-1
        ameta=ameta*dx+splaeta(:,i,k)
        amtet=amtet*dx+splatet(:,i,k)
        amphi=amphi*dx+splaphi(:,i,k)
        sqg_gpp=sqg_gpp*dx+splsqggpp(i,k)
      enddo
      cdummy=sum(ameta*expthe)
      casrr=casrr+cdummy*expon
      casrr_p=casrr_p+cdummy*dexpon
      casrr_pp=casrr_pp+cdummy*ddexpon
      cdummy=sum(ameta*expthe*dthe)
      casrr_t=casrr_t+cdummy*expon
      casrr_tp=casrr_tp+cdummy*dexpon
      cdummy=sum(ameta*expthe*ddthe)
      casrr_tt=casrr_tt+cdummy*expon
      cdummy=sum(amtet*expthe)
      cat=cat+cdummy*expon
      cat_p=cat_p+cdummy*dexpon
      cat_pp=cat_pp+cdummy*ddexpon
      cdummy=sum(amtet*expthe*dthe)
      cat_t=cat_t+cdummy*expon
      cat_tp=cat_tp+cdummy*dexpon
      cdummy=sum(amtet*expthe*ddthe)
      cat_tt=cat_tt+cdummy*expon
      cdummy=sum(amphi*expthe)
      cap=cap+cdummy*expon
      cap_p=cap_p+cdummy*dexpon
      cap_pp=cap_pp+cdummy*ddexpon
      cdummy=sum(amphi*expthe*dthe)
      cap_t=cap_t+cdummy*expon
      cap_tp=cap_tp+cdummy*dexpon
      cdummy=sum(amphi*expthe*ddthe)
      cap_tt=cap_tt+cdummy*expon
!
      ameta_s=(0.d0,0.d0)
      amtet_s=(0.d0,0.d0)
      amphi_s=(0.d0,0.d0)
      sqg_gpp_s=0.d0
      do i=ncoefs,2,-1
        ameta_s=ameta_s*dx+splaeta(:,i,k)*dfloat(i-1)
        amtet_s=amtet_s*dx+splatet(:,i,k)*dfloat(i-1)
        amphi_s=amphi_s*dx+splaphi(:,i,k)*dfloat(i-1)
        sqg_gpp_s=sqg_gpp_s*dx+splsqggpp(i,k)*dfloat(i-1)
      enddo
      cdummy=sum(ameta_s*expthe)
      casrr_s=casrr_s+cdummy*expon
      casrr_sp=casrr_sp+cdummy*dexpon
      cdummy=sum(ameta_s*expthe*dthe)
      casrr_st=casrr_st+cdummy*expon
      cdummy=sum(amtet_s*expthe)
      cat_s=cat_s+cdummy*expon
      cat_sp=cat_sp+cdummy*dexpon
      cdummy=sum(amtet_s*expthe*dthe)
      cat_st=cat_st+cdummy*expon
      cdummy=sum(amphi_s*expthe)
      cap_s=cap_s+cdummy*expon
      cap_sp=cap_sp+cdummy*dexpon
      cdummy=sum(amphi_s*expthe*dthe)
      cap_st=cap_st+cdummy*expon
!
      ameta_ss=(0.d0,0.d0)
      amtet_ss=(0.d0,0.d0)
      amphi_ss=(0.d0,0.d0)
      sqg_gpp_ss=0.d0
      do i=ncoefs,3,-1
        ameta_ss=ameta_ss*dx+splaeta(:,i,k)*dfloat(i-1)*dfloat(i-2)
        amtet_ss=amtet_ss*dx+splatet(:,i,k)*dfloat(i-1)*dfloat(i-2)
        amphi_ss=amphi_ss*dx+splaphi(:,i,k)*dfloat(i-1)*dfloat(i-2)
        sqg_gpp_ss=sqg_gpp_ss*dx+splsqggpp(i,k)*dfloat(i-1)*dfloat(i-2)
      enddo
      cdummy=sum(ameta_ss*expthe)
      casrr_ss=casrr_ss+cdummy*expon
      cdummy=sum(amtet_ss*expthe)
      cat_ss=cat_ss+cdummy*expon
      cdummy=sum(amphi_ss*expthe)
      cap_ss=cap_ss+cdummy*expon
    enddo
!
    gpp_sqg=1.d0/sqg_gpp
    gpp_sqg_s=-sqg_gpp_s*gpp_sqg**2
    gpp_sqg_ss=2.d0*gpp_sqg_s**2*sqg_gpp-sqg_gpp_ss*gpp_sqg**2
!
    casrr_ss=casrr_ss*gpp_sqg+2.d0*casrr_s*gpp_sqg_s+casrr*gpp_sqg_ss
    casrr_st=casrr_st*gpp_sqg+casrr_t*gpp_sqg_s
    casrr_sp=casrr_sp*gpp_sqg+casrr_p*gpp_sqg_s
    casrr_tt=casrr_tt*gpp_sqg
    casrr_tp=casrr_tp*gpp_sqg
    casrr_pp=casrr_pp*gpp_sqg
    casrr_s=casrr_s*gpp_sqg+casrr*gpp_sqg_s
    casrr_t=casrr_t*gpp_sqg
    casrr_p=casrr_p*gpp_sqg
    casrr=casrr*gpp_sqg
!
    rrinv=1.d0/r**2
    apsi=2.d0*real(casrr)*rrinv
    apsi_s=2.d0*real(casrr_s)*rrinv
    apsi_t=2.d0*real(casrr_t)*rrinv
    apsi_p=2.d0*real(casrr_p)*rrinv
    athe=2.d0*real(cat)
    athe_s=2.d0*real(cat_s)
    athe_t=2.d0*real(cat_t)
    athe_p=2.d0*real(cat_p)
    aphi=2.d0*real(cap)
    aphi_s=2.d0*real(cap_s)
    aphi_t=2.d0*real(cap_t)
    aphi_p=2.d0*real(cap_p)
!
    delar=(apsi-g12*athe)/g11*s_r+athe*theta_r
    delaz=(apsi-g12*athe)/g11*s_z+athe*theta_z
    delbr=((apsi_p-g12*athe_p)/g11*s_z+athe_p*theta_z)/r                      &
         -(aphi_s*s_z+aphi_t*theta_z)/r
    delbz=-((apsi_p-g12*athe_p)/g11*s_r+athe_p*theta_r)/r                     &
         +(aphi_s*s_r+aphi_t*theta_r)/r
    delbp=fcjac*( (apsi_t-g12*athe_t-g12_t*athe)/g11                          &
         -(apsi-g12*athe)*g11_t/g11**2 )                                      &
         +athe_s*(theta_r*s_z-theta_z*s_r)                                    &
         +2.d0*apsi*s_z/(g11*r)
  endif
!
  if(incore.eq.1) then
    Br=Br+delbr
    Bz=Bz+delbz
    Bp=Bp+delbp
  else
    Br=Br+delbr*plaf
    Bz=Bz+delbz*plaf
!    Bp=Bp+delbp*plaf+delar*dpladz-delaz*dpladr
    Bp=Bp+delbp*plaf
  endif
!
!
  end subroutine field_fourier
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine field_fourier_derivs(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                                 ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
! Computes the field and its derivatives using central differences
! for the field components computed by "field_fourier".
!
  implicit none
!
  double precision, parameter :: eps=1.d-5
  double precision :: rrr,ppp,zzz,r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ       &
                     ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,del              &
                     ,rm,zm,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0               &
                     ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0
!
    del=eps*r
!
    rrr=r-del
    zzz=z
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r+del
    zzz=z
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdR=0.5d0*(Br-Br0)/del
    dBpdR=0.5d0*(Bp-Bp0)/del
    dBzdR=0.5d0*(Bz-Bz0)/del
!
    rrr=r
    zzz=z-del
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r
    zzz=z+del
    ppp=phi
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdZ=0.5d0*(Br-Br0)/del
    dBpdZ=0.5d0*(Bp-Bp0)/del
    dBzdZ=0.5d0*(Bz-Bz0)/del
!
    del=eps
!
    rrr=r
    zzz=z
    ppp=phi-eps
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    rrr=r
    zzz=z
    ppp=phi+eps
    call stretch_coords(rrr,zzz,rm,zm)
    rrr=rm
    zzz=zm
    call field_eq(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(rrr,zzz)
    call field_fourier(rrr,ppp,zzz,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0      &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    dBrdp=0.5d0*(Br-Br0)/del
    dBpdp=0.5d0*(Bp-Bp0)/del
    dBzdp=0.5d0*(Bz-Bz0)/del
!
    call field_eq(r,phi,z,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(r,z)
    call field_fourier(r,phi,z,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0          &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
!
  end subroutine field_fourier_derivs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine extract_fluxcoord(phinorm,theta)
!
  use extract_fluxcoord_mod
  use input_files, only : iunit,fluxdatapath
!
  implicit none
!
  integer :: k
  double precision :: phinorm,theta,xpsif
!
  if(load_extract_fluxcoord.eq.1) then
    load_extract_fluxcoord=0
    open(iunit,file=trim(fluxdatapath)//'/phinorm_arr.dat')
    read (iunit,*) nphinorm,psifmin,hpsif
    allocate(phinorm_arr(nphinorm))
    do k=1,nphinorm
      read (iunit,*) phinorm_arr(k)
    enddo
    close(iunit)
  endif
!
  xpsif=(psif_extract-psifmin)/hpsif
  k=min(nphinorm-2,max(0,int(xpsif)))
  phinorm=phinorm_arr(k+1)*(k+1-xpsif)+phinorm_arr(k+2)*(xpsif-k)
!
  theta=theta_extract
!
  end subroutine extract_fluxcoord
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
