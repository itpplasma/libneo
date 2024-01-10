!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine vector_potentials(nr_in,np_in,nz_in,ntor_in,      &
             rmin_in,rmax_in,pmin_in,pmax_in,zmin_in,zmax_in,  &
             br,bp,bz)

  use bdivfree_mod, only : nr,nz,ntor,icp,ipoint,rmin,zmin,hr,hz,pmin,pfac,&
    rpoi,zpoi,apav,rbpav_coef,aznre,aznim,arnre,arnim
  use libneo_kinds, only : complex_kind, real_kind
  use math_constants, only : PI

  implicit none

  integer, intent(in) :: nr_in, np_in, nz_in, ntor_in
  real(kind=real_kind), intent(in) :: rmin_in, rmax_in, pmin_in, pmax_in, &
                                & zmin_in, zmax_in
  real(kind=real_kind), dimension(nr_in,np_in,nz_in), intent(out) :: br, bp, bz

  integer :: ip,np,n,ir,iz
  integer, dimension(:), allocatable :: imi,ima,jmi,jma

  integer :: nashli_rukami
  integer :: irmin, irmax, i,j
  real(kind=real_kind), dimension(4), parameter :: weight=(/-1., 13., 13., -1./)/24.

  real(kind=real_kind) :: hp,r,rm,zm,sumbz,hrm1,hzm1
  real(kind=real_kind), dimension(:),   allocatable :: dummy
  real(kind=real_kind), dimension(:,:), allocatable :: a_re, a_im, rbpav_dummy
  real(kind=real_kind), dimension(:,:), allocatable :: brm,bpm,bzm

  complex(kind=complex_kind) :: four_ampl
  complex(kind=complex_kind), dimension(:,:), allocatable :: expon

  integer, parameter :: mp=4 ! power of Lagrange's polynomial =3
  integer,          dimension(mp)    :: indx,indy
  real(kind=real_kind), dimension(mp)    :: xp,yp
  real(kind=real_kind), dimension(mp,mp) :: fp

  nr=nr_in
  nz=nz_in
  np=np_in-1
  ntor=ntor_in
  nashli_rukami=(nr_in+1)/2

  rmin=rmin_in
  zmin=zmin_in
  hr=(rmax_in-rmin_in)/(nr-1)
  hz=(zmax_in-zmin_in)/(nz-1)
  hp=2.d0*PI/np
  pmin=pmin_in
  pfac = real(nint(2.d0*PI/(pmax_in-pmin_in)), kind=real_kind)

  allocate(expon(np,ntor),a_re(nr,nz),a_im(nr,nz),rbpav_dummy(nr,nz))
  allocate(imi(nz),ima(nz),jmi(nr),jma(nr), dummy(nr))
  allocate(rpoi(nr),zpoi(nz))
  allocate(brm(nr,nz),bpm(nr,nz),bzm(nr,nz))

  imi=1
  ima=nr
  jmi=1
  jma=nz
  do ir=1,nr
    rpoi(ir)=rmin+hr*(ir-1)
  end do
  do iz=1,nz
    zpoi(iz)=zmin+hz*(iz-1)
  end do

  ! Truncation of data outside the limiting convex:
  hrm1=1.d0/hr
  hzm1=1.d0/hz
  do ip=1,np
    do ir=1,nr
      do iz=1,nz
        call stretch_coords(rpoi(ir),zpoi(iz),rm,zm)
        call indef_bdf(rm,rmin,hrm1,nr,indx)
        call indef_bdf(zm,zmin,hzm1,nz,indy)

        do i=1,mp
          xp(i) = rpoi(indx(i))
          yp(i) = zpoi(indy(i))
        end do

        do j=1,mp
          do i=1,mp
            fp(i,j) = Br(indx(i),ip,indy(j))
          end do
        end do
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,brm(ir,iz))
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bp(indx(i),ip,indy(j))
          end do
        end do
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bpm(ir,iz))
        bpm(ir,iz)=bpm(ir,iz)*rm/rpoi(ir)
        do j=1,mp
          do i=1,mp
            fp(i,j) = Bz(indx(i),ip,indy(j))
          end do
        end do
        call plag2d_bdf(rm,zm,fp,hrm1,hzm1,xp,yp,bzm(ir,iz))

      end do
    end do
    Br(:,ip,:)=brm
    Bp(:,ip,:)=bpm
    Bz(:,ip,:)=bzm
  end do
  ! End of data truncation

  allocate(ipoint(nr,nz))
  icp=nr*nz
  allocate(aznre(6,6,icp,ntor),aznim(6,6,icp,ntor))
  allocate(arnre(6,6,icp,ntor),arnim(6,6,icp,ntor))
  allocate(apav(6,6,icp),rbpav_coef(6,6,icp))

  do n=1,ntor
    do ip=1,np
      expon(ip,n)=exp(cmplx(0.d0,-n*(ip-1)*hp, kind=kind(1.0d0)))/np
    end do
  end do

  do n=1,ntor
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(br(ir,1:np,iz)*expon(:,n))*cmplx(0.d0,-r/(n*pfac), kind=kind(1.0d0))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      end do
    end do
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,aznre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,aznim(:,:,:,n),ipoint)
    do ir=1,nr
      r=rmin+hr*(ir-1)
      do iz=1,nz
        four_ampl=sum(bz(ir,1:np,iz)*expon(:,n))*cmplx(0.d0,r/(n*pfac), kind=kind(1.0d0))
        a_re(ir,iz)=dble(four_ampl)
        a_im(ir,iz)=aimag(four_ampl)
      end do
    end do
    call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,arnre(:,:,:,n),ipoint)
    call s2dcut(nr,nz,hr,hz,a_im,imi,ima,jmi,jma,icp,arnim(:,:,:,n),ipoint)
  end do

  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        dummy(ir) = sum(bz(ir,1:np,iz))*hr*r/np
     end do
     a_re(nashli_rukami,iz) = 0.
     sumbz=0.d0
     do ir=nashli_rukami+1,nr
        irmax = min(ir+1,nr)
        irmin = irmax - 3
        sumbz = sumbz + sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     end do
     sumbz=0.d0
     do ir=nashli_rukami-1,1,-1
        irmin = max(ir-1,1)
        irmax = irmin + 3
        sumbz = sumbz - sum(dummy(irmin:irmax)*weight)
        a_re(ir,iz)=sumbz
     end do
  end do

  call s2dcut(nr,nz,hr,hz,a_re,imi,ima,jmi,jma,icp,apav,ipoint)

  do iz=1,nz
     do ir=1,nr
        r=rmin+hr*(ir-1)
        rbpav_dummy(ir,iz) = r*sum(bp(ir,1:np,iz))/np
     end do
  end do

  call s2dcut(nr,nz,hr,hz,rbpav_dummy,imi,ima,jmi,jma,icp,rbpav_coef,ipoint)

  deallocate(expon,a_re,a_im,rbpav_dummy,imi,ima,jmi,jma,dummy,brm,bpm,bzm)

102 format(1000e15.7)
end subroutine vector_potentials


! Spline a single toroidal Fourier mode of the RMP coil field, as needed for MEPHIT.
! Note the explicit-shape input arrays are not checked for shape.
! Also, the truncation of data outside the limiting convex is not performed.
subroutine vector_potential_single_mode(ntor_in, nR_in, nZ_in, Rmin_in, Rmax_in, Zmin_in, Zmax_in, Bn_R, Bn_Z)
  use iso_fortran_env, only: dp => real64
  use bdivfree_mod, only: nR, nZ, Rmin, Zmin, ntor, icp, ipoint, hR, hZ, Rpoi, Zpoi, AZnRe, AZnIm, ARnRe, ARnIm
  integer, intent(in) :: nR_in, nZ_in
  real(dp), intent(in) :: Rmin_in, Rmax_in, Zmin_in, Zmax_in
  complex(dp), intent(in), dimension(nR_in, nZ_in) :: Bn_R, Bn_Z
  integer :: k, imi(nZ_in), ima(nZ_in), jmi(nR_in), jma(nR_in)
  complex(dp), dimension(nR_in, nZ_in) :: An_R, An_Z

  ntor = ntor_in
  nR = nR_in
  nZ = nZ_in
  icp = nR * nZ
  Rmin = Rmin_in
  Zmin = Zmin_in
  hR = (Rmax - Rmin) / dble(nR - 1)
  hZ = (Zmax - Zmin) / dble(nZ - 1)
  allocate(Rpoi(nR), Zpoi(nZ), ipoint(nR, nZ))
  Rpoi(:) = Rmin + [(k * hr, k = 0, nR - 1)]
  Zpoi(:) = Zmin + [(k * hz, k = 0, nZ - 1)]
  do k = 1, nZ
    An_R(:, k) = (0d0, 1d0) * Bn_Z(:, k) * Rpoi / dble(ntor)
    An_Z(:, k) = (0d0, -1d0) * Bn_R(:, k) * Rpoi / dble(ntor)
  end do
  allocate(ARnRe(6, 6, icp, ntor), ARnIm(6, 6, icp, ntor))
  allocate(AZnRe(6, 6, icp, ntor), AZnIm(6, 6, icp, ntor))
  ARnRe(:, :, :, :ntor-1) = (0d0, 0d0)
  ARnIm(:, :, :, :ntor-1) = (0d0, 0d0)
  AZnRe(:, :, :, :ntor-1) = (0d0, 0d0)
  AZnIm(:, :, :, :ntor-1) = (0d0, 0d0)
  imi(:) = 1
  ima(:) = nR
  jmi(:) = 1
  jma(:) = nZ
  call s2dcut(nR, nZ, hR, hZ, An_R%re, imi, ima, jmi, jma, icp, ARnRe(:, :, :, ntor), ipoint)
  call s2dcut(nR, nZ, hR, hZ, An_R%im, imi, ima, jmi, jma, icp, ARnIm(:, :, :, ntor), ipoint)
  call s2dcut(nR, nZ, hR, hZ, An_Z%re, imi, ima, jmi, jma, icp, AZnRe(:, :, :, ntor), ipoint)
  call s2dcut(nR, nZ, hR, hZ, An_Z%im, imi, ima, jmi, jma, icp, AZnIm(:, :, :, ntor), ipoint)
end subroutine vector_potential_single_mode


subroutine spline_vector_potential_n(n, r, z, anr,anz,anr_r,anr_z,anz_r,anz_z, &
  anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz)

  use bdivfree_mod, only : nr,nz,icp,ipoint,hr,hz,&
      & rpoi,zpoi,aznre,aznim,arnre,arnim
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  integer, intent(in) :: n
  real(kind=real_kind), intent(in) :: r, z
  complex(kind=complex_kind), intent(out) :: anr,anz,anr_r,anr_z,anz_r,anz_z
  complex(kind=complex_kind), intent(out) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz

  real(kind=real_kind) :: f,fr,fz,frr,frz,fzz
  real(kind=real_kind) :: g,gr,gz,grr,grz,gzz
  integer :: ierr

  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnim(:,:,:,n),ipoint,r,z,   &
              g,gr,gz,grr,grz,gzz,ierr)
  anr=cmplx(f,g, kind=complex_kind)
  anr_r=cmplx(fr,gr, kind=complex_kind)
  anr_z=cmplx(fz,gz, kind=complex_kind)
  anr_rr=cmplx(frr,grr, kind=complex_kind)
  anr_rz=cmplx(frz,grz, kind=complex_kind)
  anr_zz=cmplx(fzz,gzz, kind=complex_kind)
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznre(:,:,:,n),ipoint,r,z,   &
              f,fr,fz,frr,frz,fzz,ierr)
  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznim(:,:,:,n),ipoint,r,z,   &
              g,gr,gz,grr,grz,gzz,ierr)
  anz=cmplx(f,g, kind=complex_kind)
  anz_r=cmplx(fr,gr, kind=complex_kind)
  anz_z=cmplx(fz,gz, kind=complex_kind)
  anz_rr=cmplx(frr,grr, kind=complex_kind)
  anz_rz=cmplx(frz,grz, kind=complex_kind)
  anz_zz=cmplx(fzz,gzz, kind=complex_kind)

end subroutine spline_vector_potential_n


subroutine spline_bpol_n(n, r, z, B_Rn, B_Zn)
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  integer, intent(in) :: n
  real(kind=real_kind), intent(in) :: r, z
  complex(kind=complex_kind), intent(out) :: B_Rn, B_Zn

  complex(kind=complex_kind) :: anr,anz,anr_r,anr_z,anz_r,anz_z
  complex(kind=complex_kind) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz

  call spline_vector_potential_n(n, r, z, anr,anz,anr_r,anr_z,anz_r,anz_z, &
    anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz)
  B_Rn = (0.d0,1.d0)*dble(n)*anz/r
  B_Zn = -(0.d0,1.d0)*dble(n)*anr/r
end subroutine spline_bpol_n


subroutine spline_bn(n, r, z, Bn_R, Bn_phi, Bn_Z)
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  integer, intent(in) :: n
  real(kind=real_kind), intent(in) :: r, z
  complex(kind=complex_kind), intent(out) :: Bn_R, Bn_phi, Bn_Z

  complex(kind=complex_kind) :: anr,anz,anr_r,anr_z,anz_r,anz_z
  complex(kind=complex_kind) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz

  call spline_vector_potential_n(n, r, z, anr,anz,anr_r,anr_z,anz_r,anz_z, &
    anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz)
  Bn_R = cmplx(0.0,1.0, kind=complex_kind) * real(n, kind=real_kind) * anz / r
  Bn_phi = anr_z - anz_r
  Bn_Z = -cmplx(0.0,1.0, kind=complex_kind) * real(n, kind=real_kind) * anr / r
end subroutine spline_bn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field_divfree(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

  use bdivfree_mod, only : nr,nz,ntor,icp,ipoint,hr,hz,pfac,&
      & rpoi,zpoi,apav,rbpav_coef
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  double precision, intent(in) :: r, phi, z
  double precision, intent(out) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,   &
                              & dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

  integer :: n,ierr
  real(kind=real_kind) :: f,fr,fz,frr,frz,fzz
  real(kind=real_kind) :: delbr,delbz,delbp
  real(kind=real_kind) :: deldBrdR,deldBrdp,deldBrdZ
  real(kind=real_kind) :: deldBpdR,deldBpdp,deldBpdZ
  real(kind=real_kind) :: deldBzdR,deldBzdp,deldBzdZ
  complex(kind=complex_kind) :: expon,anr,anz,anr_r,anr_z,anz_r,anz_z
  complex(kind=complex_kind) :: anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz
  complex(kind=complex_kind) :: pfac_imaginary


  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,rbpav_coef,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)
  Bp = f/r
  dBpdR = fr/r - Bp/r
  dBpdZ = fz/r
  dBpdp = 0.d0

  call spline(nr,nz,rpoi,zpoi,hr,hz,icp,apav,ipoint,r,z,         &
              f,fr,fz,frr,frz,fzz,ierr)

  Br=-fz/r
  Bz=fr/r
  dBrdR=fz/r**2-frz/r
  dBrdZ=-fzz/r
  dBzdR=-fr/r**2+frr/r
  dBzdZ=frz/r
  dBrdp=0.d0
  dBzdp=0.d0

  do n=1,ntor
    call spline_vector_potential_n(n, r, z, anr,anz,anr_r,anr_z,anz_r,anz_z, &
      anr_rr,anr_rz,anr_zz,anz_rr,anz_rz,anz_zz)

    pfac_imaginary = cmplx(0.0, pfac, kind=complex_kind)

    expon=exp(cmplx(0.d0,n*pfac*phi, kind=complex_kind))
    delbr=2.d0*real(pfac_imaginary*n*anz*expon/r, kind=real_kind)
    delbz=-2.d0*real(pfac_imaginary*n*anr*expon/r, kind=real_kind)
    delbp=2.d0*real((anr_z-anz_r)*expon, kind=real_kind)
    deldBrdR=-delbr/r+2.d0*real(pfac_imaginary*n*anz_r*expon/r, kind=real_kind)
    deldBrdZ=2.d0*real(pfac_imaginary*n*anz_z*expon/r, kind=real_kind)
    deldBrdp=-2.d0*real((pfac*n)**2*anz*expon/r, kind=real_kind)
    deldBzdR=-delbz/r-2.d0*real(pfac_imaginary*n*anr_r*expon/r, kind=real_kind)
    deldBzdZ=-2.d0*real(pfac_imaginary*n*anr_z*expon/r, kind=real_kind)
    deldBzdp=2.d0*real((pfac*n)**2*anr*expon/r, kind=real_kind)
    deldBpdR=2.d0*real((anr_rz-anz_rr)*expon, kind=real_kind)
    deldBpdZ=2.d0*real((anr_zz-anz_rz)*expon, kind=real_kind)
    deldBpdp=2.d0*real(pfac_imaginary*n*(anr_z-anz_r)*expon, kind=real_kind)

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

  end do

end subroutine field_divfree

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

  implicit double precision (a-h,o-z)

  ! the power 3 of polinomial is fixed strictly:
  parameter(mp=4)
  integer indu(mp)

  indu(1) = int((u-umin)*dum1)
  if( indu(1) .le. 0 ) indu(1) = 1
  indu(mp) = indu(1) + mp - 1
  if( indu(mp) .gt. nup ) then
    indu(mp) = nup
    indu(1) = indu(mp) - mp + 1
  end if
  do i=2,mp-1
    indu(i) = indu(i-1) + 1
  end do

end subroutine indef_bdf

!---------------------------------------------------------------------
subroutine indsmp_bdf(index_,nup,indu)
  ! defines interval for 1D interpolation on uniform mesh
  ! by known index.
  ! Normally looks for the central interval of stencil, but
  ! stops moving of stencil at the boundary (works for mp=4 only!)
  ! Input:
  !    index - number of a cell on the mesh
  !    nup - total number of mesh points
  ! Output:
  !    indu(mp) - relative index of stencil points

  implicit none

  ! the power 3 of polinomial is fixed strictly:
  integer, parameter :: mp=4

  integer, intent(in) :: index_, nup
  integer, intent(out) :: indu(mp)

  integer :: i

  indu(1) = index_ - 1
  if( indu(1) .le. 0 ) indu(1) = 1
  indu(mp) = indu(1) + mp - 1
  if( indu(mp) .gt. nup ) then
    indu(mp) = nup
    indu(1) = indu(mp) - mp + 1
  end if
  do i=2,mp-1
    indu(i) = indu(i-1) + 1
  end do

end subroutine indsmp_bdf

!---------------------------------------------------------------------
subroutine plag2d_bdf(x,y,fp,dxm1,dym1,xp,yp,polyl2d)
  ! 2D interpolation by means of Lagrange polynomial
  ! uniform mesh (increasingly ordered) in all dimensions is implied
  !
  ! Input parameters:
  !   x,y - coordinates of the point for interpolation
  !   dxm1,dym1 - inverse steps in each direction
  !   xp,yp - vertices of stencil
  !
  ! Output parameters:
  ! polyl2d - polynomial itself

  implicit double precision (a-h,o-z)

  ! the power 3 is fixed strictly:
  parameter(mp=4)

  dimension cx(mp),cy(mp),fp(mp,mp),xp(mp),yp(mp)

  call coefs_bdf(x,xp,dxm1,cx)
  call coefs_bdf(y,yp,dym1,cy)

  polyl2d = 0.d0
  do j=1,mp
    do i=1,mp
      polyl2d = polyl2d + fp(i,j)*cx(i)*cy(j)
    end do
  end do

end subroutine plag2d_bdf

!---------------------------------------------------------------------
subroutine coefs_bdf(u,up,dum1,cu)

  implicit double precision (a-h,o-z)

  parameter(mp=4)
  dimension up(mp),cu(mp)
  data one6/0.16666666666667d0/
  du3 = dum1**3
  cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-one6*du3)
  cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5d0*du3)
  cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5d0*du3)
  cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (one6*du3)

end subroutine coefs_bdf

!---------------------------------------------------------------------
subroutine invert_mono_reg(nx,arry,xmin,xmax,ny,arrx,ymin,ymax)
  ! Inverts the monotonous function y(x) given on the equidistant grid
  ! of x values on the interval [xmin,xmax] by the array y_i=arry(i).
  ! The result, function x(y), is given on the equidistant grid of y values
  ! at the interval [ymin,ymax] by the array x_i=arrx(i).

  use libneo_kinds, only : real_kind

  implicit none

  integer, intent(in) :: nx, ny
  real(kind=real_kind), intent(in) :: xmin, xmax
  real(kind=real_kind), intent(out) :: ymin, ymax
  real(kind=real_kind), dimension(0:nx), intent(in) :: arry
  real(kind=real_kind), dimension(0:ny), intent(out) :: arrx

  integer :: iy,ix,ixfix,ix1,ix2,ix3,ix4

  real(kind=real_kind) :: hy,y,hx,x1,x2,x3,x4,y1,y2,y3,y4

  ixfix = -10

  ymin=arry(0)
  ymax=arry(nx)

  hy=(ymax-ymin)/ny
  hx=(xmax-xmin)/nx

  arrx(0)=xmin
  arrx(ny)=xmax

  do iy=1,ny-1
    y=ymin+iy*hy
    do ix=0,nx
      if(arry(ix).gt.y) then
        ixfix=ix-3
        exit
      end if
    end do
    ixfix=max(ixfix,-1)
    ixfix=min(ixfix,nx-4)
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
  end do

end subroutine invert_mono_reg

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine invert_mono_per(nx,arry_in,xmin,xmax,ny,arrx,ymin,ymax)
  ! Inverts the monotonous function y(x) given on the equidistant grid
  ! of x values on the interval [xmin,xmax] by the array y_i=arry(i).
  ! Special case: y(x) is a sum of linear and periodic functions.
  ! The result, function x(y), is given on the equdistant grid of y values
  ! at the interval [ymin,ymax] by the array x_i=arrx(i).

  use libneo_kinds, only : real_kind

  implicit none

  integer, intent(in) :: nx, ny
  real(kind=real_kind), intent(in) :: xmin, xmax
  real(kind=real_kind), intent(out) :: ymin, ymax
  real(kind=real_kind), dimension(0:nx), intent(in) :: arry_in
  real(kind=real_kind), dimension(0:ny), intent(out) :: arrx

  integer :: iy,ix,ixfix,ix1,ix2,ix3,ix4

  real(kind=real_kind) :: hy,y,hx,x1,x2,x3,x4,y1,y2,y3,y4
  real(kind=real_kind), dimension(:), allocatable :: arry

  ixfix = -10

  allocate(arry(-1:nx+1))
  arry(0:nx)=arry_in

  ymin=arry(0)
  ymax=arry(nx)
  arry(-1)=arry(nx-1)-ymax+ymin
  arry(nx+1)=arry(1)+ymax-ymin

  hy=(ymax-ymin)/ny
  hx=(xmax-xmin)/nx

  arrx(0)=xmin
  arrx(ny)=xmax

  do iy=1,ny-1
    y=ymin+iy*hy
    do ix=0,nx
      if(arry(ix).gt.y) then
        ixfix=ix-3
        exit
      end if
    end do

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
  end do

end subroutine invert_mono_per

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine spl_five_per_bdivfree(n,h,a,b,c,d,e,f)
  ! Periodic spline of the 5-th order. First and last values of function must
  ! be the same.

  use libneo_kinds, only : real_kind

  implicit none

  integer, intent(in) :: n
  real(kind=real_kind), intent(in) :: h
  real(kind=real_kind), dimension(n), intent(in) :: a
  real(kind=real_kind), dimension(n), intent(out) :: b, c, d, e, f

  integer :: i,ip1
  real(kind=real_kind) :: rhop,rhom,fac,xplu,xmin,gammao_m,gammao_p
  real(kind=real_kind) :: c_gammao_m,c_gammao_p

  real(kind=real_kind), dimension(:), allocatable :: alp,bet,gam

  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)

  allocate(alp(n),bet(n),gam(n))

  alp(1)=0.0d0
  bet(1)=0.0d0

  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  end do
  alp(n-2)=-1.d0/(rhop+alp(n-3))
  bet(n-2)=alp(n-2)*(bet(n-3)- &
           5.d0*(a(2)-4.d0*a(1)+6.d0*a(n-1)-4.d0*a(n-2)+a(n-3)))
  alp(n-1)=-1.d0/(rhop+alp(n-2))
  bet(n-1)=alp(n-1)*(bet(n-2)- &
           5.d0*(a(3)-4.d0*a(2)+6.d0*a(1)-4.d0*a(n-1)+a(n-2)))
  alp(n)=-1.d0/(rhop+alp(n-1))
  bet(n)=alp(n)*(bet(n-1)- &
           5.d0*(a(4)-4.d0*a(3)+6.d0*a(2)-4.d0*a(1)+a(n-1)))

  gam(n)=bet(n)
  do i=n-1,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  end do

  xplu=sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
  xmin=-sqrt(0.25d0*rhop**2-1.d0)-0.5d0*rhop
  gammao_m=(gam(2)+xplu*gam(n))/(xmin-xplu)
  gammao_p=(gam(2)+xmin*gam(n))/(xplu-xmin)
  if(abs(xmin).lt.1) then
    c_gammao_m=gammao_m/(xmin**(n-1)-1.d0)
  else
    c_gammao_m=gammao_m*(1.d0/xmin)**(n-1)/(1.d0-(1.d0/xmin)**(n-1))
  end if
  if(abs(xplu).lt.1) then
    c_gammao_p=gammao_p/(xplu**(n-1)-1.d0)
  else
    c_gammao_p=gammao_p*(1.d0/xplu)**(n-1)/(1.d0-(1.d0/xplu)**(n-1))
  end if
  gam(1)=gam(1)+c_gammao_m+c_gammao_p
  do i=2,n
    if(abs(xmin).lt.1) then
      c_gammao_m=gammao_m*xmin**(i-1)/(xmin**(n-1)-1.d0)
    else
      c_gammao_m=gammao_m*(1.d0/xmin)**(n-i)/(1.d0-(1.d0/xmin)**(n-1))
    end if
    if(abs(xplu).lt.1) then
      c_gammao_p=gammao_p*xplu**(i-1)/(xplu**(n-1)-1.d0)
    else
      c_gammao_p=gammao_p*(1.d0/xplu)**(n-i)/(1.d0-(1.d0/xplu)**(n-1))
    end if
    gam(i)=gam(i)+c_gammao_m+c_gammao_p
  end do

  alp(1)=0.0d0
  bet(1)=0.d0

  do i=1,n-1
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  end do

  e(n)=bet(n)
  do i=n-1,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
  end do

  xplu=sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
  xmin=-sqrt(0.25d0*rhom**2-1.d0)-0.5d0*rhom
  gammao_m=(e(2)+xplu*e(n))/(xmin-xplu)
  gammao_p=(e(2)+xmin*e(n))/(xplu-xmin)
  if(abs(xmin).lt.1) then
    c_gammao_m=gammao_m/(xmin**(n-1)-1.d0)
  else
    c_gammao_m=gammao_m*(1.d0/xmin)**(n-1)/(1.d0-(1.d0/xmin)**(n-1))
  end if
  if(abs(xplu).lt.1) then
    c_gammao_p=gammao_p/(xplu**(n-1)-1.d0)
  else
    c_gammao_p=gammao_p*(1.d0/xplu)**(n-1)/(1.d0-(1.d0/xplu)**(n-1))
  end if
  e(1)=e(1)+c_gammao_m+c_gammao_p
  do i=2,n
    if(abs(xmin).lt.1) then
      c_gammao_m=gammao_m*xmin**(i-1)/(xmin**(n-1)-1.d0)
    else
      c_gammao_m=gammao_m*(1.d0/xmin)**(n-i)/(1.d0-(1.d0/xmin)**(n-1))
    end if
    if(abs(xplu).lt.1) then
      c_gammao_p=gammao_p*xplu**(i-1)/(xplu**(n-1)-1.d0)
    else
      c_gammao_p=gammao_p*(1.d0/xplu)**(n-i)/(1.d0-(1.d0/xplu)**(n-1))
    end if
    e(i)=e(i)+c_gammao_m+c_gammao_p
  end do

  do i=n-1,1,-1
    f(i)=(e(i+1)-e(i))/5.d0
  end do
  f(n)=f(1)

  d(n-1)=(a(3)-3.d0*a(2)+3.d0*a(1)-a(n-1))/6.d0 &
      -(e(3)+27.d0*e(2)+93.d0*e(1)+59.d0*e(n-1))/30.d0
  d(n-2)=(a(2)-3.d0*a(1)+3.d0*a(n-1)-a(n-2))/6.d0 &
      -(e(2)+27.d0*e(1)+93.d0*e(n-1)+59.d0*e(n-2))/30.d0
  do i=n-3,1,-1
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
  end do
  d(n)=d(1)
  c(n-1)=0.5d0*(a(2)+a(n-1))-a(1)-0.5d0*d(1)-2.5d0*d(n-1) &
      -0.1d0*(e(2)+18.d0*e(1)+31.d0*e(n-1))
  b(n-1)=a(1)-a(n-1)-c(n-1)-d(n-1)-0.2d0*(4.d0*e(n-1)+e(1))

  do i=n-2,1,-1
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  end do
  b(n)=b(1)
  c(n)=c(1)

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

  deallocate(alp,bet,gam)

end subroutine spl_five_per_bdivfree

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine s2dring(nx,ny,hx,hy,f,icount,spl,ipoint)
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

  implicit double precision (a-h,o-z)

  dimension f(nx,ny),spl(6,6,icount),ipoint(nx,ny)

  integer,          dimension(:), allocatable :: imi,ima,jmi,jma
  double precision, dimension(:), allocatable :: ai,bi,ci,di,ei,fi

  nmax=max(nx,ny)

  allocate( ai(nmax),bi(nmax),ci(nmax),di(nmax),ei(nmax),fi(nmax) )
  allocate(imi(ny),ima(ny),jmi(nx),jma(nx))

  imi=1
  ima=nx
  jmi=1
  jma=ny

  spl=0.d0
  ipoint=-1

  !  spline along Y-axis

  ic = 0
  do i=1,nx
    if(jmi(i).gt.0) then
      nsi=jma(i)-jmi(i)+1
      do j=jmi(i),jma(i)
        ai(j-jmi(i)+1)=f(i,j)
      end do
      call spl_five_per_bdivfree(nsi,hy,ai,bi,ci,di,ei,fi)
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
      end do
    end if
  end do

  if (ic .ne. icount) then
    write (6,*) 'Warning, ic, icount:  ',ic,icount
  endif

  !  spline along X-axis

  do j=1,ny
    if(imi(j).gt.0) then
      nsi=ima(j)-imi(j)+1
      do l=1,6
        do i=imi(j),ima(j)
          ai(i-imi(j)+1)=spl(1,l,ipoint(i,j))
        end do
        call spl_five_reg(nsi,hx,ai,bi,ci,di,ei,fi)
        do i=imi(j),ima(j)
          ii=i-imi(j)+1
          spl(2,l,ipoint(i,j))=bi(ii)
          spl(3,l,ipoint(i,j))=ci(ii)
          spl(4,l,ipoint(i,j))=di(ii)
          spl(5,l,ipoint(i,j))=ei(ii)
          spl(6,l,ipoint(i,j))=fi(ii)
        end do
      end do
    end if
  end do

  deallocate( ai,bi,ci,di,ei,fi )

end subroutine s2dring

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine load_theta

  use theta_rz_mod, only : nsqp,nlab,nthe,icp_pt,ipoint_pt,hsqpsi,hlabel,htheqt,&
    psiaxis,sigma_qt,raxis,zaxis,spllabel,splthet,sqpsi,flab,theqt
  use input_files, only : iunit,fluxdatapath
  use libneo_kinds, only : real_kind
  use math_constants, only : TWOPI

  implicit none

  integer :: nsqpsi,nlabel,ntheta,i
  real(kind=real_kind) :: sqpsimin,sqpsimax
  real(kind=real_kind) :: flabel_min,flabel_max

  real(kind=real_kind), dimension(:),   allocatable :: flabel
  real(kind=real_kind), dimension(:,:), allocatable :: theta_of_theta_qt

  open(iunit,form='unformatted',                                 &
       file=trim(fluxdatapath)//'/theta_of_theta_qt_flabel.dat')
  read (iunit) nsqpsi,nlabel,ntheta,sqpsimin,sqpsimax,flabel_min,flabel_max &
              ,raxis,zaxis,psiaxis,sigma_qt
  allocate(theta_of_theta_qt(nlabel,0:ntheta),flabel(0:nsqpsi))
  read (iunit) theta_of_theta_qt
  read (iunit) flabel
  close(iunit)

  nsqp=nsqpsi
  nlab=nlabel
  nthe=ntheta+1

  hsqpsi=(sqpsimax-sqpsimin)/(nsqp-1)
  hlabel=(flabel_max-flabel_min)/(nlab-1)
  htheqt=TWOPI/ntheta

  allocate(sqpsi(nsqp),flab(nlab),theqt(nthe))

  do i=1,nsqp
    sqpsi(i)=sqpsimin+hsqpsi*(i-1)
  end do

  do i=1,nlab
    flab(i)=flabel_min+hlabel*(i-1)
  end do

  do i=1,nthe
    theqt(i)=htheqt*(i-1)
  end do

  icp_pt=nthe*nlab
  allocate( splthet(6,6,icp_pt), ipoint_pt(nlab,nthe) )

  call s2dring(nlab,nthe,hlabel,htheqt,theta_of_theta_qt(:,0:ntheta), &
               icp_pt,splthet,ipoint_pt)

  allocate(spllabel(6,nsqpsi))
  spllabel(1,:)=flabel(1:nsqpsi)
  call spl_five_reg(nsqpsi,hsqpsi,spllabel(1,:),spllabel(2,:),spllabel(3,:) &
                   ,              spllabel(4,:),spllabel(5,:),spllabel(6,:))

end subroutine load_theta

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine psithet_rz(rrr,zzz,                                          &
                        theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                        flabel,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)

  use theta_rz_mod, only : icall, nsqp,nlab,nthe,icp_pt,ipoint_pt,hsqpsi,hlabel,&
    htheqt,psiaxis,sigma_qt,raxis,zaxis,spllabel,splthet,flab,theqt
  use field_eq_mod, only : nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint  &
                         , psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use extract_fluxcoord_mod, only : psif_extract,theta_extract
  use libneo_kinds, only : real_kind
  use math_constants, only : TWOPI

  implicit none

  real(kind=real_kind), intent(in) :: rrr, zzz
  real(kind=real_kind), intent(out) :: theta, theta_r, theta_z, theta_rr, &
                            & theta_rz, theta_zz
  real(kind=real_kind), intent(out) :: flabel
  real(kind=real_kind), intent(out) :: s_r, s_z, s_rr, s_rz, s_zz
  real(kind=real_kind), intent(out) :: s0, ds0ds, dds0ds

  integer :: ierr,k
  real(kind=real_kind) :: theta_s,theta_t,theta_ss,theta_st,theta_tt
  real(kind=real_kind) :: sqpsi_qt
  real(kind=real_kind) :: theta_qt,t_r,t_z,t_rr,t_rz,t_zz
  real(kind=real_kind) :: rho2,rho4,dr,dz,dflabel,ddflabel,dx,dfl_dpsi,ddfl_dpsi

  if(icall.eq.0) then
    icall=1
    call load_theta
  endif

  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz, &
              psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

  sqpsi_qt=sqrt(abs(psif-psiaxis))

  k=min(int(sqpsi_qt/hsqpsi),nsqp)
  dx=sqpsi_qt-hsqpsi*k
  flabel=spllabel(1,k)+dx*(spllabel(2,k)+dx*(spllabel(3,k)           &
        +dx*(spllabel(4,k)+dx*(spllabel(5,k)+dx*spllabel(6,k)))))
  dflabel=spllabel(2,k)+dx*(2.d0*spllabel(3,k)                       &
          +dx*(3.d0*spllabel(4,k)+dx*(4.d0*spllabel(5,k)             &
          +dx*5.d0*spllabel(6,k))))
  ddflabel=2.d0*spllabel(3,k)+dx*(6.d0*spllabel(4,k)                 &
           +dx*(12.d0*spllabel(5,k)+dx*20.d0*spllabel(6,k)))

  dfl_dpsi=sign(0.5d0,psif-psiaxis)*dflabel/sqpsi_qt
  ddfl_dpsi=0.25d0*(ddflabel-dflabel/sqpsi_qt)/abs(psif-psiaxis)
  s_r=dpsidr*dfl_dpsi
  s_z=dpsidz*dfl_dpsi
  s_rr=d2psidr2*dfl_dpsi+dpsidr**2*ddfl_dpsi
  s_rz=d2psidrdz*dfl_dpsi+dpsidr*dpsidz*ddfl_dpsi
  s_zz=d2psidz2*dfl_dpsi+dpsidz**2*ddfl_dpsi

  s0=sqpsi_qt
  ds0ds=1.d0/dflabel
  dds0ds=-ds0ds**3*ddflabel

  dr=rrr-raxis
  dz=zzz-zaxis
  rho2=dr**2+dz**2
  rho4=rho2**2
  theta_qt=mod(sigma_qt*atan2(dz,dr)+TWOPI, TWOPI)
  t_r=-sigma_qt*dz/rho2
  t_z=sigma_qt*dr/rho2
  t_rr=2.d0*sigma_qt*dr*dz/rho4
  t_zz=-t_rr
  t_rz=sigma_qt*(dz**2-dr**2)/rho4

  call spline(nlab,nthe,flab,theqt,hlabel,htheqt,icp_pt,splthet,ipoint_pt, &
              flabel,theta_qt,                                             &
              theta,theta_s,theta_t,theta_ss,theta_st,theta_tt,ierr)

  theta=theta+theta_qt
  theta_r=theta_s*s_r+(theta_t+1.d0)*t_r
  theta_z=theta_s*s_z+(theta_t+1.d0)*t_z
  theta_rr=theta_ss*s_r**2+2.d0*theta_st*s_r*t_r+theta_tt*t_r**2 &
          +theta_s*s_rr+(theta_t+1.d0)*t_rr
  theta_rz=theta_ss*s_r*s_z+theta_st*(s_r*t_z+s_z*t_r)+theta_tt*t_r*t_z &
          +theta_s*s_rz+(theta_t+1.d0)*t_rz
  theta_zz=theta_ss*s_z**2+2.d0*theta_st*s_z*t_z+theta_tt*t_z**2 &
          +theta_s*s_zz+(theta_t+1.d0)*t_zz

  psif_extract=psif
  theta_extract=theta

end subroutine psithet_rz

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cspl_five_reg(n,h,a,b,c,d,e,f)
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  integer, intent(in) :: n
  real(kind=real_kind), intent(in) :: h
  complex(kind=complex_kind), dimension(n), intent(in) :: a
  complex(kind=complex_kind), dimension(n), intent(out) :: b, c, d, e, f

  integer :: i,ip1
  real(kind=real_kind) :: rhop,rhom,fac
  real(kind=real_kind) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,det
  complex(kind=complex_kind) :: abeg,bbeg,cbeg,dbeg,ebeg,fbeg
  complex(kind=complex_kind) :: aend,bend,cend,dend,eend,fend
  complex(kind=complex_kind) :: b1,b2,b3
  complex(kind=complex_kind), dimension(:), allocatable :: alp,bet,gam

  rhop=13.d0+sqrt(105.d0)
  rhom=13.d0-sqrt(105.d0)

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

  allocate(alp(n),bet(n),gam(n))

  alp(1)=0.0d0
  bet(1)=ebeg*(2.d0+rhom)-5.d0*fbeg*(3.d0+1.5d0*rhom) !gamma1

  do i=1,n-4
    ip1=i+1
    alp(ip1)=-1.d0/(rhop+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)- &
             5.d0*(a(i+4)-4.d0*a(i+3)+6.d0*a(i+2)-4.d0*a(ip1)+a(i)))
  end do

  gam(n-2)=eend*(2.d0+rhom)+5.d0*fend*(3.d0+1.5d0*rhom) !gamma
  do i=n-3,1,-1
    gam(i)=gam(i+1)*alp(i)+bet(i)
  end do

  alp(1)=0.0d0
  bet(1)=ebeg-2.5d0*5.d0*fbeg !e1

  do i=1,n-2
    ip1=i+1
    alp(ip1)=-1.d0/(rhom+alp(i))
    bet(ip1)=alp(ip1)*(bet(i)-gam(i))
  end do

  e(n)=eend+2.5d0*5.d0*fend
  e(n-1)=e(n)*alp(n-1)+bet(n-1)
  f(n-1)=(e(n)-e(n-1))/5.d0
  e(n-2)=e(n-1)*alp(n-2)+bet(n-2)
  f(n-2)=(e(n-1)-e(n-2))/5.d0
  d(n-2)=dend+1.5d0*4.d0*eend+1.5d0**2*10.d0*fend

  do i=n-3,1,-1
    e(i)=e(i+1)*alp(i)+bet(i)
    f(i)=(e(i+1)-e(i))/5.d0
    d(i)=(a(i+3)-3.d0*a(i+2)+3.d0*a(i+1)-a(i))/6.d0 &
        -(e(i+3)+27.d0*e(i+2)+93.d0*e(i+1)+59.d0*e(i))/30.d0
    c(i)=0.5d0*(a(i+2)+a(i))-a(i+1)-0.5d0*d(i+1)-2.5d0*d(i) &
        -0.1d0*(e(i+2)+18.d0*e(i+1)+31.d0*e(i))
    b(i)=a(i+1)-a(i)-c(i)-d(i)-0.2d0*(4.d0*e(i)+e(i+1))
  end do

  do i=n-3,n
    b(i)=b(i-1)+2.d0*c(i-1)+3.d0*d(i-1)+4.d0*e(i-1)+5.d0*f(i-1)
    c(i)=c(i-1)+3.d0*d(i-1)+6.d0*e(i-1)+10.d0*f(i-1)
    d(i)=d(i-1)+4.d0*e(i-1)+10.d0*f(i-1)
    if(i.ne.n) f(i)= a(i+1)-a(i)-b(i)-c(i)-d(i)-e(i)
  end do
  f(n)=f(n-1)

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

  deallocate(alp,bet,gam)

end subroutine cspl_five_reg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field_fourier(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ              &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  ! Caution: derivatives are not computed, for derivatives call
  ! a driver routine "field_fourier_derivs"

  use amn_mod, only : ntor_amn,mpol,ntor_ff,mpol_ff,nsqpsi,icall,&
    sqpsimin,sqpsimax,hsqpsi,splapsi,splatet,amnpsi,amntet,&
    amnpsi_s,amntet_s,amnpsi_ss,amntet_ss,expthe,expphi,&
    nsqpsi_ff,nmodes_ff,sqpsimin_ff,sqpsimax_ff,hsqpsi_ff,&
    ipoi_ff,splffp,splfft,fmnpsi,fmntet,fmnpsi_s, fmntet_s,&
    fmnpsi_ss,fmntet_ss
  use input_files,           only : iunit,fluxdatapath
  use inthecore_mod, only : incore,psi_sep                                 &
                          , plaf,dpladr,dpladz
  use field_eq_mod,  only : dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use theta_rz_mod,  only : psiaxis
  use bdivfree_mod,  only : pfac
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  real(kind=real_kind), intent(in) :: r, phi, z
  real(kind=real_kind), intent(out) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ,    &
                       &  dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

  integer :: m,n,i,k,ntor
  real(kind=real_kind) :: sqpsi,dx,g11,g12,g11_r,g11_z,g12_r,g12_z
  real(kind=real_kind) :: theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                      s_r,s_z,s_rr,s_rz,s_zz
  real(kind=real_kind) :: apsi,apsi_s,apsi_t,apsi_p
  real(kind=real_kind) :: athe,athe_s,athe_t,athe_p
  real(kind=real_kind) :: delbr,delbz,delbp,delar,delaz
  real(kind=real_kind) :: fcjac,g11_t,g12_t,s0,ds0ds,dds0ds,sqpsi_sep

  integer, dimension(:,:), allocatable :: idummy2

  complex(kind=complex_kind) :: expon
  complex(kind=complex_kind), dimension(:), allocatable :: a,b,c,d,e,f
  complex(kind=complex_kind), dimension(:,:,:), allocatable :: apsimn,athetmn

  integer, save :: nper

  ! Initialization ------------------------------------------------------------

  if(icall.eq.0) then
    icall=1

    nper=nint(pfac)
    print *,'nper = ',nper
    ! Toroidally symetric part of the vacuum perturbation field - comes now
    ! from the cylindrical vacuum field routine


    ! Fourier ampitudes of the non-axisymmetric vacuum perturbation field:

    open(iunit,form='unformatted',file=trim(fluxdatapath)//'/amn.dat')
    read (iunit) ntor,mpol,nsqpsi,sqpsimin,sqpsimax
    allocate(apsimn(-mpol:mpol,ntor,nsqpsi))
    allocate(athetmn(-mpol:mpol,ntor,nsqpsi))
    read (iunit) apsimn,athetmn
    close(iunit)

    call psithet_rz(r,z,                                              &
                    theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                    sqpsi,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)

    nsqpsi_ff=0
    open(iunit,form='unformatted',file=trim(fluxdatapath)//'/formfactors.dat')
    read (iunit,end=1) nmodes_ff,nsqpsi_ff,mpol_ff,mpol_ff,ntor_ff,ntor_ff
    read (iunit) sqpsimin_ff,sqpsimax_ff
    print *,'nsqpsi_ff = ',nsqpsi_ff
    if(ntor_ff.ne.ntor.or.mpol_ff.ne.mpol) then
      print *,'Number of harmonics in formfactors differs from original field'
      print *,'ntor = ',ntor,'ntor_ff = ',ntor_ff, &
              'mpol = ',mpol,'mpol_ff = ',mpol_ff
    end if
    allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
    allocate(splfft(nmodes_ff,6,nsqpsi_ff))
    read (iunit) ipoi_ff(-mpol_ff:mpol_ff,1:ntor_ff)
    read (iunit) splfft(1:nmodes_ff,1,1:nsqpsi_ff)

1   close(iunit)

    if(nsqpsi_ff.gt.0) then

      call smear_formfactors(nmodes_ff,nsqpsi_ff,sqpsimin_ff,sqpsimax_ff, &
                             splfft(1:nmodes_ff,1,1:nsqpsi_ff))

      ! use those formfactors which available and necessary
      mpol_ff=min(mpol,mpol_ff)
      ntor_ff=min(ntor,ntor_ff)
      allocate(splffp(nmodes_ff,6,nsqpsi_ff))
      splffp(:,1,:)=splfft(:,1,:)
      allocate(idummy2(-mpol_ff:mpol_ff,ntor_ff))
      idummy2=ipoi_ff(-mpol_ff:mpol_ff,1:ntor_ff)
      deallocate(ipoi_ff)
      allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
      ipoi_ff=idummy2
      deallocate(idummy2)
    else
      print *,'Formfactor file formfactors.dat empty or absent,' &
             //' compute vacuum field'
      nmodes_ff=1
      mpol_ff=1
      ntor_ff=1
      nsqpsi_ff=10
      sqpsimin_ff=0.d0
      sqpsimax_ff=sqrt(abs(psi_sep-psiaxis))
      allocate(ipoi_ff(-mpol_ff:mpol_ff,ntor_ff))
      ipoi_ff=-1
      allocate(splffp(nmodes_ff,6,nsqpsi_ff))
      allocate(splfft(nmodes_ff,6,nsqpsi_ff))
      splffp(:,1,:)=1.d0
      splfft(:,1,:)=splffp(:,1,:)
    end if

    mpol=mpol_ff
    ntor=ntor_ff
    ntor_amn=ntor

    allocate(splapsi(-mpol:mpol,ntor,6,nsqpsi))
    allocate(splatet(-mpol:mpol,ntor,6,nsqpsi))

    allocate(amnpsi(-mpol:mpol,ntor),amntet(-mpol:mpol,ntor))
    allocate(amnpsi_s(-mpol:mpol,ntor),amntet_s(-mpol:mpol,ntor))
    allocate(amnpsi_ss(-mpol:mpol,ntor),amntet_ss(-mpol:mpol,ntor))

    allocate(expthe(-mpol:mpol),expphi(ntor))

    hsqpsi=(sqpsimax-sqpsimin)/(nsqpsi-1)

    allocate(a(nsqpsi),b(nsqpsi),c(nsqpsi),d(nsqpsi),e(nsqpsi),f(nsqpsi))

    do m=-mpol,mpol
      do n=1,ntor
        a=apsimn(m,n,:)
        call cspl_five_reg(nsqpsi,hsqpsi,a,b,c,d,e,f)
        splapsi(m,n,1,:)=a
        splapsi(m,n,2,:)=b
        splapsi(m,n,3,:)=c
        splapsi(m,n,4,:)=d
        splapsi(m,n,5,:)=e
        splapsi(m,n,6,:)=f
        a=athetmn(m,n,:)
        call cspl_five_reg(nsqpsi,hsqpsi,a,b,c,d,e,f)
        splatet(m,n,1,:)=a
        splatet(m,n,2,:)=b
        splatet(m,n,3,:)=c
        splatet(m,n,4,:)=d
        splatet(m,n,5,:)=e
        splatet(m,n,6,:)=f
      end do
    end do

    ! Formfactors:

    allocate(fmnpsi(nmodes_ff))
    allocate(fmntet(nmodes_ff))
    allocate(fmnpsi_s(nmodes_ff))
    allocate(fmntet_s(nmodes_ff))
    allocate(fmnpsi_ss(nmodes_ff))
    allocate(fmntet_ss(nmodes_ff))

    hsqpsi_ff=(sqpsimax_ff-sqpsimin_ff)/(nsqpsi_ff-1)

    deallocate(a,b,c,d,e,f)
    allocate(a(nsqpsi_ff),b(nsqpsi_ff),c(nsqpsi_ff))
    allocate(d(nsqpsi_ff),e(nsqpsi_ff),f(nsqpsi_ff))

    do i=1,nmodes_ff
      a=splffp(i,1,:)
      call cspl_five_reg(nsqpsi_ff,hsqpsi_ff,a,b,c,d,e,f)
      splffp(i,1,:)=a
      splffp(i,2,:)=b
      splffp(i,3,:)=c
      splffp(i,4,:)=d
      splffp(i,5,:)=e
      splffp(i,6,:)=f
      a=splfft(i,1,:)
      call cspl_five_reg(nsqpsi_ff,hsqpsi_ff,a,b,c,d,e,f)
      splfft(i,1,:)=a
      splfft(i,2,:)=b
      splfft(i,3,:)=c
      splfft(i,4,:)=d
      splfft(i,5,:)=e
      splfft(i,6,:)=f
    end do

    ! Normalize formfactors to 1 at the separatrix:
    sqpsi_sep=sqrt(abs(psi_sep-psiaxis))
    k=min(nsqpsi_ff,max(1,ceiling((sqpsi_sep-sqpsimin_ff)/hsqpsi_ff)))
    dx=sqpsi_sep-sqpsimin_ff-hsqpsi_ff*(k-1)

    fmnpsi=splffp(:,1,k)+dx*(splffp(:,2,k)+dx*(splffp(:,3,k)        &
          +dx*(splffp(:,4,k)+dx*(splffp(:,5,k)+dx*splffp(:,6,k)))))
    fmntet=splfft(:,1,k)+dx*(splfft(:,2,k)+dx*(splfft(:,3,k)        &
          +dx*(splfft(:,4,k)+dx*(splfft(:,5,k)+dx*splfft(:,6,k)))))
    do i=1,6
      do k=1,nsqpsi_ff
        splffp(:,i,k)=splffp(:,i,k)/fmnpsi
        splfft(:,i,k)=splfft(:,i,k)/fmntet
      end do
    end do
    splffp(:,1,:)=splffp(:,1,:)-1.d0
    splfft(:,1,:)=splfft(:,1,:)-1.d0

  end if

  ! End of initialization ------------------------------------------------------

  ! Toroidally symmetric part of the perturbation field - not computed, comes
  ! from the vacuum routine

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

  ! Asymmetric part of the perturbation field:

  call psithet_rz(r,z,                                              &
                  theta,theta_r,theta_z,theta_rr,theta_rz,theta_zz, &
                  sqpsi,s_r,s_z,s_rr,s_rz,s_zz,s0,ds0ds,dds0ds)

  g11=dpsidr**2+dpsidz**2
  g12=dpsidr*theta_r+dpsidz*theta_z
  g11_r=2.d0*dpsidr*d2psidr2+2.d0*dpsidz*d2psidrdz
  g11_z=2.d0*dpsidr*d2psidrdz+2.d0*dpsidz*d2psidz2
  g12_r=d2psidr2*theta_r+dpsidr*theta_rr+d2psidrdz*theta_z+dpsidz*theta_rz
  g12_z=d2psidrdz*theta_r+dpsidr*theta_rz+d2psidz2*theta_z+dpsidz*theta_zz
  fcjac=dpsidr*theta_z-dpsidz*theta_r
  g11_t=(dpsidr*g11_z-dpsidz*g11_r)/fcjac
  g12_t=(dpsidr*g12_z-dpsidz*g12_r)/fcjac

  k=min(nsqpsi,max(1,ceiling((sqpsi-sqpsimin)/hsqpsi)))
  dx=sqpsi-sqpsimin-hsqpsi*(k-1)

  amnpsi=splapsi(:,:,1,k)+dx*(splapsi(:,:,2,k)+dx*(splapsi(:,:,3,k)       &
        +dx*(splapsi(:,:,4,k)+dx*(splapsi(:,:,5,k)+dx*splapsi(:,:,6,k)))))
  amntet=splatet(:,:,1,k)+dx*(splatet(:,:,2,k)+dx*(splatet(:,:,3,k)       &
        +dx*(splatet(:,:,4,k)+dx*(splatet(:,:,5,k)+dx*splatet(:,:,6,k)))))
  amnpsi_s=splapsi(:,:,2,k)+dx*(2.d0*splapsi(:,:,3,k)                     &
          +dx*(3.d0*splapsi(:,:,4,k)+dx*(4.d0*splapsi(:,:,5,k)            &
          +dx*5.d0*splapsi(:,:,6,k))))
  amntet_s=splatet(:,:,2,k)+dx*(2.d0*splatet(:,:,3,k)                     &
          +dx*(3.d0*splatet(:,:,4,k)+dx*(4.d0*splatet(:,:,5,k)            &
          +dx*5.d0*splatet(:,:,6,k))))
  amnpsi_ss=2.d0*splapsi(:,:,3,k)+dx*(6.d0*splapsi(:,:,4,k)               &
           +dx*(12.d0*splapsi(:,:,5,k)+dx*20.d0*splapsi(:,:,6,k)))
  amntet_ss=2.d0*splatet(:,:,3,k)+dx*(6.d0*splatet(:,:,4,k)               &
           +dx*(12.d0*splatet(:,:,5,k)+dx*20.d0*splatet(:,:,6,k)))

  ! Formfactors:

  k=min(nsqpsi_ff,max(1,ceiling((s0-sqpsimin_ff)/hsqpsi_ff)))
  dx=s0-sqpsimin_ff-hsqpsi_ff*(k-1)

  fmnpsi=splffp(:,1,k)+dx*(splffp(:,2,k)+dx*(splffp(:,3,k)        &
        +dx*(splffp(:,4,k)+dx*(splffp(:,5,k)+dx*splffp(:,6,k)))))
  fmntet=splfft(:,1,k)+dx*(splfft(:,2,k)+dx*(splfft(:,3,k)        &
        +dx*(splfft(:,4,k)+dx*(splfft(:,5,k)+dx*splfft(:,6,k)))))
  fmnpsi_s=splffp(:,2,k)+dx*(2.d0*splffp(:,3,k)                     &
          +dx*(3.d0*splffp(:,4,k)+dx*(4.d0*splffp(:,5,k)            &
          +dx*5.d0*splffp(:,6,k))))
  fmntet_s=splfft(:,2,k)+dx*(2.d0*splfft(:,3,k)                     &
          +dx*(3.d0*splfft(:,4,k)+dx*(4.d0*splfft(:,5,k)            &
          +dx*5.d0*splfft(:,6,k))))
  fmnpsi_ss=2.d0*splffp(:,3,k)+dx*(6.d0*splffp(:,4,k)               &
           +dx*(12.d0*splffp(:,5,k)+dx*20.d0*splffp(:,6,k)))
  fmntet_ss=2.d0*splfft(:,3,k)+dx*(6.d0*splfft(:,4,k)               &
           +dx*(12.d0*splfft(:,5,k)+dx*20.d0*splfft(:,6,k)))

  ! convert forfactor derivatives to derivatives over new label:

  fmnpsi_ss=fmnpsi_ss*ds0ds**2+fmnpsi_s*dds0ds
  fmnpsi_s=fmnpsi_s*ds0ds
  fmntet_ss=fmntet_ss*ds0ds**2+fmntet_s*dds0ds
  fmntet_s=fmntet_s*ds0ds

  ! Product:

  do m=-mpol_ff,mpol_ff
    do n=1,ntor_ff
      if(n*nper.gt.ntor_ff) cycle
      if(ipoi_ff(m,n*nper).gt.0) then
        k=ipoi_ff(m,n*nper)
        amnpsi_ss(m,n)=amnpsi_ss(m,n)*fmnpsi(k)         &
                      +2.d0*amnpsi_s(m,n)*fmnpsi_s(k)   &
                      +amnpsi(m,n)*fmnpsi_ss(k)
        amntet_ss(m,n)=amntet_ss(m,n)*fmntet(k)         &
                      +2.d0*amntet_s(m,n)*fmntet_s(k)   &
                      +amntet(m,n)*fmntet_ss(k)
        amnpsi_s(m,n) =amnpsi_s(m,n)*fmnpsi(k)          &
                      +amnpsi(m,n)*fmnpsi_s(k)
        amntet_s(m,n) =amntet_s(m,n)*fmntet(k)          &
                      +amntet(m,n)*fmntet_s(k)
        amnpsi(m,n)   =amnpsi(m,n)*fmnpsi(k)
        amntet(m,n)   =amntet(m,n)*fmntet(k)
      end if
    end do
  end do

  expthe(0)=(1.d0,0.d0)
  expthe(1)=exp(cmplx(0.d0,theta, kind=kind(1.0d0)))
  expthe(-1)=conjg(expthe(1))
  do m=2,mpol
    expthe(m)=expthe(m-1)*expthe(1)
    expthe(-m)=conjg(expthe(m))
  end do

  expphi(1)=exp(cmplx(0.d0,pfac*phi, kind=kind(1.0d0)))
  do n=2,ntor_amn
    expphi(n)=expphi(n-1)*expphi(1)
  end do

  apsi=0.d0
  apsi_s=0.d0
  apsi_t=0.d0
  apsi_p=0.d0
  athe=0.d0
  athe_s=0.d0
  athe_t=0.d0
  athe_p=0.d0
  do m=-mpol,mpol
    do n=1,ntor_amn
      if(n*nper.gt.ntor_ff) cycle
      if(ipoi_ff(m,n*nper).gt.0) then
        expon=expthe(m)*expphi(n)
        apsi=apsi+2.d0*dble(expon*amnpsi(m,n))
        apsi_s=apsi_s+2.d0*dble(expon*amnpsi_s(m,n))
        apsi_t=apsi_t+2.d0*dble((0.d0,1.d0)*m*expon*amnpsi(m,n))
        apsi_p=apsi_p+2.d0*dble((0.d0,1.d0)*n*expon*amnpsi(m,n))*pfac
        athe=athe+2.d0*dble(expon*amntet(m,n))
        athe_s=athe_s+2.d0*dble(expon*amntet_s(m,n))
        athe_t=athe_t+2.d0*dble((0.d0,1.d0)*m*expon*amntet(m,n))
        athe_p=athe_p+2.d0*dble((0.d0,1.d0)*n*expon*amntet(m,n))*pfac
      end if
    end do
  end do

  delar=(apsi-g12*athe)/g11*dpsidr+athe*theta_r
  delaz=(apsi-g12*athe)/g11*dpsidz+athe*theta_z
  delbr=((apsi_p-g12*athe_p)/g11*dpsidz+athe_p*theta_z)/r
  delbz=-((apsi_p-g12*athe_p)/g11*dpsidr+athe_p*theta_r)/r
  delbp=fcjac*( (apsi_t-g12*athe_t-g12_t*athe)/g11  &
       -(apsi-g12*athe)*g11_t/g11**2 )              &
       +athe_s*(theta_r*s_z-theta_z*s_r)

  if(incore.eq.1) then
    Br=Br+delbr
    Bz=Bz+delbz
    Bp=Bp+delbp
  else
    Br=Br+delbr*plaf
    Bz=Bz+delbz*plaf
    Bp=Bp+delbp*plaf+delar*dpladz-delaz*dpladr
  end if

end subroutine field_fourier


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field_fourier_derivs(r,phi,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ    &
                                 ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  ! Computes the field and its derivatives using central differences
  ! for the field components computed by "field_fourier".

  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind), parameter :: eps=1.d-7

  real(kind=real_kind), intent(in) :: phi, r, z
  real(kind=real_kind), intent(out) :: Br, Bp, Bz, dBrdR, dBrdp, dBrdZ
  real(kind=real_kind), intent(out) :: dBpdR, dBpdp, dBpdZ,dBzdR, dBzdp, dBzdZ

  real(kind=real_kind) :: rrr,ppp,zzz,del                              &
                     ,rm,zm,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0       &
                     ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0

    del=eps*r

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

    del=eps

    rrr=r
    zzz=z
    ppp=phi-del
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
    ppp=phi+del
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

    call field_eq(r,phi,z,Br0,Bp0,Bz0,dBrdR0,dBrdp0,dBrdZ0   &
                  ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)
    call inthecore(r,z)
    call field_fourier(r,phi,z,Br,Bp,Bz,dBrdR0,dBrdp0,dBrdZ0          &
                      ,dBpdR0,dBpdp0,dBpdZ0,dBzdR0,dBzdp0,dBzdZ0)

end subroutine field_fourier_derivs

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine extract_fluxcoord(phinorm,theta)

  use extract_fluxcoord_mod, only : load_extract_fluxcoord,nphinorm,&
    psif_extract,theta_extract,psifmin,hpsif,phinorm_arr
  use input_files, only : iunit,fluxdatapath
  use libneo_kinds, only : real_kind

  implicit none

  real(kind=real_kind), intent(out) :: phinorm,theta

  integer :: k
  real(kind=real_kind) :: xpsif

  if(load_extract_fluxcoord.eq.1) then
    load_extract_fluxcoord=0
    open(iunit,file=trim(fluxdatapath)//'/phinorm_arr.dat')
    read (iunit,*) nphinorm,psifmin,hpsif
    allocate(phinorm_arr(nphinorm))
    do k=1,nphinorm
      read (iunit,*) phinorm_arr(k)
    end do
    close(iunit)
  end if

  xpsif=(psif_extract-psifmin)/hpsif
  k=min(nphinorm-2,max(0,int(xpsif)))
  phinorm=phinorm_arr(k+1)*(k+1-xpsif)+phinorm_arr(k+2)*(xpsif-k)

  theta=theta_extract

end subroutine extract_fluxcoord

!> Since the field in the main volume and outside separatrix are given
!> by different models, it smoothly connects these two fields in some
!> narrow transition layer located in the main volume adjacent to the
!> separatrix. It should have an effect on the field only if one
!> computes something inside this transition layer which is not so easy
!> to hit by chance.
subroutine smear_formfactors(nmodes_ff,nsqpsi_ff,sqpsimin_ff,sqpsimax_ff, &
                               formfactors)

  use inthecore_mod, only : psi_sep,psi_cut
  use libneo_kinds, only : complex_kind, real_kind
  use theta_rz_mod,  only : psiaxis

  implicit none

  integer, intent(in) :: nmodes_ff,nsqpsi_ff
  real(kind=real_kind), intent(in) :: sqpsimin_ff,sqpsimax_ff
  complex(kind=complex_kind), dimension(nmodes_ff,nsqpsi_ff), intent(inout) :: formfactors

  real(kind=real_kind) :: hsqpsi_ff,apsif
  real(kind=real_kind) :: apsi_sep,apsi_cut,weight,dweight,ddweight

  integer :: i
  real(kind=real_kind) :: R=1.d0,Z=0.d0

  call inthecore(R,Z)

  hsqpsi_ff=(sqpsimax_ff-sqpsimin_ff)/dble(nsqpsi_ff-1)
  apsi_sep=abs(psi_sep-psiaxis)
  apsi_cut=abs(psi_cut-psiaxis)

  do i=1,nsqpsi_ff
    apsif=(sqpsimin_ff+hsqpsi_ff*dble(i-1))**2
    call localizer(apsi_cut,apsi_sep,apsif,weight,dweight,ddweight)
    formfactors(:,i)=weight*formfactors(:,i)+1.d0-weight
  end do

end subroutine smear_formfactors
