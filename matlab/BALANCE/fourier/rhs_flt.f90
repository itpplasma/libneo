!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module rhs1_mod
    double precision :: dz_dphi,isw_rhs1
  end module rhs1_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs1(ndim,phi,y,dery)
!
  use rhs1_mod
!
  implicit none
!
  integer :: ndim ! = 5
!
  double precision :: phi,y,dery
  double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
  dimension y(ndim),dery(ndim)
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!
  x(1)=y(1)
  x(2)=phi
  x(3)=y(2)
!  
  call mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!
  dery(1)=hctrvr(1)/hctrvr(2)
  dery(2)=hctrvr(3)/hctrvr(2)
  dery(3)=y(1)*hctrvr(3)/hctrvr(2)
  if(isw_rhs1.eq.1) then
    dery(4)=y(1)
    dery(5)=y(2)
  elseif(isw_rhs1.eq.2) then
    dery(4)=bmod*y(1)*y(2)*hctrvr(1)
    dery(5)=y(1)**2*hctrvr(3)/hctrvr(2)
  endif
!
  dz_dphi=dery(2)
!
  return
  end subroutine rhs1
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module rhs2_mod
    integer :: mpol
    double precision :: aiota
  end module rhs2_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs2(ndim,phi,y,dery)
!
  use bdivfree_mod
  use rhs2_mod
  use field_eq_mod, only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
!
  implicit none
!
!   DBLE(y(1))                     - $R$
!   DBLE(y(1))                     - $Z$
!   y(2)-y(1+4*ntor*(2*mpol+1)) - subintegrands of Fourier transform
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: ind,n,m,ndim,ierr   !ndim = 2*ntor*(2*mpol+1)+1
!
  double precision :: r,phi,z,fac
  double precision :: x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
  double precision :: f,fr,fz,frr,frz,fzz
  double precision :: g,gr,gz,grr,grz,gzz
!
  dimension x(3),bder(3),hcovar(3),hctrvr(3),hcoder(3,3),hctder(3,3)
!
  double complex :: anr,anz,expon
  double complex, dimension(ndim) :: y,dery
!
!
  x(1)=dble(y(1))
  x(2)=phi
  x(3)=aimag(y(1))
!
  call mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
!
  ind=1
  dery(1)=dcmplx(hctrvr(1)/hctrvr(2),hctrvr(3)/hctrvr(2))
!
  fac=0.5d0*aiota/pi
  r=x(1)
  z=x(3)
!
  do n=1,ntor
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,arnim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anr=dcmplx(f,g)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznre(:,:,:,n),ipoint,r,z,   &
                f,fr,fz,frr,frz,fzz,ierr)
    call spline(nr,nz,rpoi,zpoi,hr,hz,icp,aznim(:,:,:,n),ipoint,r,z,   &
                g,gr,gz,grr,grz,gzz,ierr)
    anz=dcmplx(f,g)
    do m=-mpol,mpol
      expon=exp(dcmplx(0.d0,-m*aiota*phi))
! $A_R$ - component:
      ind=ind+1
      dery(ind)=expon*anr*fac
! $A_Z$ - component:
      ind=ind+1
      dery(ind)=expon*anz*fac
! $A^\psi = \bA\cdot\nabla\psi$ - component:
      ind=ind+1
      dery(ind)=expon*(anr*dpsidr+anz*dpsidz)*fac
! $A_\theta = \bA\cdot\difp{\br}{\theta}$ - component:
      ind=ind+1
      if(ind.gt.ndim) then
        print *,'rhs2 error : dimension exceeded'
        stop
      endif
      dery(ind)=expon*(anr*hctrvr(1)+anz*hctrvr(3))/(aiota*hctrvr(2))*fac
    enddo
  enddo
!
  return
  end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
