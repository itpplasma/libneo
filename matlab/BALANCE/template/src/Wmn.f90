module params

  integer :: nsize, lmax, jmax, jpmax, mnmax
  double complex :: bx1, bx20, gam, bx4, theta

end module params

module help_igamma_mod
  double complex :: zbeg,zend
end module help_igamma_mod

!-----------------------------------------------------------------------------!

subroutine incomplete_gamma(s, z, abs_err, rel_err, G, est_err, Nterms, status)

! Evaluates G(s, z) by continued fractions method.
! abs_err, rel_err - desired tolerances
! G - G(s, z) value
! est_err - estimated relative error
! Nterms - number of terms used
! status == 0 - Ok.

implicit none;

complex(8), intent(in)   :: s, z;
real(8),    intent(in)   :: abs_err, rel_err;
complex(8), intent(out)  :: G;
real(8),    intent(out)  :: est_err;
integer,    intent(out)  :: Nterms;
integer,    intent(out)  :: status;

integer,       parameter :: prec = 8;
integer,       parameter :: Nmin = 2, Nmax = 4194304;
!complex(prec), parameter :: NaN = cmplx(1.0d100, 1.0d100, prec);
!complex(prec), parameter :: c0 = cmplx(0.0d0, 0.0d0, prec), c1 = cmplx(1.0d0, 0.0d0, prec);
double complex, parameter :: NaN = (1.0d100, 1.0d100);
double complex, parameter :: c0 = (0.0d0, 0.0d0), c1 = (1.0d0, 0.0d0);

complex(prec), dimension(Nmax) :: a, b;

complex(prec) :: term, F, Fp;

integer :: k, j, M, N;

logical :: goon;

M = 0;
N = 4;

F  = NaN;
Fp = NaN;

goon = .true.;

do while (N <= Nmax .and. goon)

    term = c0; ! tail value

    do k = N,Nmin,-1

        if (k > M) then

            j = k / 2;

            a(k) = cmplx(j, 0, prec);

            if (j * 2 .eq. k) a(k) = a(k) - s;

            b(k) = c1;

        end if

        term = a(k) / z / (b(k) + term);

    end do

    Fp = F;

    F = c1 / z / (c1 + term); ! top level term

    goon = ( abs(F - Fp) > abs_err + rel_err * abs(F) );

    if (N > M) M = N;

    N = N * 2;

end do

est_err = abs(c1 - Fp/F);

if (goon) then
    status = 1;
else
    status = 0;
end if

Nterms = N / 2;

!G = z**s * exp(-z) * F;
G = z**s * F;

end subroutine

!-----------------------------------------------------------------------------!

!--- subroutines ---------------------------------------------------------------

subroutine dwdtau(x, y, f)
   
   use params, only: bx1, bx20, gam, bx4, nsize, lmax, jmax, jpmax, mnmax, theta
   implicit none
   double complex, parameter :: imun=(0.d0,1.d0)
   integer :: j, jp, l, mn, k, jjpmax
   double precision :: x
   double complex   :: tau, t
   double precision, dimension(nsize) :: y, f
   double complex :: expon, subint, ax4l, argexp, exponpow, exponax4, exponj
   double complex :: aj, ajp, ajmn

   tau=x*theta
   t=exp(-tau)
!
   ajmn=(1.d0-t)/(1.d0+gam*(1.d0+t))
   aj  =(1.d0+t)/(1.d0+gam*(1.d0+t))
   ajp =(1.d0-t)/(1.d0+gam*(1.d0-t))
   ax4l=1.d0/sqrt(1.d0+imun*bx4*tau)
   argexp=(imun*bx20-bx1**2)*tau+bx1**2*(1.d0+2.d0*gam)*ajmn
   aj  =aj/ajmn**2
   ajp =ajp/ajmn**2
!
   expon = exp(argexp)/sqrt((1.d0+gam)**2-(gam*t)**2)*ax4l
   expon = expon*theta
!
   k=0
   do mn=0,mnmax
     jjpmax=mn/2
     exponax4=expon
     do l=1,lmax
       expon=expon*ax4l
       exponj=expon
       do j=0,jjpmax
         exponpow=exponj
         do jp=0,jjpmax-j
           k=k+1
           f(k)=real(exponpow)
           k=k+1
           f(k)=dimag(exponpow)
           exponpow=exponpow*ajp
         end do
         exponj=exponj*aj
       end do
     end do 
     expon=exponax4*ajmn
   end do  

end subroutine dwdtau

subroutine help_igamma(x, y, f)

  use params, only: lmax
  use help_igamma_mod, only : zbeg,zend
  implicit none
  integer :: l, k
  double precision :: x
  double complex   :: z,dzdx
  double precision, dimension(2*lmax) :: y, f
  double complex :: expon, sqinv
  double complex :: aj, ajp, ajmn
!
  dzdx=zend-zbeg
  z=zbeg+dzdx*x
  sqinv=(1.d0,0.d0)/sqrt(z)
  expon=exp(zbeg-z)*sqinv*dzdx
  k=0
  do l=1,lmax
    expon=expon*sqinv
    k=k+1
    f(k)=real(expon)
    k=k+1
    f(k)=dimag(expon)
  enddo

end subroutine help_igamma

subroutine evaluate_integral(wintegral)
   
  use params
  use help_igamma_mod, only : zbeg,zend
  implicit none
  double precision, parameter :: eps = 1.0d-12, taumax=30.d0, zshift=0.1d0
  logical :: switch
  integer :: k, mn, l, j, jp, jjpmax, Nterms, ierr
  double precision, dimension(:), allocatable :: y
  double precision :: x_ini, x_end, est_err 
  double complex :: argexp,expon,z,sqxbfac,sqix4fac,oneovopg,znumer,zdenom
  double complex :: G,s,sqonepix4t,Gscale
  double complex, dimension(0:mnmax,lmax,0:jmax,0:jpmax) :: wintegral
  double complex, dimension(:), allocatable :: delG
!
  external :: dwdtau,help_igamma
!    
  k=0
  do mn=0,mnmax
    jjpmax=mn/2
    do l=1,lmax
      do j=0,jjpmax
        do jp=0,jjpmax-j
          k=k+2
        end do
      end do
    end do
  end do
!
  nsize=k
  allocate(y(nsize))
!  
!
  x_ini=1.d-14
  x_end = 20.0d0/abs(bx1) + 400.0d0/abs(bx1)**2
  switch=x_end.gt.taumax
  if(switch) x_end=taumax
  y = 0.0d0 
!
  call odeint_allroutines(y, nsize, x_ini, x_end, eps, dwdtau)
!
  k=0
  do mn=0,mnmax
    jjpmax=mn/2
    do l=1,lmax
      do j=0,jjpmax
        do jp=0,jjpmax-j
          k=k+1
          wintegral(mn,l,j,jp)=dcmplx(y(k),y(k+1))
          k=k+1
        end do
      end do
    end do
  end do
   
  deallocate(y)
!
  if(switch) then
    oneovopg=1.d0/(1.d0+gam)
    sqxbfac=bx1**2-(0.d0,1.d0)*bx20
    argexp=-sqxbfac*taumax*theta+bx1**2*(1.d0+2.d0*gam)/(1.d0+gam)
    expon=exp(argexp)
    sqonepix4t=1.d0/sqrt(1.d0+(0.d0,1.d0)*bx4*taumax*theta)
    znumer=sqxbfac*(1.d0+(0.d0,1.d0)*bx4*taumax*theta)
    sqxbfac=sqrt(sqxbfac)
    zdenom=(0.d0,1.d0)*bx4
    if(abs(zdenom).gt.abs(znumer)*eps) then
      allocate(delG(lmax))
      z=znumer/zdenom
      sqix4fac=1.d0/sqrt(zdenom)
      if(real(z).lt.0.d0.and.abs(dimag(z)).lt.zshift) then
        zbeg=z
        zend=z+dcmplx(0.d0,sign(zshift,dimag(z)))
        nsize=k
        allocate(y(nsize))
        x_ini=0.d0
        x_end=1.d0
        y=0.d0
!
        call odeint_allroutines(y, nsize, x_ini, x_end, eps, help_igamma)
!
        k=0
        do l=1,lmax
          delG(l)=dcmplx(y(k+1),y(k+2))
          k=k+2
        enddo
        deallocate(y)
        z=zend
        Gscale=exp(zbeg-zend)
      else
        delG=(0.d0,0.d0)
        Gscale=(1.d0,0.d0)
      endif
!
      do l=1,lmax
        s=1.d0-0.5d0*dfloat(l+1)
!
        call incomplete_gamma(s,z,eps,eps,G,est_err,Nterms,ierr)
!
        G=G*Gscale+delG(l)
        do mn=0,mnmax
          jjpmax=mn/2
          do j=0,jjpmax
            do jp=0,jjpmax-j
              wintegral(mn,l,j,jp)=wintegral(mn,l,j,jp)+expon*G         &
                                  *oneovopg**(mn+1-j-jp)*sqxbfac**(l-1) &
                                  *sqix4fac**(l+1)
            end do
          end do
        end do
      end do
!
      deallocate(delG)
    else
!
      do l=1,lmax
        do mn=0,mnmax
          jjpmax=mn/2
          do j=0,jjpmax
            do jp=0,jjpmax-j
              wintegral(mn,l,j,jp)=wintegral(mn,l,j,jp)+expon         &
                                  *oneovopg**(mn+1-j-jp)/znumer       &
                                  *sqonepix4t**(l-1)
            end do
          end do
        end do
      end do
!
    endif
  endif
!   
end subroutine evaluate_integral

subroutine GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)
 
  use params

  implicit none
  double precision, parameter :: addimpart=0.1d0,epscontour=0.03d0
  integer :: mpar, lperp, m, k, kp, n, j, jp, l, lp
  double precision :: x1,x2,x3,x4
  double complex :: onemin4ix3, onemin4ix3pow, bx1ovi, omix3, sqomix3
  double complex :: oldarg
  double complex, dimension(0:mpar,0:mpar,1:lperp) :: Imnl

  double precision, dimension(:),       allocatable :: doubfac
  double precision, dimension(:,:),     allocatable :: bincoef
  double precision, dimension(:,:,:,:), allocatable :: symbm
  double complex,   dimension(:,:,:,:), allocatable :: calsymbm
  double complex,   dimension(:,:,:,:), allocatable :: wintegral

  mnmax=2*mpar
  lmax=lperp
  jmax=mpar
  jpmax=mpar 

  allocate(doubfac(0:mpar))
  allocate(bincoef(0:mpar,0:mpar))
  allocate(symbm(0:mpar,0:mpar,0:mpar,0:mpar))
  allocate(calsymbm(0:mpar,0:mpar,0:mpar,0:mpar))
  allocate(wintegral(0:mnmax,lmax,0:jmax,0:jpmax))
!
  doubfac(0)=1.d0
  do m=1,mpar
    doubfac(m)=doubfac(m-1)*dfloat(2*m-1)
  enddo 
!
  bincoef=0.d0
  do m=0,mpar
    bincoef(m,0)=1.d0
    do k=1,m
      bincoef(m,k)=bincoef(m,k-1)*dfloat(m-k+1)/dfloat(k)
    enddo
  enddo
!
  symbm=0.d0
!
  do m=0,mpar
    do n=0,mpar
      do j=0,mpar
        do jp=0,mpar
          symbm(m,n,j,jp)=0.d0
          do k=0,m
            do kp=0,n
              do l=0,k
                do lp=0,kp
                  if (2*j.eq.l+lp.and.2*jp.eq.k+kp-l-lp) then
                    symbm(m,n,j,jp)=symbm(m,n,j,jp)+(-1.d0)**(k-l) &
                                   *bincoef(m,k)*bincoef(n,kp)     &
                                   *bincoef(k,l)*bincoef(kp,lp)
                  endif 
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
! 
  omix3=1.d0-(0.d0,4.d0)*x3
  sqomix3=sqrt(omix3)
  gam=(0.d0,2.d0)*x3/(omix3+sqomix3)
  onemin4ix3=1.d0/sqrt(sqomix3)
  bx1=x1/sqomix3*onemin4ix3
  bx1ovi=(0.d0,-1.d0)*bx1
  bx20=(x2+2.d0*x3/(sqomix3+1.d0))/sqomix3 
  bx4=x4/sqomix3 
  oldarg=(0.d0,1.d0)*bx20-bx1**2
  theta=dcmplx(1.d0,sign(addimpart                                            &
       +max(0.d0,real(oldarg)/max(abs(dimag(oldarg)),abs(oldarg)*epscontour)) &
       ,dimag(oldarg)))
!
  do m=0,mpar
    do n=0,mpar
      onemin4ix3pow=onemin4ix3**(m+n+3)
      do j=0,mpar
        do jp=0,mpar
          calsymbm(m,n,j,jp)=symbm(m,n,j,jp)*doubfac(j)*doubfac(jp) &
                            *0.5d0**(j+jp)*bx1ovi**(m+n-2*(j+jp))   &
                            *onemin4ix3pow 
        enddo
      enddo
    enddo
  enddo
!
  call  evaluate_integral(wintegral)
!  
  Imnl=(0.d0,0.d0)
  do m=0,mpar
    do n=0,mpar
      do l=1,lperp
        do j=0,mpar
          do jp=0,mpar
            Imnl(m,n,l)=Imnl(m,n,l)+wintegral(m+n,l,j,jp)*calsymbm(m,n,j,jp)
          enddo
        enddo
      enddo
    enddo
  enddo
  
 
end subroutine GreenMom

subroutine getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,Imnkl)
!
  implicit none
!
  logical :: consenergy
  integer :: mpar,mperp,lperp,l,m,n,k
  double precision :: x1, x2, x3, x4, arggam
  double complex, dimension(0:mpar,0:mpar,0:mperp,0:mperp)  :: Imnkl
  double precision, dimension(:),     allocatable :: gamkl
  double complex,   dimension(:,:,:), allocatable :: Imnl
!
  lperp=2*mperp+1
  allocate(Imnl(0:mpar,0:mpar,1:lperp))
  if(consenergy) then
    allocate(gamkl(0:2*mperp))
    do k=0,2*mperp
      arggam=0.5d0*dfloat(k)+1
      gamkl(k)=gamma(arggam)
    enddo
  endif
!
  call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl) 
!
  do k=0,mperp
    do l=0,mperp
      do m=0,mpar
        do n=0,mpar
          if(consenergy) then
            Imnkl(m,n,k,l)=Imnl(m,n,k+l+1)+gamkl(k)*gamkl(l)/gamkl(k+l)      &
                *(Imnl(m,2,k+1)-Imnl(m,0,k+1))*(Imnl(n,2,l+1)-Imnl(n,0,l+1)) &
                /(1.d0-Imnl(0,0,1)+2.d0*Imnl(2,0,1)-Imnl(2,2,1))
          else
            Imnkl(m,n,k,l)=Imnl(m,n,k+l+1)
          endif
        enddo
      enddo
    enddo
  enddo
!
  deallocate(Imnl)
  if(consenergy) deallocate(gamkl)
end subroutine getIfunc_drift
!--------------------------------------------

subroutine test

   use params
   implicit none
   logical :: consenergy
   integer :: i, nmax, mpar, mperp, lperp, m, j, jp, l, n
   double precision :: hx, x1, x2, x3, x4
   double complex, parameter :: imun=(0.d0,1.d0)
   double complex :: z1,z2
   double complex, dimension(0:3,0:3) :: Imn
   double complex, dimension(:,:,:), allocatable  :: Imnl
   double complex, dimension(:,:,:,:), allocatable  :: Imnkl
   double precision, dimension(:), allocatable :: doubfac
!
   mpar=3
   mperp=2
   lperp=2*mperp+1
   allocate(Imnl(0:mpar,0:mpar,1:lperp))
   allocate(Imnkl(0:mpar,0:mpar,0:mperp,0:mperp))
   allocate(doubfac(-2:2*mpar-1))
!
   doubfac=0.d0
   doubfac(-1)=1.d0
   do m=1,mpar
     doubfac(2*m-1)=doubfac(2*m-3)*dfloat(2*m-1)
   enddo

   x2=1.0d0
   x3=0.0d0
   x4=0.01d0
   nmax = 100
   hx   = 1.1d1/nmax
   do i=-nmax,nmax
print *,i
     if(i.eq.0) cycle
     x1=hx*i
!
     call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)
!
     write(100,*) x1, real(Imnl(0,0,1)), dimag(Imnl(0,0,1))
     write(101,*) x1, real(Imnl(0,1,1)), dimag(Imnl(0,1,1))
     write(102,*) x1, real(Imnl(0,2,1)), dimag(Imnl(0,2,1))
     write(103,*) x1, real(Imnl(0,3,1)), dimag(Imnl(0,3,1))
     write(110,*) x1, real(Imnl(1,0,1)), dimag(Imnl(1,0,1))
     write(111,*) x1, real(Imnl(1,1,1)), dimag(Imnl(1,1,1))
     write(112,*) x1, real(Imnl(1,2,1)), dimag(Imnl(1,2,1))
     write(113,*) x1, real(Imnl(1,3,1)), dimag(Imnl(1,3,1))
     write(120,*) x1, real(Imnl(2,0,1)), dimag(Imnl(2,0,1))
     write(121,*) x1, real(Imnl(2,1,1)), dimag(Imnl(2,1,1))
     write(122,*) x1, real(Imnl(2,2,1)), dimag(Imnl(2,2,1))
     write(123,*) x1, real(Imnl(2,3,1)), dimag(Imnl(2,3,1))
     write(130,*) x1, real(Imnl(3,0,1)), dimag(Imnl(3,0,1))
     write(131,*) x1, real(Imnl(3,1,1)), dimag(Imnl(3,1,1))
     write(132,*) x1, real(Imnl(3,2,1)), dimag(Imnl(3,2,1))
     write(133,*) x1, real(Imnl(3,3,1)), dimag(Imnl(3,3,1))
!
     consenergy=.true.
!
     call getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,Imnkl)
!
     z1=Imnkl(2,1,0,0)
     z2=Imnkl(3,0,0,0)
     write(1001,*) x1,abs(z1-z2)/(abs(z1)+abs(z2)),real(z1-z2),dimag(z1-z2)
!
     write(200,*) x1, real(Imnkl(0,0,0,0)), dimag(Imnkl(0,0,0,0))
     write(201,*) x1, real(Imnkl(0,1,0,0)), dimag(Imnkl(0,1,0,0))
     write(202,*) x1, real(Imnkl(0,2,0,0)), dimag(Imnkl(0,2,0,0))
     write(203,*) x1, real(Imnkl(0,3,0,0)), dimag(Imnkl(0,3,0,0))
     write(210,*) x1, real(Imnkl(1,0,0,0)), dimag(Imnkl(1,0,0,0))
     write(211,*) x1, real(Imnkl(1,1,0,0)), dimag(Imnkl(1,1,0,0))
     write(212,*) x1, real(Imnkl(1,2,0,0)), dimag(Imnkl(1,2,0,0))
     write(213,*) x1, real(Imnkl(1,3,0,0)), dimag(Imnkl(1,3,0,0))
     write(220,*) x1, real(Imnkl(2,0,0,0)), dimag(Imnkl(2,0,0,0))
     write(221,*) x1, real(Imnkl(2,1,0,0)), dimag(Imnkl(2,1,0,0))
     write(222,*) x1, real(Imnkl(2,2,0,0)), dimag(Imnkl(2,2,0,0))
     write(223,*) x1, real(Imnkl(2,3,0,0)), dimag(Imnkl(2,3,0,0))
     write(230,*) x1, real(Imnkl(3,0,0,0)), dimag(Imnkl(3,0,0,0))
     write(231,*) x1, real(Imnkl(3,1,0,0)), dimag(Imnkl(3,1,0,0))
     write(232,*) x1, real(Imnkl(3,2,0,0)), dimag(Imnkl(3,2,0,0))
     write(233,*) x1, real(Imnkl(3,3,0,0)), dimag(Imnkl(3,3,0,0))
   end do


end subroutine test
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module get_matrix_tmp
    integer :: ispec
  end module get_matrix_tmp
!
!-----------------------------------------------------------------------------
!
  subroutine calc_transport_coeffs_ornuhl_drift(isp,ndim,D_11,D_12,D_21,D_22)
!
  use baseparam_mod, only     : Z_i,e_charge,am,p_mass,c,btor,e_mass
  use grid_mod, only          : npoib,rb,params_b
  use wave_code_data, only    : nui,nue,Es,Br,B0
  use get_matrix_tmp, only    : ispec
  use sample_matrix_mod, only : nlagr,itermax,xbeg,xend,eps,npoi,xarr,amat_arr
!
  implicit none
!
  integer, parameter :: mpar=3,mperp=2
  integer :: isp,ndim,i,np,ierr
  double precision :: rloc,anu,v_T,comfac,epm2,brm2,epbr_re,epbr_im,D_12a
  double precision, dimension(ndim) :: D_11,D_12,D_21,D_22
  double complex, dimension(0:mpar,0:mpar,0:mperp,0:mperp) :: Imnkl
!
  nlagr=3
  eps=1.d-4
  itermax=30
!
  ispec=isp
!
  xbeg=rb(1)
  xend=rb(npoib)
  if(ispec.eq.1) then
    print *,'start sample_matrix_bal, electrons'
open(8888,file='Imnkl0000_e.dat')
  else
    print *,'start sample_matrix_bal, ions'
open(8888,file='Imnkl0000_i.dat')
  endif
!
  call sample_matrix_bal(ierr)
!
  print *,'sample_matrix_bal error: ',ierr,'  number of points: ',npoi
!
  do i=1,npoib
!
    call interpolateIfunc_drift(rb(i),Imnkl)
!
write(8888,*) rb(i),real(Imnkl(0,0,0,0)),dimag(Imnkl(0,0,0,0))
    if(ispec.eq.1) then
!electrons:
      anu=nue(i)
      v_T=sqrt(params_b(3,i)/e_mass)
    else
!ions:
      anu=nui(i)
      v_T=sqrt(params_b(4,i)/p_mass/am)
    endif
    comfac=0.5d0/(anu*B0(i)**2)
    epm2=c**2*abs(Es(i))**2
    brm2=v_T**2*abs(Br(i))**2
    epbr_re=2.d0*c*v_T*real(conjg(Es(i))*Br(i))
    epbr_im=2.d0*c*v_T*dimag(conjg(Es(i))*Br(i))
!
    D_11(i)=comfac*(epm2   *real(Imnkl(0,0,0,0))                              &
           +        epbr_re*real(Imnkl(1,0,0,0))                              &
           +        brm2   *real(Imnkl(1,1,0,0)))
    D_12(i)=comfac*(epm2   *real(Imnkl(0,0,0,2)+0.5d0*Imnkl(2,0,0,0))         &
           +        epbr_re*real(Imnkl(1,0,0,2)                               &
           +                     0.25d0*(Imnkl(3,0,0,0)+Imnkl(2,1,0,0)))      &
           +        brm2   *real(Imnkl(1,1,0,2)+0.5d0*Imnkl(3,1,0,0)))
    D_21(i)=D_12(i)
    D_22(i)=comfac*(epm2   *real(2.d0*Imnkl(0,0,2,2)+Imnkl(2,0,0,2)           &
           +                     0.25d0*Imnkl(2,2,0,0))                       &
           +        epbr_re*real(2.d0*Imnkl(1,0,2,2)                          &
           +                     0.5d0*(Imnkl(3,0,0,2)+Imnkl(2,1,0,2))        &
           +                     0.25d0*Imnkl(3,2,0,0))                       &
           +        brm2   *real(2.d0*Imnkl(1,1,2,2)+Imnkl(3,1,0,2)           &
           +                     0.25d0*Imnkl(3,3,0,0)))
!
    D_12a=comfac*epbr_im*0.25d0*dimag(Imnkl(2,1,0,0)-Imnkl(3,0,0,0))
!
    D_12(i) = D_12(i)+D_12a
    D_21(i) = D_21(i)-D_12a
!
  enddo
close(8888)
!
  deallocate(xarr,amat_arr)
!
  end subroutine calc_transport_coeffs_ornuhl_drift
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_matrix
!
  use baseparam_mod, only     : Z_i,e_charge,am,p_mass,c,e_mass,ev
  use grid_mod, only          : npoib,rb,params_b
  use wave_code_data, only    : kp,ks,om_E,B0t,B0,nui,nue
  use get_matrix_tmp, only    : ispec
  use sample_matrix_mod, only : n1,n2,x,amat
!
  implicit none
!
  integer, parameter :: npoi_der=4,nder=1
  integer, parameter :: mpar=3,mperp=2
!
  logical :: consenergy
  integer :: i,ipb,ipe
  double precision :: r,bprime,bmod,bthe,akpar,akperp,omega_E,T,anu,v_T,om_c
  double precision :: x1,x2,x3,x4
  double precision, dimension(0:nder,npoi_der) :: coef
  double complex, dimension(0:mpar,0:mpar,0:mperp,0:mperp) :: Imnkl
!
  r=x
!
  call binsrc(rb,1,npoib,r,i)
!
  ipb=i-npoi_der/2
  ipe=ipb+npoi_der-1
  if(ipb.lt.1) then
    ipb=1
    ipe=ipb+npoi_der-1
  elseif(ipe.gt.npoib) then
    ipe=npoib
    ipb=ipe-npoi_der+1
  endif
!
  call plag_coeff(npoi_der,nder,r,rb(ipb:ipe),coef)
!
  bprime=sum(coef(1,:)*B0(ipb:ipe))
  bmod=sum(coef(0,:)*B0(ipb:ipe))
  bthe=sum(coef(0,:)*B0t(ipb:ipe))
  akpar=sum(coef(0,:)*kp(ipb:ipe))
  akperp=sum(coef(0,:)*ks(ipb:ipe))
  omega_E=sum(coef(0,:)*om_E(ipb:ipe))
  if(ispec.eq.1) then
!electrons:
    T=sum(coef(0,:)*params_b(3,ipb:ipe))
    anu=sum(coef(0,:)*nue(ipb:ipe))
    v_T=sqrt(T/e_mass)
    om_c=-e_charge*bmod/(e_mass*c)
  else
!ions:
    T=sum(coef(0,:)*params_b(4,ipb:ipe))
    anu=sum(coef(0,:)*nui(ipb:ipe))
    v_T=sqrt(T/p_mass/am)
    om_c=Z_i*e_charge*bmod/(am*p_mass*c)
  endif
  x1=akpar*v_T/anu
  x2=-omega_E/anu
  x3=akperp*(bthe/bmod)**2*v_T**2/(r*anu*om_c)
  x4=akperp*bprime*v_T**2/(anu*bmod*om_c)
!
  consenergy=.true.
!
  call getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,Imnkl)
!
  if(allocated(amat)) deallocate(amat)
  n1=2
  n2=10
  allocate(amat(n1,n2))
!
  amat(1,1)=Imnkl(0,0,0,0)
  amat(1,2)=Imnkl(0,1,0,0)
  amat(1,3)=Imnkl(0,2,0,0)
  amat(1,4)=Imnkl(0,3,0,0)
!
  amat(1,5)=Imnkl(1,1,0,0)
  amat(1,6)=Imnkl(1,2,0,0)
  amat(1,7)=Imnkl(1,3,0,0)
!
  amat(1,8)=Imnkl(2,2,0,0)
  amat(1,9)=Imnkl(2,3,0,0)
!
  amat(1,10)=Imnkl(3,3,0,0)
!
  amat(2,1)=Imnkl(0,0,0,2)
  amat(2,2)=Imnkl(0,1,0,2)
  amat(2,3)=Imnkl(0,2,0,2)
  amat(2,4)=Imnkl(0,3,0,2)
!
  amat(2,5)=Imnkl(1,1,0,2)
  amat(2,6)=Imnkl(1,2,0,2)
  amat(2,7)=Imnkl(1,3,0,2)
!
  amat(2,8)=Imnkl(0,0,2,2)
  amat(2,9)=Imnkl(0,1,2,2)
!
  amat(2,10)=Imnkl(1,1,2,2)
!
  end subroutine get_matrix
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine interpolateIfunc_drift(r,Imnkl)
!
  use sample_matrix_mod, only : nlagr,n1,n2,npoi,xarr,amat_arr,amat
!
  integer, parameter :: mpar=3,mperp=2,nder=0
!
  integer :: i,j,npoilag,nshift,ibeg,iend
  double precision :: r
  double complex, dimension(0:mpar,0:mpar,0:mperp,0:mperp) :: Imnkl
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef
!
  npoilag=nlagr+1
  nshift=npoilag/2
  ALLOCATE(coef(0:nder,npoilag))
!
  call binsrc(xarr,1,npoi,r,i)
!
  ibeg=i-nshift
  iend=ibeg+nlagr
  if(ibeg.lt.1) then
    ibeg=1
    iend=ibeg+nlagr
  elseif(iend.gt.npoi) then
    iend=npoi
    ibeg=iend-nlagr
  endif
!
  CALL plag_coeff(npoilag,nder,r,xarr(ibeg:iend),coef)
!
  do i=1,n1
    do j=1,n2
      amat(i,j)=sum(coef(0,:)*amat_arr(i,j,ibeg:iend))
    enddo
  enddo
!
  deallocate(coef)
!
  Imnkl(0,0,0,0)=amat(1,1)
  Imnkl(0,1,0,0)=amat(1,2)
  Imnkl(0,2,0,0)=amat(1,3)
  Imnkl(0,3,0,0)=amat(1,4)
!
  Imnkl(1,1,0,0)=amat(1,5)
  Imnkl(1,2,0,0)=amat(1,6)
  Imnkl(1,3,0,0)=amat(1,7)
!
  Imnkl(2,2,0,0)=amat(1,8)
  Imnkl(2,3,0,0)=amat(1,9)
!
  Imnkl(3,3,0,0)=amat(1,10)
!
  Imnkl(0,0,0,2)=amat(2,1)
  Imnkl(0,1,0,2)=amat(2,2)
  Imnkl(0,2,0,2)=amat(2,3)
  Imnkl(0,3,0,2)=amat(2,4)
!
  Imnkl(1,1,0,2)=amat(2,5)
  Imnkl(1,2,0,2)=amat(2,6)
  Imnkl(1,3,0,2)=amat(2,7)
!
  Imnkl(0,0,2,2)=amat(2,8)
  Imnkl(0,1,2,2)=amat(2,9)
!
  Imnkl(1,1,2,2)=amat(2,10)
!
  Imnkl(1,0,0,0:2:2)=Imnkl(0,1,0,0:2:2)
  Imnkl(2,0,0,0:2:2)=Imnkl(0,2,0,0:2:2)
  Imnkl(3,0,0,0:2:2)=Imnkl(0,3,0,0:2:2)
!
  Imnkl(2,1,0,0:2:2)=Imnkl(1,2,0,0:2:2)
  Imnkl(3,1,0,0:2:2)=Imnkl(1,3,0,0:2:2)
!
  Imnkl(3,2,0,0)=Imnkl(2,3,0,0)
!
  end subroutine interpolateIfunc_drift
