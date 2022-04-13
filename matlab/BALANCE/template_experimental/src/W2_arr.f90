!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine W2_arr (x1_in,x2_in,Imn)

implicit none;

integer, parameter :: dp = 8
integer, parameter :: dpc = 8
real(dp), parameter :: pi    = 3.141592653589793238462643383279502884197_dp
real(dp) :: x1_in,x2_in;
complex(dpc), dimension(0:3,0:3) :: Imn
complex(dpc) :: t1, t2, F11m, x1, x2
!complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), one = cmplx(1.0d0, 0.0d0, 8);
double complex, parameter :: I = (0.0d0, 1.0d0), one = (1.0d0, 0.0d0)

integer :: l,nmax,m,n;

real(dp) :: F_im, F_re;
complex(dpc), allocatable, dimension(:,:,:) :: W2;

nmax=3

!allocate (W2(0:1, 0:3, -Nmax:Nmax))
allocate (W2(0:3, 0:3, 0:nmax))
W2 = (0.0d0, 0.0d0)

x1 = dcmplx(x1_in,0.d0)
t1 = x1**2

do l = 0, nmax ! does not work for 3, 2 and 1!


    x2 = x2_in + I*l

    t2 = - I*x2 + t1;

!      call Hypergeometric1F1_kummer_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);
!
!      F11m = F_re + I*F_im;

     call hypergeometric1f1_cont_fract_1_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);
!
     F11m = F_re + I*F_im;

!asymptotic form of expressions: has no fake singularity at kp=0

W2(0,0,l) =  &
(-I)*(-x2**2 + 2 + 5*x1**2 -  &
   (3*I)*x2*(1 + x1**2) + (3 + F11m)*x1**4)/ &
 ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
  (x2 + I*(2 + x1**2)))

W2(0,1,l) =  &
(-I)*x1*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2+ x1**4))/ &
 ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
  (x2 + I*(2 + x1**2)))

W2(0,2,l) =  &
-((x2 + I)*(2 + 3*x1**2 -  &
    I*x2*(1 - F11m*x1**2) + x1**4))/ &
  ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
   (x2 + I*(2 + x1**2)))

W2(0,3,l) =  &
(-I)*x1*(F11m*x2**3- (3 + 2*F11m)*x2 +  &
   I*x2**2*(3*F11m - x1**2) -  &
   I*(6 + (7 + 2*F11m)*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,0,l) =  &
(-I)*x1*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,1,l) =  &
(-I)*x2*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,2,l) =  &
(-I)*x1*(F11m*x2**3+ I*x2**2*(F11m - x1**2) -  &
   x2*(3 + 2*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,3,l) =  &
x2*((-I)*F11m*x2**3+ I*(3 + 2*F11m)*x2 -  &
   6 - (7 + 2*F11m)*x1**2 +  &
   x2**2*(3*F11m - x1**2) - x1**4)/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))


end do

W2(2,2,0) = I*x1*W2(1,2,1) + 2.d0*W2(1,1,1) - I*x1*W2(1,2,0) + W2(0,2,0)
W2(2,2,1) = I*x1*W2(1,2,2) + 2.d0*W2(1,1,2) - I*x1*W2(1,2,1) + W2(0,2,1)
W2(2,3,0) = I*x1*W2(2,2,1) + 2.d0*W2(1,2,1) - I*x1*W2(2,2,0) + 2.0d0*W2(1,2,0)
W2(2,3,1) = I*x1*W2(1,3,2) + 3.d0*W2(1,2,2) - I*x1*W2(1,3,1) + W2(0,3,1)
W2(3,3,0) = I*x1*W2(2,3,1) + 3.d0*W2(2,2,1) - I*x1*W2(2,3,0) + 2.0d0*W2(1,3,0)

Imn = W2(:,:,0)
do m=0,3
  do n=0,m-1
    Imn(m,n)=Imn(n,m)
  enddo
enddo

end subroutine

!------------------------------------------------------------------------------
