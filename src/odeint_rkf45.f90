module odeint_sub

implicit none

contains

subroutine odeint_allroutines(y, nvar, x1, x2, eps, derivs)
  use libneo_kinds, only : dp

  implicit none

  external :: derivs
  integer, intent(in) :: nvar
  real(dp), intent(in) :: x1, x2, eps
  real(dp), dimension(nvar) :: y
  real(dp), dimension(nvar) :: yp

  real(dp) :: epsrel, epsabs
  integer :: flag

  flag = 1
  epsrel = eps
  epsabs = 1d-31

  call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )

  if (any(flag.eq.[3,4,5,8])) print*, "the error flag in r8_rkf45 is: ", flag

  if (flag == 6) then
    epsrel = 10*epsrel
    epsabs = 10*epsabs
    flag = 2
    call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )
  elseif (flag == 7) then
    flag = 2
    call r8_rkf45 ( derivs, nvar, y, yp, x1, x2, epsrel, epsabs, flag )
  endif

end subroutine odeint_allroutines

end module odeint_sub
