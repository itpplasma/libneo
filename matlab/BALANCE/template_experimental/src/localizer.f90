!
  subroutine localizer(sig,x1,x2,x,weight)
!
  implicit none
!
  double precision :: sig,x1,x2,x,t,weight
!
  if(sig.gt.0.d0) then
! from 1 to 0:
    t=(x-x1)/(x2-x1)
  else
! from 0 to 1:
    t=(x2-x)/(x2-x1)
  endif
!
  if(t.le.0.d0) then
    weight=1.d0
  elseif(t.ge.1.d0) then
    weight=0.d0
  else
    weight=exp(-6.283185307179586d0/(1.d0-t)*exp(-1.414213562373095d0/t))
  endif
!
  return
  end subroutine localizer
