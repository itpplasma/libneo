!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE plag_coeff(npoi,nder,x,xp,coef)
    !
    ! npoi - number of points (determines the order of Lagrange
    ! polynomial
    ! which is equal npoi-1)
    ! nder - number of derivatives computed 0 - function only, 1 - first
    ! derivative
    ! x - actual point where function and derivatives are evaluated
    ! xp(npoi) - array of points where function is known
    ! coef(0:nder,npoi) - weights for computation of function and
    ! derivatives,
    ! f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
    ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value value
    !
    !
    INTEGER, INTENT(in)                                :: npoi,nder
    double precision, INTENT(in)                          :: x
    double precision, DIMENSION(npoi), INTENT(in)         :: xp
    double precision, DIMENSION(0:nder,npoi), INTENT(out) :: coef
    double precision, DIMENSION(:), ALLOCATABLE           :: dummy
    !
    INTEGER                                            :: i,k,j
    double precision                                      :: fac
    !
    DO i=1,npoi
       coef(0,i)=1.d0
       DO k=1,npoi
          IF(k.EQ.i) CYCLE
          coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
       ENDDO
    ENDDO
    !
    IF(nder.EQ.0) RETURN
    !
    ALLOCATE(dummy(npoi))
    !
    DO i=1,npoi
       dummy=1.d0
       dummy(i)=0.d0
       DO k=1,npoi
          IF(k.EQ.i) CYCLE
          fac=(x-xp(k))/(xp(i)-xp(k))
          DO j=1,npoi
             IF(j.EQ.k) THEN
                dummy(j)=dummy(j)/(xp(i)-xp(k))
             ELSE
                dummy(j)=dummy(j)*fac
             ENDIF
          ENDDO
       ENDDO
       coef(1,i)=SUM(dummy)
    ENDDO
    !
    DEALLOCATE(dummy)
    !
    RETURN
  END SUBROUTINE plag_coeff
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
!
