module plag_coeff_sub
  implicit none

  contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE plag_coeff(npoi,nder,x,xp,coef)
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
    ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value

    use libneo_kinds, only : dp

    INTEGER, INTENT(in) :: npoi,nder
    real(dp), INTENT(in) :: x
    real(dp), DIMENSION(npoi), INTENT(in) :: xp
    real(dp), DIMENSION(0:nder,npoi), INTENT(out) :: coef
    real(dp), DIMENSION(:), ALLOCATABLE :: dummy

    INTEGER :: i,k,j
    real(dp) :: fac

    DO i=1,npoi
      coef(0,i) = 1.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        coef(0,i) = coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
      END DO
    END DO

    IF (nder.EQ.0) RETURN

    ALLOCATE(dummy(npoi))

    DO i=1,npoi
      dummy = 1.d0
      dummy(i) = 0.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        fac = (x-xp(k))/(xp(i)-xp(k))
        DO j=1,npoi
          IF (j.EQ.k) THEN
            dummy(j) = dummy(j)/(xp(i)-xp(k))
          ELSE
            dummy(j) = dummy(j)*fac
          END IF
        END DO
      END DO
      coef(1,i) = SUM(dummy)
    END DO

    DEALLOCATE(dummy)

  END SUBROUTINE plag_coeff

end module plag_coeff_sub
