!
  module rk4d_mod
    integer :: nmax=0
    double precision, dimension(:), allocatable :: DYDX,YT,DYT,DYM
  end module rk4d_mod
!
  module rk4c_mod
    integer :: nmax=0
    double complex,   dimension(:), allocatable :: DYDX,YT,DYT,DYM
  end module rk4c_mod
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE RK4D(Y,N,X,H,DERIVS)
!
      use rk4d_mod
!
      DOUBLE PRECISION Y,X,H,HH,H6,XH
      EXTERNAL DERIVS
      DIMENSION Y(N)
!
      if(n.ne.nmax) then
        if(allocated(dydx)) deallocate(dydx)
        if(allocated(yt)) deallocate(yt)
        if(allocated(dyt)) deallocate(dyt)
        if(allocated(dym)) deallocate(dym)
        nmax=n
        allocate(DYDX(NMAX),YT(NMAX),DYT(NMAX),DYM(NMAX))
      endif
!
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      CALL DERIVS(N,X,Y,DYDX)
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(N,XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(N,XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(N,X+H,YT,DYT)
      DO 14 I=1,N
        Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
14    CONTINUE
        X=X+H
      RETURN
      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      SUBROUTINE RK4C(Y,N,X,H,DERIVS)
!
      use rk4c_mod
!
      DOUBLE PRECISION X,H,HH,H6,XH
      EXTERNAL DERIVS
      double complex, dimension(n) :: Y
!
      if(n.ne.nmax) then
        if(allocated(dydx)) deallocate(dydx)
        if(allocated(yt)) deallocate(yt)
        if(allocated(dyt)) deallocate(dyt)
        if(allocated(dym)) deallocate(dym)
        nmax=n
        allocate(DYDX(NMAX),YT(NMAX),DYT(NMAX),DYM(NMAX))
      endif
!
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      CALL DERIVS(N,X,Y,DYDX)
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(N,XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(N,XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(N,X+H,YT,DYT)
      DO 14 I=1,N
        Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
14    CONTINUE
        X=X+H
      RETURN
      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

