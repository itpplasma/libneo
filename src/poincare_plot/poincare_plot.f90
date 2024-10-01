module neo_poincare_plot
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
    

contains


subroutine make_poincare_plot(field, rmn, rmx, zmn, zmx, idir, nplot, rmnplot, &
                              rmxplot, phi0, zet0, ah, hn1, nphi)
    use ode_integration, only: odeint
    use neo_field_base, only: field_t

    class(field_t), intent(in) :: field

    real(dp), parameter :: pi = 3.14159265358979d0
    real(dp), parameter :: per_phi = 2*pi
    real(dp), parameter :: neqn = 2
    integer :: idir,nplot,nphi,iunit,npoint,i,j,nok,nbad
    
    real(dp) :: y(2)
    real(dp) :: rmn,rmx,zmn,zmx,raxis,zaxis,rmnplot,rmxplot,phi0,zet0
    real(dp) :: ah,hn1,hn1min,relerr,phiout,divB,epsB
    real(dp) :: phi,rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

    phi0 = phi0*2.*pi
    ah = (rmxplot - rmnplot)/(nplot-1)
    if(nplot .gt. 50) STOP 'too large number of plots'

    ! parameters for ODEINT
    hn1 = 0.1*idir
    hn1min = 0.
    !  relerr = 1.e-12
    relerr = 1.e-8
    !  nphi = 1000
    nphi = 3000
    !  nphi = 10000
    do i=1, nplot
        iunit = 50 + (i-1)*idir
        npoint = 0
        phi = phi0

        y(1) = rmnplot + (i-1)*ah
        y(2) = zet0
        write(*,*)'surface #',i
        do j=1,nphi
    print *,j,nphi

            phiout = phi + per_phi*idir

            call odeint(y,neqn,phi,phiout,relerr,hn1,hn1min,nok,nbad,fieldline_derivative) 

            npoint = npoint + 1
            phi = phiout
            rrr = y(1)  ! cylindic R
            ppp = phi   ! cylindic $\phi$
            zzz = y(2)  ! cylindic Z
            ppp = set_to_period(ppp)
            call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ                &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

            divB = Br/rrr + dBrdR + dBzdZ + dBpdp/rrr
            epsB = divB*rrr/Bp
            write(iunit,*)rrr,zzz,phi !,epsB,Br,Bp,Bz,rad_old,theta_old

            if( y(1).lt.rmn .or. y(1).gt.rmx .or.            &
                        y(2).lt.zmn .or. y(2).gt.zmx ) exit
        enddo
        close(iunit)
        write(*,*) npoint, 'points within the domain'
    enddo

end subroutine make_poincare_plot

subroutine fieldline_derivative(phi,y,dy_dphi)
    real(dp), intent(in) :: phi
    real(dp), intent(in) :: y(2)
    real(dp), intent(out) :: dy_dphi(2)

    real(dp) :: rrr,ppp,zzz
    real(dp) :: Br,Bp, Bz

    rrr = y(1)
    zzz = y(2)
    ppp = phi
    ppp = set_to_period(ppp)

    
    call dummy_field(rrr, ppp, zzz, Br, Bp, Bz)
    dy_dphi(1) = Br*rrr/Bp
    dy_dphi(2) = Bz*rrr/Bp

    return
end subroutine fieldline_derivative

subroutine dummy_field(rrr, ppp, zzz, Br, Bp, Bz)
    real(dp), intent(in) :: rrr, ppp, zzz
    real(dp), intent(out) :: Br, Bp, Bz

    Br = 0.0_dp
    Bp = 0.0_dp
    Bz = 0.0_dp

    return
end subroutine dummy_field

function set_to_period(x)
    real(dp), intent(in) :: x
    real(dp) :: set_to_period

    set_to_period = x
    
    return
end function set_to_period

end module neo_poincare_plot
