module neo_poincare_plot
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none
    

contains


subroutine make_poincare_plot(field, Rmin, Rmax, Zmin, Zmax, idir, &
                              R_start, phi_start, Z_start, n_periods)
    use ode_integration, only: odeint_allroutines
    use neo_field_base, only: field_t

    class(field_t), intent(in) :: field
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax
    real(dp), intent(in) :: R_start, phi_start, Z_start
    integer, intent(in) :: idir, n_periods

    real(dp), parameter :: pi = 3.14159265358979d0
    real(dp), parameter :: per_phi = 2*pi
    integer, parameter :: neqn = 2
   
    real(dp) :: relerr
    integer :: file_id
    real(dp) :: phi, phiout
    integer :: period
    real(dp) :: y(2)
    real(dp) :: R, Z
    real(dp) :: Br, Bp, Bz

    relerr = 1.e-8

    open(newunit=file_id, file='poincare_plot.dat', status='replace')
    y(1) = R_start
    phi = phi_start
    y(2) = Z_start
    do period = 1, n_periods
        phiout = phi + per_phi*idir
        call odeint_allroutines(y, neqn, phi , phiout , relerr, fieldline_derivative) 
        phi = phiout
        R = y(1)
        phi = set_to_period(phi)
        Z = y(2)
        write(file_id,*) R, phi, Z
        if( R.lt.Rmin .or. R.gt.Rmax .or. Z.lt.Zmin .or. Z.gt.Zmax ) exit
    enddo
    close(file_id)
    write(*,*) period, 'finished periods of' , n_periods
end subroutine make_poincare_plot

subroutine fieldline_derivative(phi,y,dy_dphi,ierr)
    real(dp), intent(in) :: phi
    real(dp), intent(in), dimension(:) :: y
    real(dp), intent(out), dimension(:) :: dy_dphi
    integer, intent(out) :: ierr

    real(dp) :: R, Z
    real(dp) :: Br,Bp, Bz

    R = y(1)
    Z = y(2)
    call dummy_field(R, set_to_period(phi), Z, Br, Bp, Bz)
    dy_dphi(1) = Br*R/Bp
    dy_dphi(2) = Bz*R/Bp
    ierr = 0
end subroutine fieldline_derivative

subroutine dummy_field(R, phi, Z, Br, Bp, Bz)
    real(dp), intent(in) :: R, phi, Z
    real(dp), intent(out) :: Br, Bp, Bz

    Br = 0.0_dp
    Bp = 0.0_dp
    Bz = 1.0_dp
end subroutine dummy_field

function set_to_period(x)
    real(dp), intent(in) :: x
    real(dp) :: set_to_period

    set_to_period = x
end function set_to_period

end module neo_poincare_plot
