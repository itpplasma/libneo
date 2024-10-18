module neo_poincare
use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_field_base, only: field_t
implicit none

type :: poincare_config_t
    integer :: n_fieldlines
    real(dp) :: fieldline_start_Rmin, fieldline_start_Rmax
    real(dp) :: fieldline_start_phi, fieldline_start_Z
    integer :: n_periods
    real(dp) :: period_length
    real(dp) :: integrate_err
    real(dp) :: plot_Rmin, plot_Rmax, plot_Zmin, plot_Zmax
end type poincare_config_t

real(dp), parameter :: pi = 3.14159265358979d0
integer, parameter :: poincare_dim = 2

contains


subroutine make_poincare(field, config, output_filename)

    class(field_t), intent(in) :: field
    type(poincare_config_t), intent(in) :: config
    character(len=*), intent(in), optional :: output_filename

    integer :: fieldline
    real(dp), dimension(:,:), allocatable :: R, Z
    real(dp), dimension(:), allocatable :: start_R, start_Z

    allocate(R(config%n_fieldlines, config%n_periods))
    allocate(Z(config%n_fieldlines, config%n_periods))
    allocate(start_R(config%n_fieldlines))
    allocate(start_Z(config%n_fieldlines))
    call get_fieldlines_startpoints(config, start_R, start_Z)
    do fieldline = 1, config%n_fieldlines
        R(fieldline, 1) = start_R(fieldline)
        Z(fieldline, 1) = start_Z(fieldline)
        write(*,*) 'fieldline #', fieldline
        call get_poincare_RZ(field, config, R(fieldline, :), Z(fieldline, :))
    enddo
    call write_poincare_RZ_to_file(R, Z, output_filename)
    deallocate(R, Z)
    deallocate(start_R, start_Z)
end subroutine make_poincare

subroutine get_fieldlines_startpoints(config, start_R, start_Z)
    type(poincare_config_t), intent(in) :: config
    real(dp), dimension(:), intent(out) :: start_R, start_Z

    real(dp) :: delta_R
    integer :: fieldline

    if (config%n_fieldlines .gt. 1) then
        delta_R = (config%fieldline_start_Rmax - config%fieldline_start_Rmin) &
                  / real(config%n_fieldlines - 1, kind=dp)
            do fieldline = 1, config%n_fieldlines
                start_R(fieldline) = config%fieldline_start_Rmin + &
                            (fieldline - 1) * delta_R
                start_Z(fieldline) = config%fieldline_start_Z
            end do
    elseif (config%n_fieldlines .eq. 1) then
        start_R(1) = config%fieldline_start_Rmin
        start_Z(1) = config%fieldline_start_Z
    else
        write(*,*) 'Error: n_fieldlines must be greater than 0'
        stop
    end if

end subroutine get_fieldlines_startpoints

subroutine get_poincare_RZ(field, config, R, Z)

    class(field_t), intent(in) :: field
    type(poincare_config_t), intent(in) :: config
    real(dp), dimension(:), intent(inout) :: R, Z

    real(dp) :: RZ(poincare_dim), phi, phi_end
    integer :: period

    RZ(1) = R(1)
    RZ(2) = Z(1)
    phi = config%fieldline_start_phi
    do period = 2, config%n_periods
        phi_end = phi + config%period_length
        call integrate_RZ_along_fieldline(field, RZ, phi, phi_end, config%integrate_err) 
        phi = phi_end
        if(is_in_plot_region(RZ, config)) then
            R(period) = RZ(1)
            Z(period) = RZ(2)
        else
            R(period:) = R(period-1)
            Z(period:) = Z(period-1)
            exit
        end if
        write(*,*) period, '/', config%n_periods, 'periods finished'
    enddo
end subroutine get_poincare_RZ

subroutine integrate_RZ_along_fieldline(field, RZ, phi_start, phi_end, relerr)
    use ode_integration, only: odeint_allroutines

    class(field_t), intent(in) :: field
    real(dp), intent(inout) :: RZ(poincare_dim)
    real(dp), intent(in) :: phi_start, phi_end
    real(dp), intent(in) :: relerr

    call odeint_allroutines(RZ, poincare_dim, phi_start, phi_end, &
                            relerr, fieldline_derivative, initial_stepsize=0.1_dp)

    contains 

        subroutine fieldline_derivative(phi, RZ, dRZ_dphi, ierr)
            real(dp), intent(in) :: phi
            real(dp), intent(in), dimension(:) :: RZ
            real(dp), intent(out), dimension(:) :: dRZ_dphi
            integer, intent(out) :: ierr
        
            real(dp), parameter :: tol = 1.0e-10_dp
            real(dp) :: RphiZ(3)
            real(dp) :: B(3)
        
            RphiZ(1) = RZ(1)
            RphiZ(2) = phi
            RphiZ(3) = RZ(2)
            call field%compute_bfield(RphiZ, B)
            if (B(2) .lt. tol) then
                print *, 'Error: B(2) vanishes'
                ierr = 1
                return
            endif
            dRZ_dphi(1) = B(1) * RphiZ(1) / B(2)
            dRZ_dphi(2) = B(3) * RphiZ(1) / B(2)
            ierr = 0
        end subroutine fieldline_derivative
end subroutine integrate_RZ_along_fieldline

function is_in_plot_region(RZ, config)
    real(dp), intent(in) :: RZ(poincare_dim)
    type(poincare_config_t), intent(in) :: config
    logical :: is_in_plot_region

    is_in_plot_region = .false.
    if (RZ(1) .ge. config%plot_Rmin .and. &
        RZ(1) .le. config%plot_Rmax .and. &
        RZ(2) .ge. config%plot_Zmin .and. &
        RZ(2) .le. config%plot_Zmax) is_in_plot_region = .true.
end function is_in_plot_region

subroutine write_poincare_RZ_to_file(R, Z, filename)
    real(dp), dimension(:,:), intent(in) :: R, Z
    character(*), intent(in), optional :: filename

    integer :: file_id
    integer :: n_fieldlines, fieldline, period, n_periods
    character(len=256) :: output_file


    if (present(filename)) then
        output_file = filename
    else
        output_file = 'poincare.dat'
    endif
    n_fieldlines = size(R, 1)
    n_periods = size(R, 2)
    open(newunit=file_id, file=output_file, status='replace')
    do fieldline = 1, n_fieldlines
        do period = 1, n_periods
            write(file_id,*) R(fieldline, period), Z(fieldline, period)
        end do
        write(file_id,*) ' '
    end do
    close(file_id)
end subroutine write_poincare_RZ_to_file


subroutine read_config_file(config, config_file)
    type(poincare_config_t), intent(out) :: config
    character(len=*), intent(in) :: config_file

    integer :: file_id
    integer :: n_periods, n_fieldlines
    real(dp) :: integrate_err, period_length
    real(dp) :: fieldline_start_Rmin, fieldline_start_Rmax
    real(dp) :: fieldline_start_phi, fieldline_start_Z
    real(dp) :: plot_Rmin, plot_Rmax, plot_Zmin, plot_Zmax

    namelist /poincare/ &
                n_fieldlines, &
                fieldline_start_Rmin, &
                fieldline_start_Rmax, &
                fieldline_start_phi, &
                fieldline_start_Z, &
                n_periods, &
                period_length, &
                plot_Rmin, &
                plot_Rmax, &
                plot_Zmin, &
                plot_Zmax, &
                integrate_err
    open(newunit=file_id, file=config_file, status='old')
    read(file_id, nml=poincare)
    close(file_id)
    config%n_fieldlines = n_fieldlines
    config%fieldline_start_Rmin = fieldline_start_Rmin
    config%fieldline_start_Rmax = fieldline_start_Rmax
    config%fieldline_start_phi = fieldline_start_phi
    config%fieldline_start_Z = fieldline_start_Z
    config%n_periods = n_periods
    config%period_length = period_length
    config%integrate_err = integrate_err
    config%plot_Rmin = plot_Rmin
    config%plot_Rmax = plot_Rmax
    config%plot_Zmin = plot_Zmin
    config%plot_Zmax = plot_Zmax
end subroutine read_config_file

end module neo_poincare
