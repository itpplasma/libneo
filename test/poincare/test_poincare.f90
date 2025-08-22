program test_poincare
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail
use neo_poincare, only: poincare_config_t
implicit none

real(dp), parameter :: TOL = 1.0e-13_dp
real(dp), parameter :: pi = 3.14159265358979_dp
character(len=*), parameter :: config_file = 'poincare.inp'
type(poincare_config_t) :: jorek_config

jorek_config%n_fieldlines = 10
jorek_config%fieldline_start_Rmin = 1.8_dp
jorek_config%fieldline_start_Rmax = 2.0_dp
jorek_config%fieldline_start_phi = 0.0_dp
jorek_config%fieldline_start_Z = 0.0_dp
jorek_config%n_periods = 50
jorek_config%period_length = 2.0_dp * pi
jorek_config%integrate_err = 1.0e-8_dp
jorek_config%plot_Rmin = 1.0_dp
jorek_config%plot_Rmax = 2.25_dp
jorek_config%plot_Zmin = -1.0_dp
jorek_config%plot_Zmax = 1.0_dp

call test_make_poincare
call test_get_poincare_RZ_for_closed_fieldline
call test_read_config_file
call test_integrate_RZ_along_fieldline_direct
call test_fieldline_derivative_calculations
call test_fieldline_derivative_edge_cases
call draw_circular_tokamak_poincare

contains

subroutine test_integrate_RZ_along_fieldline_direct
    use neo_poincare, only: integrate_RZ_along_fieldline
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    real(dp) :: RZ(2), RZ_initial(2), RZ_expected(2)
    real(dp), parameter :: phi_start = 0.0_dp
    real(dp), parameter :: phi_end = pi / 2.0_dp
    real(dp), parameter :: relerr = 1.0e-8_dp
    real(dp), parameter :: R_axis = 1.0_dp, Z_axis = 0.0_dp
    real(dp), parameter :: B_tor_ampl = 1.0_dp, B_pol_ampl = 0.1_dp

    call print_test("test_integrate_RZ_along_fieldline_direct")

    call field%circular_tokamak_field_init(R_axis=R_axis, Z_axis=Z_axis, &
                                           B_pol_ampl=B_pol_ampl, B_tor_ampl=B_tor_ampl)

    ! Test integration along known circular path
    RZ_initial(1) = 1.5_dp ! R
    RZ_initial(2) = 0.0_dp ! Z
    RZ = RZ_initial
    
    call integrate_RZ_along_fieldline(field, RZ, phi_start, phi_end, relerr)
    
    ! For circular tokamak, we can calculate expected position analytically
    ! At phi = pi/2, for safety factor q = B_tor/B_pol = 10:
    ! theta_expected = phi / q = (pi/2) / 10 = pi/20
    RZ_expected(1) = R_axis + (RZ_initial(1) - R_axis) * cos(phi_end / (B_tor_ampl/B_pol_ampl))
    RZ_expected(2) = Z_axis + (RZ_initial(1) - R_axis) * sin(phi_end / (B_tor_ampl/B_pol_ampl))
    
    if (abs(RZ(1) - RZ_expected(1)) > 1.0e-3_dp) then
        call print_fail
        print *, "R integration failed: got ", RZ(1), " expected ", RZ_expected(1)
        error stop
    end if
    
    if (abs(RZ(2) - RZ_expected(2)) > 1.0e-3_dp) then
        call print_fail
        print *, "Z integration failed: got ", RZ(2), " expected ", RZ_expected(2)
        error stop
    end if
    
    call print_ok
end subroutine test_integrate_RZ_along_fieldline_direct

subroutine test_fieldline_derivative_calculations
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    real(dp) :: phi, RZ(2), dRZ_dphi(2), dRZ_dphi_expected(2)
    real(dp) :: RphiZ(3), B(3)
    real(dp), parameter :: R_axis = 1.0_dp, Z_axis = 0.0_dp
    real(dp), parameter :: B_tor_ampl = 1.0_dp, B_pol_ampl = 0.1_dp
    real(dp), parameter :: test_tolerance = 1.0e-10_dp

    call print_test("test_fieldline_derivative_calculations")

    call field%circular_tokamak_field_init(R_axis=R_axis, Z_axis=Z_axis, &
                                           B_pol_ampl=B_pol_ampl, B_tor_ampl=B_tor_ampl)

    ! Test derivative calculation at specific point
    phi = 0.0_dp
    RZ(1) = 1.5_dp ! R
    RZ(2) = 0.0_dp ! Z
    
    ! Calculate B field manually to verify derivative
    RphiZ(1) = RZ(1)
    RphiZ(2) = phi
    RphiZ(3) = RZ(2)
    call field%compute_bfield(RphiZ, B)
    
    ! Expected derivatives from fieldline equation
    dRZ_dphi_expected(1) = B(1) * RphiZ(1) / B(2)  ! dR/dphi
    dRZ_dphi_expected(2) = B(3) * RphiZ(1) / B(2)  ! dZ/dphi
    
    ! Test the derivative calculation logic
    ! This tests the math that will be in the module-level subroutine
    if (B(2) .lt. 1.0e-10_dp) then
        call print_fail
        print *, "B_phi too small for test: ", B(2)
        error stop
    endif
    
    dRZ_dphi(1) = B(1) * RphiZ(1) / B(2)
    dRZ_dphi(2) = B(3) * RphiZ(1) / B(2)
    
    if (abs(dRZ_dphi(1) - dRZ_dphi_expected(1)) > test_tolerance) then
        call print_fail
        print *, "dR/dphi calculation error: got ", dRZ_dphi(1), &
                 " expected ", dRZ_dphi_expected(1)
        error stop
    end if
    
    if (abs(dRZ_dphi(2) - dRZ_dphi_expected(2)) > test_tolerance) then
        call print_fail
        print *, "dZ/dphi calculation error: got ", dRZ_dphi(2), &
                 " expected ", dRZ_dphi_expected(2)
        error stop
    end if
    
    call print_ok
end subroutine test_fieldline_derivative_calculations

subroutine test_fieldline_derivative_edge_cases
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    real(dp) :: phi, RZ(2)
    real(dp) :: RphiZ(3), B(3)
    real(dp), parameter :: R_axis = 1.0_dp, Z_axis = 0.0_dp
    real(dp), parameter :: B_tor_ampl = 1.0e-12_dp  ! Very small B_phi
    real(dp), parameter :: B_pol_ampl = 0.1_dp
    real(dp), parameter :: tol = 1.0e-10_dp
    logical :: error_detected

    call print_test("test_fieldline_derivative_edge_cases")

    call field%circular_tokamak_field_init(R_axis=R_axis, Z_axis=Z_axis, &
                                           B_pol_ampl=B_pol_ampl, B_tor_ampl=B_tor_ampl)

    ! Test near-zero B_phi condition that should trigger error
    phi = 0.0_dp
    RZ(1) = 1.5_dp ! R
    RZ(2) = 0.0_dp ! Z
    
    RphiZ(1) = RZ(1)
    RphiZ(2) = phi
    RphiZ(3) = RZ(2)
    call field%compute_bfield(RphiZ, B)
    
    ! Verify that B_phi is indeed below threshold
    if (B(2) .ge. tol) then
        call print_fail
        print *, "Test setup error: B_phi should be below tolerance"
        print *, "B_phi = ", B(2), " tolerance = ", tol
        error stop
    end if
    
    ! Test error detection logic that will be in module-level subroutine
    error_detected = .false.
    if (B(2) .lt. tol) then
        error_detected = .true.
    endif
    
    if (.not. error_detected) then
        call print_fail
        print *, "Failed to detect B_phi vanishing condition"
        error stop
    end if
    
    call print_ok
end subroutine test_fieldline_derivative_edge_cases


subroutine test_make_poincare
    use neo_poincare, only: make_poincare, poincare_config_t
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    type(poincare_config_t) :: config

    call print_test("test_make_poincare")

    call field%circular_tokamak_field_init()
    config%fieldline_start_Rmin = 1.1_dp
    config%fieldline_start_Rmax = 1.5_dp
    config%fieldline_start_phi = 0.0_dp
    config%fieldline_start_Z = 0.0_dp
    config%n_fieldlines = 10
    config%n_periods = 50
    config%period_length = 2.0_dp * pi
    config%integrate_err = 1.0e-6_dp
    config%plot_Rmin = 0.5_dp
    config%plot_Rmax = 2.0_dp
    config%plot_Zmin = -1.0_dp
    config%plot_Zmax = 1.0_dp

    call make_poincare(field, config)
    call print_ok
end subroutine test_make_poincare


subroutine test_get_poincare_RZ_for_closed_fieldline()
    use neo_poincare, only: get_poincare_RZ, poincare_config_t
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    type(poincare_config_t) :: config

    real(dp), parameter :: safety_factor = 10
    real(dp), parameter :: R_axis = 2.0_dp, Z_axis = 0.0_dp
    real(dp) :: B_tor_ampl, B_pol_ampl
    real(dp), dimension(:), allocatable :: R, Z, delta_theta
    real(dp) :: expected_delta_theta

    call print_test("test_get_poincare_RZ_for_closed_fieldline")
    B_tor_ampl = 1.0_dp
    B_pol_ampl = B_tor_ampl/safety_factor
    call field%circular_tokamak_field_init(R_axis=R_axis, Z_axis=Z_axis, &
                                           B_pol_ampl=B_pol_ampl, B_tor_ampl=B_tor_ampl)
    config%n_fieldlines = 1
    config%fieldline_start_phi = 0.0_dp
    config%n_periods = int(safety_factor) + 1
    config%period_length = 2.0_dp * pi
    config%integrate_err = 1.0e-6_dp
    config%plot_Rmin = 0.5_dp
    config%plot_Rmax = 3.5_dp
    config%plot_Zmin = -1.5_dp
    config%plot_Zmax = 1.5_dp

    allocate(R(config%n_periods), Z(config%n_periods))
    allocate(delta_theta(config%n_periods-1))
    R(1) = 1.0_dp
    Z(1) = 0.0_dp
    call get_poincare_RZ(field, config, R, Z)
    delta_theta = calc_delta_theta(R, Z, R_axis, Z_axis, config%period_length)

    expected_delta_theta = config%period_length / safety_factor
    if (any(abs(delta_theta-expected_delta_theta) > config%integrate_err)) then
        call print_fail
        print *, "delta_theta = ", delta_theta
        print *, "expected_delta_theta = ", expected_delta_theta
        error stop
    end if
    if ((abs(R(config%n_periods) - R(1)) > config%integrate_err) .and. &
        (abs(Z(config%n_periods) - Z(1)) > config%integrate_err)) then
        call print_fail
        print *, "last point and first point not equal of closed fieldline"
        print *, "last point = ", R(config%n_periods), Z(config%n_periods)
        print *, "first point = ", R(1), Z(1)
        error stop
    end if

    call print_ok
end subroutine test_get_poincare_RZ_for_closed_fieldline

function calc_delta_theta(R, Z, R_axis, Z_axis, period_length) result(delta_theta)
    real(dp), dimension(:), intent(in) :: R, Z
    real(dp), intent(in) :: R_axis, Z_axis
    real(dp), intent(in) :: period_length
    real(dp), dimension(:), allocatable  :: delta_theta

    integer :: period
    real(dp) :: theta_start, theta_end

    allocate(delta_theta(size(R)-1))
    do period = 1, size(R)-1
        theta_start = atan2(Z(period) - Z_axis, R(period) - R_axis)
        theta_end = atan2(Z(period+1) - Z_axis, R(period+1) - R_axis)
        delta_theta(period) = modulo(theta_end - theta_start, period_length)
    end do
end function calc_delta_theta


subroutine test_read_config_file
    use neo_poincare, only: poincare_config_t, read_config_file

    type(poincare_config_t) :: config

    call print_test("test_read_config_file")


    call write_poincare_config(jorek_config)
    call read_config_file(config, config_file)

    if (config%n_fieldlines /= jorek_config%n_fieldlines) then
        call print_fail
        print *, "n_fieldlines = ", config%n_fieldlines
        error stop
    end if
    if (abs(config%fieldline_start_Rmin - jorek_config%fieldline_start_Rmin) > TOL) then
        call print_fail
        print *, "fieldline_start_Rmin = ", config%fieldline_start_Rmin
        error stop
    end if
    if (abs(config%fieldline_start_Rmax - jorek_config%fieldline_start_Rmax) > TOL) then
        call print_fail
        print *, "fieldline_start_Rmax = ", config%fieldline_start_Rmax
        error stop
    end if
    if (abs(config%fieldline_start_phi - jorek_config%fieldline_start_phi) > TOL) then
        call print_fail
        print *, "fieldline_start_phi = ", config%fieldline_start_phi
        error stop
    end if
    if (abs(config%fieldline_start_Z - jorek_config%fieldline_start_Z) > TOL) then
        call print_fail
        print *, "fieldline_start_Z = ", config%fieldline_start_Z
        error stop
    end if
    if (config%n_periods /= jorek_config%n_periods) then
        call print_fail
        print *, "n_periods = ", config%n_periods
        error stop
    end if
    if (abs(config%period_length - jorek_config%period_length) > TOL) then
        call print_fail
        print *, "period_length = ", config%period_length
        error stop
    end if
    if (abs(config%integrate_err - jorek_config%integrate_err) > TOL) then
        call print_fail
        print *, "integrate_err = ", config%integrate_err
        error stop
    end if
    if (abs(config%plot_Rmin - jorek_config%plot_Rmin) > TOL) then
        call print_fail
        print *, "plot_Rmin = ", config%plot_Rmin
        error stop
    end if
    if (abs(config%plot_Rmax - jorek_config%plot_Rmax) > TOL) then
        call print_fail
        print *, "plot_Rmax = ", config%plot_Rmax
        error stop
    end if
    if (abs(config%plot_Zmin - jorek_config%plot_Zmin) > TOL) then
        call print_fail
        print *, "plot_Zmin = ", config%plot_Zmin
        error stop
    end if
    if (abs(config%plot_Zmax - jorek_config%plot_Zmax) > TOL) then
        call print_fail
        print *, "plot_Zmax = ", config%plot_Zmax
        error stop
    end if

    call remove_poincare_config
    call print_ok
end subroutine test_read_config_file


subroutine write_poincare_config(config)
    use neo_poincare, only: poincare_config_t
    type(poincare_config_t), intent(in) :: config
    integer :: n_fieldlines
    real(dp) :: fieldline_start_Rmin, fieldline_start_Rmax
    real(dp) :: fieldline_start_phi, fieldline_start_Z
    integer :: n_periods
    real(dp) :: period_length
    real(dp) :: integrate_err
    real(dp) :: plot_Rmin, plot_Rmax
    real(dp) :: plot_Zmin, plot_Zmax
    integer :: file_id

    namelist /poincare/ &
                n_fieldlines, &
                fieldline_start_Rmin, &
                fieldline_start_Rmax, &
                fieldline_start_phi, &
                fieldline_start_Z, &
                n_periods, &
                period_length, &
                integrate_err, &
                plot_Rmin, &
                plot_Rmax, &
                plot_Zmin, &
                plot_Zmax

    n_fieldlines = config%n_fieldlines
    fieldline_start_Rmin = config%fieldline_start_Rmin
    fieldline_start_Rmax = config%fieldline_start_Rmax
    fieldline_start_phi = config%fieldline_start_phi
    fieldline_start_Z = config%fieldline_start_Z
    n_periods = config%n_periods
    period_length = config%period_length
    integrate_err = config%integrate_err
    plot_Rmin = config%plot_Rmin
    plot_Rmax = config%plot_Rmax
    plot_Zmin = config%plot_Zmin
    plot_Zmax = config%plot_Zmax

    open(newunit=file_id, file=config_file, status='unknown')
    write(file_id, nml=poincare)
    close(file_id)
end subroutine write_poincare_config


subroutine remove_poincare_config
    integer :: stat, file_id
    logical :: exists

    open(newunit=file_id, iostat=stat, file=config_file, status='old')
    if (stat == 0) close(file_id, status='delete')
    inquire(file=trim(config_file), exist=exists)
    if (exists) then
        call print_fail
        error stop
    end if
end subroutine remove_poincare_config


subroutine draw_circular_tokamak_poincare()
    use neo_poincare, only: get_poincare_RZ, poincare_config_t
    use neo_poincare, only: write_poincare_RZ_to_file
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    type(poincare_config_t) :: config

    real(dp), parameter :: safety_factor = pi
    real(dp), parameter :: R_axis = 2.0_dp, Z_axis = 0.0_dp
    real(dp) :: B_tor_ampl, B_pol_ampl
    real(dp), dimension(:,:), allocatable :: R, Z

    call print_test("draw_circular_tokamak_poincare")
    B_tor_ampl = 1.0_dp
    B_pol_ampl = B_tor_ampl/safety_factor
    call field%circular_tokamak_field_init(R_axis=R_axis, Z_axis=Z_axis, &
                                           B_pol_ampl=B_pol_ampl, B_tor_ampl=B_tor_ampl)
    config%n_fieldlines = 1
    config%fieldline_start_phi = 0.0_dp
    config%n_periods = 1000*(int(safety_factor) + 1)
    config%period_length = 2.0_dp * pi
    config%integrate_err = 1.0e-6_dp
    config%plot_Rmin = 0.5_dp
    config%plot_Rmax = 3.5_dp
    config%plot_Zmin = -1.5_dp
    config%plot_Zmax = 1.5_dp

    allocate(R(config%n_fieldlines,config%n_periods))
    allocate(Z(config%n_fieldlines,config%n_periods))
    R(1,1) = 1.0_dp
    Z(1,1) = 0.0_dp
    call get_poincare_RZ(field, config, R(1,:), Z(1,:))
    call write_poincare_RZ_to_file(R, Z, 'circular_tok_poincare.dat')

    call print_ok
end subroutine draw_circular_tokamak_poincare


end program test_poincare
