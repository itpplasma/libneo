program test_poincare
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

real(dp), parameter :: pi = 3.14159265358979_dp
character(len=*), parameter :: config_file = 'poincare.inp'

!call write_poincare_config_jorek_field
!call test_make_poincare
!call test_make_poincare_flux_pumping
!call remove_poincare_input
call write_poincare_config_circular_tokamak_field
call test_make_poincare_circular_tokamak_field
call remove_poincare_input

contains


subroutine write_poincare_config_jorek_field
    use neo_poincare, only: poincare_config_t

    type(poincare_config_t) :: config

    call print_test("write_poincare_config_jorek_field")
    config%n_fieldlines = 10
    config%fieldline_start_Rmin = 1.75_dp
    config%fieldline_start_Rmax = 2.0_dp
    config%fieldline_start_phi = 0.0_dp
    config%fieldline_start_Z = -0.2_dp
    config%n_periods = 10
    config%period_length = 2.0_dp * pi
    config%integrate_err = 1.0e-8_dp
    config%plot_Rmin = 1.0_dp
    config%plot_Rmax = 2.25_dp
    config%plot_Zmin = -1.0_dp
    config%plot_Zmax = 1.0_dp
    call write_poincare_config (config)
    call print_ok
end subroutine write_poincare_config_jorek_field

subroutine write_poincare_config_circular_tokamak_field
    use neo_poincare, only: poincare_config_t

    type(poincare_config_t) :: config

    call print_test("write_poincare_config_circular_tokamak_field")
    config%n_fieldlines = 10
    config%fieldline_start_Rmin = 1.1_dp
    config%fieldline_start_Rmax = 1.5_dp
    config%fieldline_start_phi = 0.0_dp
    config%fieldline_start_Z = 0.0_dp
    config%n_periods = 50
    config%period_length = 2.0_dp * pi
    config%integrate_err = 1.0e-8_dp
    config%plot_Rmin = 0.2_dp
    config%plot_Rmax = 2.25_dp
    config%plot_Zmin = -0.8_dp
    config%plot_Zmax = 0.8_dp

    call write_poincare_config(config)
    call print_ok
end subroutine write_poincare_config_circular_tokamak_field

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


subroutine test_make_poincare
    use neo_poincare, only: make_poincare, read_config_file, poincare_config_t
    use neo_jorek_field, only: jorek_field_t
    use util_for_test_jorek_field, only: get_filename, filename_len

    character(len=filename_len) :: jorek_file
    type(jorek_field_t) :: field
    type(poincare_config_t) :: config

    call print_test("test_make_poincare")

    call get_filename(jorek_file)
    call field%jorek_field_init(jorek_file)
    call read_config_file(config, config_file)
    call make_poincare(field, config)
    call print_ok
end subroutine test_make_poincare


subroutine test_make_poincare_flux_pumping
    use neo_poincare, only: make_poincare, read_config_file, poincare_config_t
    use neo_jorek_field, only: jorek_field_t

    character(len=512) :: jorek_file
    type(jorek_field_t) :: field
    type(poincare_config_t) :: config

    call print_test("test_make_poincare_flux_pumping")

    jorek_file ="/proj/plasma/DATA/AUG/JOREK/2024-05_test_haowei_flux_pumping/" // &
    "exprs_Rmin1.140_Rmax2.130_Zmin-0.921_Zmax0.778_phimin0.000_phimax6.283_s40000.h5"

    call field%jorek_field_init(jorek_file)
    call read_config_file(config, config_file)
    call make_poincare(field, config)
    call print_ok
end subroutine test_make_poincare_flux_pumping


subroutine test_make_poincare_circular_tokamak_field()
    use neo_poincare, only: make_poincare, read_config_file, poincare_config_t
    use neo_circular_tokamak_field, only: circular_tokamak_field_t

    type(circular_tokamak_field_t) :: field
    type(poincare_config_t) :: config

    call print_test("test_make_poincare_circular_tokamak_field")
    call field%circular_tokamak_field_init()
    call read_config_file(config, config_file)
    call make_poincare(field, config)

    call print_ok
end subroutine test_make_poincare_circular_tokamak_field


subroutine remove_poincare_input
    integer :: stat, file_id
    logical :: exists

    call print_test("remove_poincare_input")
    
    open(newunit=file_id, iostat=stat, file=config_file, status='old')
    if (stat == 0) close(file_id, status='delete')
    inquire(file=trim(config_file), exist=exists)
    if (exists) then
        call print_fail
        error stop
    else
        call print_ok
    end if
end subroutine remove_poincare_input

end program test_poincare