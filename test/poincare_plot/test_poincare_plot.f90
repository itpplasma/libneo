program test_poincare_plot
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
implicit none

real(dp), parameter :: pi = 3.14159265358979_dp
character(len=*), parameter :: config_file = 'poincare_plot.inp'

call make_poincare_plot_input
call test_poincare_plot_call
call remove_poincare_plot_input

contains

subroutine make_poincare_plot_input
    integer :: n_fieldlines = 10
    real(dp) :: fieldline_start_Rmin = 1.75_dp, fieldline_start_Rmax = 2.0_dp
    real(dp) :: fieldline_start_phi = 0.0_dp, fieldline_start_Z = 0.0_dp
    integer :: n_periods = 100
    real(dp) :: period_length = 2.0_dp * pi
    real(dp) :: integrate_err = 1.0e-6_dp
    real(dp) :: plot_Rmin = 1.0_dp, plot_Rmax = 2.25_dp
    real(dp) :: plot_Zmin = -1.0_dp, plot_Zmax = 1.0_dp

    integer :: file_id

    namelist /poincare_plot/ &
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

    call print_test("make_poincare_plot_input")
    open(newunit=file_id, file=config_file, status='unknown')
    write(file_id, nml=poincare_plot)
    close(file_id)
    call print_ok
end subroutine make_poincare_plot_input

subroutine test_poincare_plot_call
    use neo_poincare_plot, only: make_poincare_plot

    logical :: exists

    call print_test("test_poincare_plot_call")

    inquire(file=trim(config_file), exist=exists)
    if (.not. exists) then
        call print_fail
        error stop
    else
        call print_ok
    endif
end subroutine test_poincare_plot_call

subroutine remove_poincare_plot_input
    integer :: stat, file_id
    logical :: exists

    call print_test("remove_poincare_plot_input")
    
    open(newunit=file_id, iostat=stat, file=config_file, status='old')
    if (stat == 0) close(file_id, status='delete')
    inquire(file=trim(config_file), exist=exists)
    if (exists) then
        call print_fail
        error stop
    else
        call print_ok
    end if
end subroutine remove_poincare_plot_input

end program test_poincare_plot