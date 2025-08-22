program benchmark_poincare
    use neo_poincare
    use neo_field_base, only: field_t
    use neo_example_field, only: example_field_t
    use libneo_kinds, only: dp
    implicit none

    type(example_field_t) :: field
    type(poincare_config_t) :: config
    integer :: i, j, n_runs
    real(dp) :: start_time, end_time, total_time
    real(dp), dimension(:,:), allocatable :: R, Z
    real(dp), dimension(:), allocatable :: start_R, start_Z

    ! Initialize field
    call field%example_field_init()

    ! Configure test parameters
    config%n_fieldlines = 5
    config%fieldline_start_Rmin = 1.0_dp
    config%fieldline_start_Rmax = 1.5_dp
    config%fieldline_start_phi = 0.0_dp
    config%fieldline_start_Z = 0.0_dp
    config%n_periods = 10
    config%period_length = 0.1_dp
    config%integrate_err = 1.0e-8_dp
    config%plot_Rmin = 0.5_dp
    config%plot_Rmax = 2.0_dp
    config%plot_Zmin = -1.0_dp
    config%plot_Zmax = 1.0_dp

    ! Allocate arrays
    allocate(R(config%n_fieldlines, config%n_periods))
    allocate(Z(config%n_fieldlines, config%n_periods))
    allocate(start_R(config%n_fieldlines))
    allocate(start_Z(config%n_fieldlines))

    ! Number of benchmark runs
    n_runs = 100

    write(*,*) 'Running Poincare benchmark...'
    write(*,*) 'Number of runs:', n_runs
    write(*,*) 'Fieldlines:', config%n_fieldlines
    write(*,*) 'Periods:', config%n_periods

    ! Warm-up run
    call get_fieldlines_startpoints(config, start_R, start_Z)
    do j = 1, config%n_fieldlines
        R(j, 1) = start_R(j)
        Z(j, 1) = start_Z(j)
        call get_poincare_RZ(field, config, R(j, :), Z(j, :))
    end do

    ! Timed benchmark
    call cpu_time(start_time)
    do i = 1, n_runs
        call get_fieldlines_startpoints(config, start_R, start_Z)
        do j = 1, config%n_fieldlines
            R(j, 1) = start_R(j)
            Z(j, 1) = start_Z(j)
            call get_poincare_RZ(field, config, R(j, :), Z(j, :))
        end do
    end do
    call cpu_time(end_time)

    total_time = end_time - start_time
    write(*,*) 'Total time (seconds):', total_time
    write(*,*) 'Average time per run (ms):', (total_time / n_runs) * 1000.0_dp
    write(*,*) 'Throughput (runs/second):', n_runs / total_time

    deallocate(R, Z, start_R, start_Z)

end program benchmark_poincare