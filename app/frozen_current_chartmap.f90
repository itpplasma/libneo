!> Frozen-current field on the nodes of a chartmap.
!>
!> Usage: frozen_current_chartmap.x <field.h5> <chartmap.nc> <output.nc>
!>
!> Evaluates the field saved by frozen_current.x (field_file) at every node
!> of the chart and writes the chart with A and B (libneo_chartmap_field),
!> so orbit codes spline it directly instead of evaluating coils at start.
program frozen_current_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime
    use neo_frozen_current_field, only: frozen_current_field_t
    use neo_frozen_current_io, only: load_frozen_current_field
    use libneo_chartmap_field, only: chartmap_field_t, sample_chartmap_field, &
                                     write_chartmap_field
    implicit none

    character(len=1024) :: field_file, chartmap_file, output
    type(frozen_current_field_t) :: field
    type(chartmap_field_t) :: f
    integer :: status
    real(dp) :: t0

    if (command_argument_count() /= 3) then
        print *, 'usage: frozen_current_chartmap.x <field.h5> <chartmap.nc> <output.nc>'
        error stop
    end if
    call get_command_argument(1, field_file)
    call get_command_argument(2, chartmap_file)
    call get_command_argument(3, output)

    call load_frozen_current_field(trim(field_file), field, status)
    if (status /= 0) then
        print *, 'frozen_current_chartmap: cannot load ', trim(field_file)
        error stop
    end if
    t0 = omp_get_wtime()
    call sample_chartmap_field(field, trim(chartmap_file), trim(field_file), f)
    call write_chartmap_field(trim(chartmap_file), f, trim(output))
    print '(a,3i6,a,f9.2,a)', ' nodes (rho, theta, zeta):', shape(f%b(1, :, :, :)), &
        '   time', omp_get_wtime() - t0, ' s'
    print '(a,es10.3,a,es10.3)', ' spline error at cell centres: |dB|/|B| ', &
        f%spline_error_b, '   |dA|/max|A| ', f%spline_error_a
end program frozen_current_chartmap
