program export_boozer_signgs_probe
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates
    use boozer_chartmap, only: export_boozer_chartmap

    implicit none

    character(len=1024) :: wout_file, chartmap_file

    if (command_argument_count() /= 2) then
        print *, 'Usage: export_boozer_signgs_probe.x <wout.nc> <chartmap.nc>'
        error stop
    end if

    call get_command_argument(1, wout_file)
    call get_command_argument(2, chartmap_file)

    use_B_r = .false.

    call get_boozer_coordinates(trim(wout_file), radial_spline_order=5, &
        angular_spline_order=5, grid_refinment=3)
    call export_boozer_chartmap(trim(chartmap_file))
end program export_boozer_signgs_probe
