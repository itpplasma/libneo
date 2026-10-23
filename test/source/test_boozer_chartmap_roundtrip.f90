!> Pin field equivalence across the Boozer chartmap export/import roundtrip.
!>
!> Builds the Boozer field from the LandremanPaul2021 QA reference wout via
!> get_boozer_coordinates, records |B| at interior test points through
!> splint_boozer_coord, exports a chartmap NetCDF file, reimports it with
!> load_boozer_from_chartmap, and re-evaluates |B| at the same points. The
!> maximum relative |B| error over the interior points must stay below 8e-5.
program test_boozer_chartmap_roundtrip
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use new_vmec_stuff_mod, only: nper
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          sfl_quality_config_t, boozer_regular_radial_factor
    use netcdf, only: nf90_open, nf90_nowrite, nf90_global, nf90_get_att, &
                      nf90_close, nf90_noerr, nf90_enotatt, nf90_create, nf90_clobber, &
                      nf90_put_att, nf90_enddef, nf90_def_dim, nf90_def_var, &
                      nf90_double, nf90_int, nf90_put_var
    use boozer_chartmap_io, only: read_boozer_chartmap
    use boozer_chartmap_types, only: boozer_chartmap_data_t
    use boozer_chartmap, only: export_boozer_chartmap, load_boozer_from_chartmap, &
        set_boozer_chartmap_spline_orders
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  chartmap_from_cyl_ok
    use cylindrical_cartesian, only: cart_to_cyl
    implicit none

    real(dp), parameter :: PI = 3.14159265358979_dp
    real(dp), parameter :: TWOPI = 2.0_dp*PI
    real(dp), parameter :: field_tol = 8.0e-5_dp
    character(len=*), parameter :: wout_file = "wout.nc"
    character(len=*), parameter :: chartmap_file = "roundtrip_test.nc"

    integer, parameter :: n_test = 50
    real(dp) :: s_test(n_test), th_test(n_test), ph_test(n_test)
    real(dp) :: bmod_ref(n_test), bmod_new(n_test)

    real(dp) :: A_theta, A_phi_val, dA_theta_dr, dA_phi_dr
    real(dp) :: d2A_phi_dr2, d3A_phi_dr3
    real(dp) :: Bth, dBth, d2Bth, Bph, dBph, d2Bph
    real(dp) :: Bmod, dBmod(3), d2Bmod(6), Br, dBr(3), d2Br(6)
    real(dp) :: phi_period, rel_err, max_err_bmod
    integer :: i
    integer :: ncid, nc_status, overridden, radial_factor, exit_status
    character(len=256) :: override_reason
    character(len=8) :: construction
    character(len=1024) :: executable, argument
    type(sfl_quality_config_t) :: quality
    type(boozer_chartmap_data_t) :: imported

    integer, parameter :: n_coord_chk = 3
    class(coordinate_system_t), allocatable :: coord_cs
    real(dp) :: u_ref_chk(3, n_coord_chk), u_back_before(3, n_coord_chk)

    if (command_argument_count() /= 0) then
        call get_command_argument(1, argument)
        if (trim(argument) /= '--invalid-factor') error stop 'unexpected test argument'
        call read_boozer_chartmap("roundtrip_policy.nc", imported)
        stop
    end if
    use_B_r = .false.
    quality%tol(1) = 0.0_dp
    quality%allow_poor_equilibrium = .true.
    quality%override_reason = 'chartmap roundtrip tests recorded quality override'
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3, quality=quality)

    phi_period = TWOPI/real(nper, dp)
    call init_test_points(phi_period)

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
                                 A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, Bmod, dBmod, d2Bmod, &
                                 Br, dBr, d2Br)
        bmod_ref(i) = Bmod
    end do

    call export_boozer_chartmap(chartmap_file)
    nc_status = nf90_open(chartmap_file, nf90_nowrite, ncid)
    if (nc_status /= nf90_noerr) error stop 'cannot reopen chartmap'
    nc_status = nf90_get_att(ncid, nf90_global, &
                             'sfl_allow_poor_equilibrium', overridden)
    if (nc_status /= nf90_noerr) error stop 'missing chartmap quality override'
    if (overridden /= 1) error stop 'chartmap quality override was not recorded'
    nc_status = nf90_get_att(ncid, nf90_global, &
                             'sfl_override_reason', override_reason)
    if (nc_status /= nf90_noerr) error stop 'missing chartmap override reason'
    if (trim(override_reason) /= trim(quality%override_reason)) &
        error stop 'chartmap override reason changed'
    nc_status = nf90_get_att(ncid, nf90_global, 'regular_radial_factor', radial_factor)
    if (nc_status /= nf90_noerr .or. radial_factor /= 2) &
        error stop 'hybrid chartmap lost regular radial policy'
    nc_status = nf90_close(ncid)
    if (nc_status /= nf90_noerr) error stop 'cannot close chartmap'
    call check_coordinate_inverse_before_load(chartmap_file)
    call set_boozer_chartmap_spline_orders(3, 3)
    call load_boozer_from_chartmap(chartmap_file)
    if (boozer_regular_radial_factor /= 2) error stop 'chartmap load lost regular radial policy'

    do i = 1, n_test
        call splint_boozer_coord(s_test(i), th_test(i), ph_test(i), 0, &
                                 A_theta, A_phi_val, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, Bth, dBth, d2Bth, &
                                 Bph, dBph, d2Bph, Bmod, dBmod, d2Bmod, &
                                 Br, dBr, d2Br)
        bmod_new(i) = Bmod
    end do

    max_err_bmod = 0.0_dp
    do i = 1, n_test
        if (abs(bmod_ref(i)) > 0.0_dp) then
            rel_err = abs(bmod_new(i) - bmod_ref(i))/abs(bmod_ref(i))
        else
            rel_err = abs(bmod_new(i))
        end if
        max_err_bmod = max(max_err_bmod, rel_err)
    end do

    print *, "max relative |B| error over chartmap roundtrip:", max_err_bmod
    if (max_err_bmod > field_tol) then
        error stop "Bmod chartmap roundtrip error too large"
    end if
    print *, "Boozer chartmap export/import preserves |B| within tolerance"

    call check_coordinate_inverse_after_load()

    ! Importing a chartmap replaces the field without importing VMEC quality.
    ! Its export must not retain the prior field's explicit override.
    call export_boozer_chartmap('boozer_reexport.chartmap.nc')
    nc_status = nf90_open('boozer_reexport.chartmap.nc', nf90_nowrite, ncid)
    if (nc_status /= nf90_noerr) error stop 'cannot reopen imported chartmap export'
    nc_status = nf90_get_att(ncid, nf90_global, 'sfl_construction', construction)
    if (nc_status /= nf90_noerr .or. construction /= 'chartmap') &
        error stop 'imported field retained stale VMEC construction'
    nc_status = nf90_get_att(ncid, nf90_global, &
                             'sfl_allow_poor_equilibrium', overridden)
    if (nc_status /= nf90_enotatt) error stop 'imported field retained stale override'
    nc_status = nf90_get_att(ncid, nf90_global, 'regular_radial_factor', radial_factor)
    if (nc_status /= nf90_noerr .or. radial_factor /= 2) &
        error stop 'chartmap reexport lost regular radial policy'
    nc_status = nf90_close(ncid)
    if (nc_status /= nf90_noerr) error stop 'cannot close imported chartmap export'

    ! A native synthetic chartmap validates older files and rejects malformed
    ! policy attributes through the same reader used by the roundtrip.
    call write_policy_fixture()
    call read_boozer_chartmap('roundtrip_policy.nc', imported)
    if (imported%regular_radial_factor /= 1) error stop 'old chartmap policy default changed'
    call get_command_argument(0, executable)
    call write_policy_fixture(3)
    call execute_command_line('"'//trim(executable)//'" --invalid-factor'// &
                              ' > roundtrip_policy_rejection.log 2>&1', exitstat=exit_status)
    if (exit_status == 0) error stop 'invalid regular radial policy was accepted'
    call write_policy_fixture(2, floating=.true.)
    call execute_command_line('"'//trim(executable)//'" --invalid-factor'// &
                              ' > roundtrip_policy_rejection.log 2>&1', exitstat=exit_status)
    if (exit_status == 0) error stop 'floating regular radial policy was accepted'
    call write_policy_fixture(2, vector=.true.)
    call execute_command_line('"'//trim(executable)//'" --invalid-factor'// &
                              ' > roundtrip_policy_rejection.log 2>&1', exitstat=exit_status)
    if (exit_status == 0) error stop 'vector regular radial policy was accepted'

contains

    subroutine write_policy_fixture(factor, floating, vector)
        integer, intent(in), optional :: factor
        logical, intent(in), optional :: floating, vector
        integer :: file_id, dims(3), vars(11), status, j, k, v
        real(dp) :: nodes(3), angles(3), xyz(3,3,3,3)
        character(len=7), parameter :: names(11) = [character(len=7) :: &
            'rho', 'theta', 'zeta', 's', 'A_phi', 'B_theta', 'B_phi', 'x', 'y', 'z', 'Bmod']

        nodes = [0.0_dp, 0.5_dp, 1.0_dp]
        angles = TWOPI*[0.0_dp, 1.0_dp, 2.0_dp]/3.0_dp
        do k = 1, 3
            do j = 1, 3
                xyz(:,j,k,1) = 10.0_dp*cos(angles(k))
                xyz(:,j,k,2) = 10.0_dp*sin(angles(k))
                xyz(:,j,k,3) = 0.0_dp
            end do
        end do
        status = nf90_create('roundtrip_policy.nc', nf90_clobber, file_id)
        call nc_check(status)
        do v = 1, 3
            call nc_check(nf90_def_dim(file_id, trim(names(v)), 3, dims(v)))
            call nc_check(nf90_def_var(file_id, trim(names(v)), nf90_double, [dims(v)], vars(v)))
        end do
        call nc_check(nf90_def_dim(file_id, 's', 3, v))
        call nc_check(nf90_def_var(file_id, 's', nf90_double, [v], vars(4)))
        call nc_check(nf90_def_var(file_id, 'A_phi', nf90_double, [v], vars(5)))
        call nc_check(nf90_put_att(file_id, vars(5), 'radial_abscissa', 's'))
        do v = 6, 7
            call nc_check(nf90_def_var(file_id, trim(names(v)), nf90_double, [dims(1)], vars(v)))
        end do
        do v = 8, 11
            call nc_check(nf90_def_var(file_id, trim(names(v)), nf90_double, dims, vars(v)))
        end do
        call nc_check(nf90_def_var(file_id, 'num_field_periods', nf90_int, v))
        call nc_check(nf90_put_att(file_id, nf90_global, 'torflux', 1.0_dp))
        call nc_check(nf90_put_att(file_id, nf90_global, 'aminor_m', 1.0_dp))
        if (present(factor)) then
            if (present(vector)) then
                call nc_check(nf90_put_att(file_id, nf90_global, 'regular_radial_factor', [factor, factor]))
            else if (present(floating)) then
                call nc_check(nf90_put_att(file_id, nf90_global, 'regular_radial_factor', real(factor, dp)))
            else
                call nc_check(nf90_put_att(file_id, nf90_global, 'regular_radial_factor', factor))
            end if
        end if
        call nc_check(nf90_enddef(file_id))
        call nc_check(nf90_put_var(file_id, v, 1))
        call nc_check(nf90_put_var(file_id, vars(1), nodes))
        call nc_check(nf90_put_var(file_id, vars(2), angles))
        call nc_check(nf90_put_var(file_id, vars(3), angles))
        do v = 4, 7
            call nc_check(nf90_put_var(file_id, vars(v), nodes))
        end do
        do v = 8, 10
            call nc_check(nf90_put_var(file_id, vars(v), xyz(:,:,:,v-7)))
        end do
        call nc_check(nf90_put_var(file_id, vars(11), 1.0_dp+xyz(:,:,:,3)))
        call nc_check(nf90_close(file_id))
    end subroutine write_policy_fixture

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) error stop 'cannot write synthetic chartmap policy fixture'
    end subroutine nc_check

    !> Interior test points: s in (0.1, 0.8), full theta period, scattered phi.
    subroutine init_test_points(phi_per)
        real(dp), intent(in) :: phi_per
        integer :: j
        real(dp) :: frac

        do j = 1, n_test
            frac = real(j, dp)/real(n_test + 1, dp)
            s_test(j) = 0.1_dp + 0.7_dp*frac
            th_test(j) = TWOPI*frac
            ph_test(j) = phi_per*mod(real(2*j + 1, dp), real(n_test + 1, dp)) &
                         /real(n_test + 1, dp)
        end do
    end subroutine init_test_points

    !> Issue #395 probe sequence on the exported production chartmap: an exact
    !> forward image at interior rho (including rho = 0.5 and toroidal angles
    !> outside the stored wedge) must be recovered by from_cyl. The result is
    !> stored so the after-load check can pin invariance against the
    !> module-level Boozer spline construction.
    subroutine check_coordinate_inverse_before_load(path)
        character(len=*), intent(in) :: path

        real(dp) :: x(3), xcyl(3), u_back(3), x_back(3)
        integer :: j, ierr

        u_ref_chk(:, 1) = [0.5_dp, 1.9_dp, 0.35_dp*TWOPI]
        u_ref_chk(:, 2) = [0.5_dp, 4.0_dp, 0.80_dp*TWOPI]
        u_ref_chk(:, 3) = [0.3_dp, 0.7_dp, 0.55_dp*TWOPI]

        call make_chartmap_coordinate_system(coord_cs, path)

        do j = 1, n_coord_chk
            call coord_cs%evaluate_cart(u_ref_chk(:, j), x)
            call cart_to_cyl(x, xcyl)
            call coord_cs%from_cyl(xcyl, u_back, ierr)
            if (ierr /= chartmap_from_cyl_ok) then
                print *, "coordinate inverse failed, ierr=", ierr, &
                    " at u_ref=", u_ref_chk(:, j)
                error stop "chartmap coordinate inverse failed on exported chart"
            end if
            call coord_cs%evaluate_cart(u_back, x_back)
            if (sqrt(sum((x_back - x)**2)) > 1.0e-6_dp) then
                print *, "roundtrip residual=", sqrt(sum((x_back - x)**2)), &
                    " at u_ref=", u_ref_chk(:, j)
                error stop "chartmap coordinate roundtrip residual too large"
            end if
            if (abs(u_back(1) - u_ref_chk(1, j)) > 1.0e-6_dp) then
                print *, "recovered rho=", u_back(1), " want=", u_ref_chk(1, j)
                error stop "chartmap coordinate inverse lost the radial label"
            end if
            u_back_before(:, j) = u_back
        end do
        print *, "Coordinate inverse recovers exported-chart points before load"
    end subroutine check_coordinate_inverse_before_load

    !> Building the module-level Boozer field splines must not change the
    !> coordinate inverse of the standalone chartmap coordinate system.
    subroutine check_coordinate_inverse_after_load()
        real(dp) :: x(3), xcyl(3), u_back(3)
        integer :: j, ierr

        do j = 1, n_coord_chk
            call coord_cs%evaluate_cart(u_ref_chk(:, j), x)
            call cart_to_cyl(x, xcyl)
            call coord_cs%from_cyl(xcyl, u_back, ierr)
            if (ierr /= chartmap_from_cyl_ok) then
                error stop "coordinate inverse regressed after Boozer spline load"
            end if
            if (any(u_back /= u_back_before(:, j))) then
                print *, "u before=", u_back_before(:, j)
                print *, "u after =", u_back
                error stop "Boozer spline load changed the coordinate inverse"
            end if
        end do
        print *, "Boozer spline load leaves the coordinate inverse unchanged"
    end subroutine check_coordinate_inverse_after_load

end program test_boozer_chartmap_roundtrip
