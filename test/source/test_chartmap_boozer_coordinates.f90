program test_chartmap_boozer_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, validate_chartmap_file, &
                                  chartmap_from_cyl_ok, BOOZER
    use math_constants, only: TWOPI
    use netcdf
    implicit none

    character(len=*), parameter :: filename = "chartmap_boozer.nc"
    character(len=*), parameter :: filename_reactor = "chartmap_boozer_reactor.nc"
    integer :: ierr, nerrors
    character(len=2048) :: message

    nerrors = 0

    call write_boozer_file(filename)
    call validate_chartmap_file(filename, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validator rejected Boozer chartmap: ", trim(message)
        error stop 1
    end if

    call run_roundtrip_checks(filename, nerrors)

    call write_reactor_boozer_file(filename_reactor)
    call run_full_torus_cold_inverse_checks(filename_reactor, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in Boozer chartmap tests"
        error stop 1
    end if

    print *, "All Boozer chartmap coordinate tests passed!"

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    subroutine run_roundtrip_checks(path, nerrors)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        real(dp) :: u(3), u_back(3), u_cyl(3), x(3), xcyl(3)
        real(dp) :: du_theta, du_zeta, max_dphi
        real(dp) :: zeta_period, dphi_geom
        integer :: ir, it, iz, ierr_local
        real(dp), parameter :: rho_vals(3) = [0.12_dp, 0.45_dp, 0.82_dp]
        real(dp), parameter :: theta_vals(4) = [0.1_dp, 1.4_dp, 3.0_dp, 4.9_dp]
        real(dp), parameter :: zeta_fracs(4) = [0.0_dp, 0.17_dp, 0.41_dp, 0.73_dp]
        real(dp), parameter :: tol_u = 5.0e-6_dp

        max_dphi = 0.0_dp
        call make_chartmap_coordinate_system(cs, path)

        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            if (ccs%zeta_convention /= BOOZER) then
                print *, "  FAIL: expected BOOZER zeta convention, got ", &
                    ccs%zeta_convention
                nerrors = nerrors + 1
                return
            end if

            zeta_period = TWOPI/real(ccs%num_field_periods, dp)

            do iz = 1, size(zeta_fracs)
                do it = 1, size(theta_vals)
                    do ir = 1, size(rho_vals)
                        u = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*zeta_period]
                        call ccs%evaluate_cart(u, x)

                        xcyl(1) = sqrt(x(1)**2 + x(2)**2)
                        xcyl(2) = atan2(x(2), x(1))
                        xcyl(3) = x(3)

                        dphi_geom = angular_diff(xcyl(2), u(3), TWOPI)
                        max_dphi = max(max_dphi, dphi_geom)

                        call ccs%from_cart(x, u_back, ierr_local)
                        if (ierr_local /= 0) then
                            print *, "  FAIL: from_cart ierr=", ierr_local, " at u=", u
                            nerrors = nerrors + 1
                            cycle
                        end if

                        du_theta = angular_diff(u_back(2), u(2), TWOPI)
                        du_zeta = angular_diff(u_back(3), u(3), zeta_period)
                        if (abs(u_back(1) - u(1)) > tol_u .or. du_theta > tol_u .or. &
                            du_zeta > tol_u) then
                            print *, "  FAIL: from_cart roundtrip mismatch"
                            print *, "    u     =", u
                            print *, "    u_back=", u_back
                            nerrors = nerrors + 1
                            cycle
                        end if

                        call ccs%from_cyl(xcyl, u_cyl, ierr_local)
                        if (ierr_local /= chartmap_from_cyl_ok) then
                            print *, "  FAIL: from_cyl ierr=", ierr_local, " at u=", u
                            nerrors = nerrors + 1
                            cycle
                        end if

                        du_theta = angular_diff(u_cyl(2), u(2), TWOPI)
                        du_zeta = angular_diff(u_cyl(3), u(3), zeta_period)
                        if (abs(u_cyl(1) - u(1)) > tol_u .or. du_theta > tol_u .or. &
                            du_zeta > tol_u) then
                            print *, "  FAIL: from_cyl roundtrip mismatch"
                            print *, "    u     =", u
                            print *, "    u_cyl =", u_cyl
                            nerrors = nerrors + 1
                        end if
                    end do
                end do
            end do

            if (max_dphi < 5.0e-2_dp) then
                print *, "  FAIL: synthetic Boozer chartmap is too close to cylindrical"
                nerrors = nerrors + 1
            else
                print *, "  PASS: Boozer chartmap uses non-cylindrical toroidal parameter"
            end if
        class default
            print *, "  FAIL: chartmap file did not load as chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_roundtrip_checks

    pure real(dp) function angular_diff(a, b, period)
        real(dp), intent(in) :: a, b, period
        angular_diff = abs(modulo(a - b + 0.5_dp*period, period) - 0.5_dp*period)
    end function angular_diff

    ! Regression for issue #395: the cold inverse of a rotating-frame Boozer chart
    ! must recover interior points whose toroidal angle lies anywhere on the full
    ! torus, not only inside the stored field-period wedge. The lab-frame position
    ! of the co-rotating chart is 2pi-periodic in zeta, so zeta is compared modulo
    ! TWOPI and the recovered u must reproduce the Cartesian target.
    subroutine run_full_torus_cold_inverse_checks(path, nerrors)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs

        call make_chartmap_coordinate_system(cs, path)

        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            call check_cold_interior_points(ccs, nerrors)
            call check_cold_outside_target(ccs, nerrors)
            call check_cold_near_axis_target(ccs, nerrors)
        class default
            print *, "  FAIL: reactor chartmap did not load as chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_full_torus_cold_inverse_checks

    subroutine check_cold_interior_points(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), u_back(3), x(3), xcyl(3), x_back(3)
        real(dp) :: du_theta, du_zeta
        integer :: ir, it, iz, ierr_local, npoints, nfail
        real(dp), parameter :: rho_vals(4) = [0.12_dp, 0.45_dp, 0.5_dp, 0.82_dp]
        real(dp), parameter :: theta_vals(5) = [0.0_dp, 1.4_dp, 3.0_dp, 4.9_dp, &
                                                6.2_dp]
        real(dp), parameter :: zeta_fracs(6) = [0.0_dp, 0.13_dp, 0.29_dp, 0.48_dp, &
                                                0.77_dp, 0.999_dp]
        real(dp), parameter :: tol_u = 1.0e-8_dp
        real(dp), parameter :: tol_x = 1.0e-8_dp

        npoints = 0
        nfail = 0
        do iz = 1, size(zeta_fracs)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*TWOPI]
                    ! Every third point enters many field periods outside the
                    ! stored wedge; the inverse must recover zeta modulo TWOPI.
                    if (mod(npoints, 3) == 0) u(3) = u(3) + 3.0_dp*TWOPI
                    npoints = npoints + 1

                    call ccs%evaluate_cart(u, x)
                    xcyl(1) = sqrt(x(1)**2 + x(2)**2)
                    xcyl(2) = atan2(x(2), x(1))
                    xcyl(3) = x(3)

                    call ccs%from_cyl(xcyl, u_back, ierr_local)
                    if (ierr_local /= chartmap_from_cyl_ok) then
                        print *, "  FAIL: cold from_cyl ierr=", ierr_local, &
                            " at u=", u
                        nfail = nfail + 1
                        cycle
                    end if

                    du_theta = angular_diff(u_back(2), u(2), TWOPI)
                    du_zeta = angular_diff(u_back(3), u(3), TWOPI)
                    call ccs%evaluate_cart(u_back, x_back)
                    if (abs(u_back(1) - u(1)) > tol_u .or. du_theta > tol_u .or. &
                        du_zeta > tol_u .or. &
                        sqrt(sum((x_back - x)**2)) > tol_x) then
                        print *, "  FAIL: cold from_cyl roundtrip mismatch"
                        print *, "    u     =", u
                        print *, "    u_back=", u_back
                        print *, "    |dx|  =", sqrt(sum((x_back - x)**2))
                        nfail = nfail + 1
                    end if
                end do
            end do
        end do

        nerrors = nerrors + nfail
        if (nfail == 0) then
            print *, "  PASS: ", npoints, " cold full-torus roundtrips recovered"
        end if
    end subroutine check_cold_interior_points

    ! A target genuinely outside the mapped volume must return a nonzero status
    ! instead of masquerading as an interior root.
    subroutine check_cold_outside_target(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u_edge(3), x_edge(3), x(3), u_back(3)
        real(dp) :: e_cov(3, 3)
        integer :: ierr_local

        u_edge = [1.0_dp, 0.9_dp, 0.7_dp*TWOPI]
        call ccs%evaluate_cart(u_edge, x_edge)
        call ccs%covariant_basis(u_edge, e_cov)
        x = x_edge + 0.5_dp*e_cov(:, 1)

        call ccs%from_cart(x, u_back, ierr_local)
        if (ierr_local == chartmap_from_cyl_ok) then
            print *, "  FAIL: past-edge target accepted as interior root, u=", u_back
            nerrors = nerrors + 1
        else
            print *, "  PASS: past-edge target returns nonzero status"
        end if
    end subroutine check_cold_outside_target

    ! A near-axis target below the chart resolution must either locate with a tiny
    ! radius that reproduces the target, or fail explicitly; it must never return
    ! success on an unrelated radial branch.
    subroutine check_cold_near_axis_target(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), u_back(3), x(3), x_back(3)
        integer :: ierr_local

        u = [1.0e-7_dp, 2.3_dp, 0.6_dp*TWOPI]
        call ccs%evaluate_cart(u, x)
        call ccs%from_cart(x, u_back, ierr_local)
        if (ierr_local == chartmap_from_cyl_ok) then
            call ccs%evaluate_cart(u_back, x_back)
            if (u_back(1) > 1.0e-4_dp .or. sqrt(sum((x_back - x)**2)) > 1.0e-6_dp) then
                print *, "  FAIL: near-axis target accepted on wrong branch"
                print *, "    u_back=", u_back
                nerrors = nerrors + 1
                return
            end if
        end if
        print *, "  PASS: near-axis target not accepted on an unrelated branch"
    end subroutine check_cold_near_axis_target

    ! Reactor-scale endpoint-excluded rotating-frame Boozer chart, nfp = 5, with a
    ! field-period-symmetric Boozer toroidal shift, sized so the inverse exercises
    ! the production probe/refine/solve path at reactor dimensions.
    subroutine write_reactor_boozer_file(path)
        character(len=*), intent(in) :: path

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_nfp
        integer, parameter :: nrho = 25, ntheta = 48, nzeta = 32, nfp = 5
        real(dp), parameter :: r0 = 1000.0_dp
        real(dp), parameter :: a = 250.0_dp
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
                    z(nrho, ntheta, nzeta)
        real(dp) :: rho_val, theta_val, zeta_val
        real(dp) :: phi_geom, rmaj, zpos, period
        integer :: ir, it, iz

        period = TWOPI/real(nfp, dp)

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
        end do
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = period*real(iz - 1, dp)/real(nzeta, dp)
        end do

        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    phi_geom = zeta_val + 0.12_dp*rho_val*sin(theta_val) + &
                               0.03_dp*sin(real(nfp, dp)*zeta_val)
                    rmaj = r0 + a*rho_val*cos(theta_val) + &
                           0.05_dp*a*rho_val*cos(theta_val - real(nfp, dp)*zeta_val)
                    zpos = a*rho_val*sin(theta_val) + &
                           0.03_dp*a*rho_val*sin(real(nfp, dp)*zeta_val)

                    x(ir, it, iz) = rmaj*cos(phi_geom)
                    y(ir, it, iz) = rmaj*sin(phi_geom)
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

        call nc_check(nf90_create(trim(path), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "boozer"))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "rho_convention", "rho_tor"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, var_nfp))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_nfp, nfp))
        call nc_check(nf90_close(ncid))
    end subroutine write_reactor_boozer_file

    subroutine write_boozer_file(path)
        character(len=*), intent(in) :: path

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_nfp
        integer, parameter :: nrho = 9, ntheta = 11, nzeta = 10, nfp = 3
        real(dp), parameter :: r0 = 180.0_dp
        real(dp), parameter :: a = 45.0_dp
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)
        real(dp) :: rho_val, theta_val, zeta_val
        real(dp) :: phi_geom, rmaj, zpos, period
        integer :: ir, it, iz

        period = TWOPI/real(nfp, dp)

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
        end do
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = period*real(iz - 1, dp)/real(nzeta, dp)
        end do

        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    phi_geom = zeta_val + 0.08_dp*rho_val*sin(theta_val) + &
                               0.02_dp*sin(2.0_dp*zeta_val)
                    rmaj = r0 + a*rho_val*cos(theta_val) + &
                           0.04_dp*a*rho_val*cos(theta_val - 2.0_dp*zeta_val)
                    zpos = a*rho_val*sin(theta_val) + &
                           0.02_dp*a*rho_val*sin(zeta_val)

                    x(ir, it, iz) = rmaj*cos(phi_geom)
                    y(ir, it, iz) = rmaj*sin(phi_geom)
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

        call nc_check(nf90_create(trim(path), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "boozer"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, var_nfp))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_nfp, nfp))
        call nc_check(nf90_close(ncid))
    end subroutine write_boozer_file

end program test_chartmap_boozer_coordinates
