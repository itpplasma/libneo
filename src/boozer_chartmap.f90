module boozer_chartmap
    !> File I/O for extended Boozer chartmaps. Keeps the NetCDF read/write off
    !> the boozer_sub converter: load parses a chartmap and hands the record to
    !> boozer_sub's build_boozer_from_chartmap; export gathers via the public
    !> splint entry points and writes the NetCDF.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_sub, only: build_boozer_from_chartmap, splint_boozer_coord, &
                          boozer_to_vmec
    use boozer_chartmap_io, only: read_boozer_chartmap
    use boozer_chartmap_types, only: boozer_chartmap_data_t

    implicit none
    private
    public :: load_boozer_from_chartmap, export_boozer_chartmap

contains

    subroutine load_boozer_from_chartmap(filename)
        !> Parse an extended chartmap NetCDF file and build the module-level
        !> Boozer batch splines from it.
        character(len=*), intent(in) :: filename

        type(boozer_chartmap_data_t) :: d

        call read_boozer_chartmap(filename, d)
        call build_boozer_from_chartmap(d)
        print *, 'Loaded Boozer splines from chartmap: ', trim(filename)
    end subroutine load_boozer_from_chartmap

    subroutine export_boozer_chartmap(filename)
        !> Export Boozer coordinate data computed by get_boozer_coordinates()
        !> to an extended chartmap NetCDF file. Must be called after
        !> get_boozer_coordinates() and while VMEC splines are still active
        !> (needed for X, Y, Z geometry evaluation).
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper, raxis_cc, zaxis_cs, ntor_axis
        use boozer_coordinates_mod, only: ns_B, n_theta_B, n_phi_B, &
                                          h_theta_B, h_phi_B
        use spline_vmec_sub, only: splint_vmec_data
        use netcdf

        character(len=*), intent(in) :: filename

        integer :: ncid, status
        integer :: dim_rho, dim_s, dim_theta, dim_zeta
        integer :: var_rho, var_s, var_theta, var_zeta
        integer :: var_x, var_y, var_z
        integer :: var_aphi, var_btheta, var_bphi, var_bmod, var_nfp
        integer :: i_rho, i_theta, i_phi
        integer :: n_theta_out, n_phi_out
        real(dp) :: s, theta_B, phi_B, theta_V, phi_V
        real(dp) :: R, Zval, alam
        real(dp) :: A_phi_dum, A_theta_dum, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: d2A_phi_ds2, d3A_phi_ds3, B_theta_val, B_phi_val
        real(dp) :: dB_theta, d2B_theta, dB_phi, d2B_phi, Bmod_val, B_r_val
        real(dp) :: dBmod(3), d2Bmod(6), dB_r(3), d2B_r(6)
        real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
        real(dp) :: dl_ds, dl_dt, dl_dp
        real(dp) :: R_axis, Z_axis, phi_V0, axis_geom_err, bmod_axis
        ! Chartmap now reaches the exact magnetic axis (rho=0, s=0). The axis is
        ! the m=0 limit: theta-independent geometry from the exact VMEC axis
        ! curve raxis_cc/zaxis_cs (no flux surface exists at s=0). boozer_to_vmec
        ! and splint_boozer_coord are singular at s=0, so the toroidal angle and
        ! the field profiles are taken at the regular limit s_axis_eval (the
        ! Boozer-data spline floor), while the axis position is exact.
        real(dp), parameter :: rho_min = 0.0_dp
        real(dp), parameter :: s_axis_eval = 1.0e-6_dp
        real(dp), allocatable :: rho_arr(:), s_arr(:), theta_arr(:), zeta_arr(:)
        real(dp), allocatable :: A_phi_arr(:), B_theta_arr(:), B_phi_arr(:)
        real(dp), allocatable :: x_arr(:, :, :), y_arr(:, :, :), z_arr(:, :, :)
        real(dp), allocatable :: bmod_arr(:, :, :)

        ! Chartmap files store endpoint-excluded angular grids. The reader
        ! reconstructs endpoint planes before building periodic splines.
        n_theta_out = n_theta_B - 1
        n_phi_out = n_phi_B - 1

        allocate (rho_arr(ns_B))
        allocate (s_arr(ns_B))
        allocate (theta_arr(n_theta_out), zeta_arr(n_phi_out))
        allocate (A_phi_arr(ns_B), B_theta_arr(ns_B), B_phi_arr(ns_B))
        allocate (x_arr(ns_B, n_theta_out, n_phi_out))
        allocate (y_arr(ns_B, n_theta_out, n_phi_out))
        allocate (z_arr(ns_B, n_theta_out, n_phi_out))
        allocate (bmod_arr(ns_B, n_theta_out, n_phi_out))

        ! Radial grid
        do i_rho = 1, ns_B
            rho_arr(i_rho) = rho_min + (1.0_dp - rho_min) &
                             *real(i_rho - 1, dp)/real(ns_B - 1, dp)
            s_arr(i_rho) = rho_min**2 + (1.0_dp - rho_min**2) &
                           *real(i_rho - 1, dp)/real(ns_B - 1, dp)
        end do
        ! Angular grids (endpoint excluded, for chartmap geometry)
        do i_theta = 1, n_theta_out
            theta_arr(i_theta) = real(i_theta - 1, dp)*h_theta_B
        end do
        do i_phi = 1, n_phi_out
            zeta_arr(i_phi) = real(i_phi - 1, dp)*h_phi_B
        end do

        ! A_phi is a flux profile on s. B_theta/B_phi stay on rho for now.
        do i_rho = 1, ns_B
            call splint_boozer_coord(max(s_arr(i_rho), s_axis_eval), 0.0_dp, 0.0_dp, 0, &
                                     A_theta_dum, A_phi_arr(i_rho), dA_theta_ds, &
                                     dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3, &
                                     B_theta_val, dB_theta, d2B_theta, &
                                     B_phi_val, dB_phi, d2B_phi, &
                                     Bmod_val, dBmod, d2Bmod, &
                                     B_r_val, dB_r, d2B_r)
            s = max(rho_arr(i_rho)**2, s_axis_eval)
            call splint_boozer_coord(s, 0.0_dp, 0.0_dp, 0, &
                                     A_theta_dum, A_phi_dum, dA_theta_ds, &
                                     dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3, &
                                     B_theta_arr(i_rho), dB_theta, d2B_theta, &
                                     B_phi_arr(i_rho), dB_phi, d2B_phi, &
                                     Bmod_val, dBmod, d2Bmod, &
                                     B_r_val, dB_r, d2B_r)
        end do

        ! Axis slice (rho=0): the exact magnetic axis is m=0, theta-independent.
        ! boozer_to_vmec degenerates poloidally at s=0, so take the regular
        ! toroidal limit phi_V0(phi_B) at s_axis_eval and average over theta to
        ! cancel the O(rho) poloidal contamination (the m=1 mean is zero on the
        ! full-period grid), leaving phi_V0 to O(rho**2). The axis position is
        ! the exact VMEC axis curve, identical for every theta -> regular at rho=0.
        axis_geom_err = 0.0_dp
        do i_phi = 1, n_phi_out
            phi_B = zeta_arr(i_phi)
            phi_V0 = 0.0_dp
            do i_theta = 1, n_theta_out
                call boozer_to_vmec(s_axis_eval, theta_arr(i_theta), phi_B, theta_V, phi_V)
                phi_V0 = phi_V0 + phi_V
            end do
            phi_V0 = phi_V0/real(n_theta_out, dp)
            call eval_vmec_axis(phi_V0, R_axis, Z_axis)
            x_arr(1, :, i_phi) = R_axis*cos(phi_V0)
            y_arr(1, :, i_phi) = R_axis*sin(phi_V0)
            z_arr(1, :, i_phi) = Z_axis

            ! Convention + smoothness gate: at s=0 splint_vmec_data returns the
            ! m=0 (theta-independent) limit; with the m=0 healing anchored to
            ! the exact axis it must equal the exact axis.
            call splint_vmec_data(0.0_dp, 0.0_dp, phi_V0, &
                                  A_phi_dum, A_theta_dum, dA_phi_ds, &
                                  dA_theta_ds, aiota, &
                                  R, Zval, alam, dR_ds, dR_dt, dR_dp, &
                                  dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
            axis_geom_err = max(axis_geom_err, abs(R - R_axis) + abs(Zval - Z_axis))
        end do

        ! Interior (rho>0): Boozer angles -> VMEC angles -> VMEC geometry.
        do i_phi = 1, n_phi_out
            do i_theta = 1, n_theta_out
                do i_rho = 2, ns_B
                    s = rho_arr(i_rho)**2
                    theta_B = theta_arr(i_theta)
                    phi_B = zeta_arr(i_phi)
                    call boozer_to_vmec(s, theta_B, phi_B, theta_V, phi_V)
                    call splint_vmec_data(s, theta_V, phi_V, &
                                          A_phi_dum, A_theta_dum, dA_phi_ds, &
                                          dA_theta_ds, aiota, &
                                          R, Zval, alam, dR_ds, dR_dt, dR_dp, &
                                          dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
                    x_arr(i_rho, i_theta, i_phi) = R*cos(phi_V)
                    y_arr(i_rho, i_theta, i_phi) = R*sin(phi_V)
                    z_arr(i_rho, i_theta, i_phi) = Zval
                end do
            end do
        end do
        print *, 'export_boozer_chartmap: axis vs VMEC m=0 limit max|dR|+|dZ| =', &
            axis_geom_err

        ! Bmod axis slice (rho=0): theta-independent m=0 axis field strength,
        ! the theta-average at s_axis_eval (the m=1 mean is zero on the grid).
        do i_phi = 1, n_phi_out
            phi_B = zeta_arr(i_phi)
            bmod_axis = 0.0_dp
            do i_theta = 1, n_theta_out
                call splint_boozer_coord(s_axis_eval, theta_arr(i_theta), phi_B, 0, &
                                         A_theta_dum, A_phi_dum, dA_theta_ds, &
                                         dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3, &
                                         B_theta_val, dB_theta, d2B_theta, &
                                         B_phi_val, dB_phi, d2B_phi, &
                                         Bmod_val, dBmod, d2Bmod, &
                                         B_r_val, dB_r, d2B_r)
                bmod_axis = bmod_axis + Bmod_val
            end do
            bmod_arr(1, :, i_phi) = bmod_axis/real(n_theta_out, dp)
        end do

        ! Bmod interior (rho>0)
        do i_phi = 1, n_phi_out
            phi_B = real(i_phi - 1, dp)*h_phi_B
            do i_theta = 1, n_theta_out
                theta_B = real(i_theta - 1, dp)*h_theta_B
                do i_rho = 2, ns_B
                    s = rho_arr(i_rho)**2
                    call splint_boozer_coord(s, theta_B, phi_B, 0, &
                                             A_theta_dum, A_phi_dum, dA_theta_ds, &
                                             dA_phi_ds, d2A_phi_ds2, d3A_phi_ds3, &
                                             B_theta_val, dB_theta, d2B_theta, &
                                             B_phi_val, dB_phi, d2B_phi, &
                                             bmod_arr(i_rho, i_theta, i_phi), &
                                             dBmod, d2Bmod, &
                                             B_r_val, dB_r, d2B_r)
                end do
            end do
        end do

        ! Write NetCDF file
        status = nf90_create(trim(filename), nf90_clobber, ncid)
        call nc_assert(status, "create "//trim(filename))

        ! Dimensions: one endpoint-excluded angular grid for geometry and fields.
        call nc_assert(nf90_def_dim(ncid, "rho", ns_B, dim_rho), "def_dim rho")
        call nc_assert(nf90_def_dim(ncid, "s", ns_B, dim_s), "def_dim s")
        call nc_assert(nf90_def_dim(ncid, "theta", n_theta_out, dim_theta), &
                       "def_dim theta")
        call nc_assert(nf90_def_dim(ncid, "zeta", n_phi_out, dim_zeta), &
                       "def_dim zeta")

        ! Coordinate variables
        call nc_assert(nf90_def_var(ncid, "rho", nf90_double, [dim_rho], var_rho), &
                       "def_var rho")
        call nc_assert(nf90_def_var(ncid, "s", nf90_double, [dim_s], var_s), &
                       "def_var s")
        call nc_assert(nf90_def_var(ncid, "theta", nf90_double, [dim_theta], &
                                    var_theta), "def_var theta")
        call nc_assert(nf90_def_var(ncid, "zeta", nf90_double, [dim_zeta], &
                                    var_zeta), "def_var zeta")

        ! Geometry (NF90 reverses dims: Fortran (rho,theta,zeta) -> NetCDF order)
        call nc_assert(nf90_def_var(ncid, "x", nf90_double, &
                                    [dim_rho, dim_theta, dim_zeta], var_x), &
                       "def_var x")
        call nc_assert(nf90_put_att(ncid, var_x, "units", "cm"), "att x units")
        call nc_assert(nf90_def_var(ncid, "y", nf90_double, &
                                    [dim_rho, dim_theta, dim_zeta], var_y), &
                       "def_var y")
        call nc_assert(nf90_put_att(ncid, var_y, "units", "cm"), "att y units")
        call nc_assert(nf90_def_var(ncid, "z", nf90_double, &
                                    [dim_rho, dim_theta, dim_zeta], var_z), &
                       "def_var z")
        call nc_assert(nf90_put_att(ncid, var_z, "units", "cm"), "att z units")

        ! Boozer field data
        call nc_assert(nf90_def_var(ncid, "A_phi", nf90_double, [dim_s], &
                                    var_aphi), "def_var A_phi")
        call nc_assert(nf90_put_att(ncid, var_aphi, "radial_abscissa", "s"), &
                       "att A_phi radial_abscissa")
        call nc_assert(nf90_def_var(ncid, "B_theta", nf90_double, [dim_rho], &
                                    var_btheta), "def_var B_theta")
        call nc_assert(nf90_def_var(ncid, "B_phi", nf90_double, [dim_rho], &
                                    var_bphi), "def_var B_phi")
        call nc_assert(nf90_def_var(ncid, "Bmod", nf90_double, &
                                    [dim_rho, dim_theta, dim_zeta], var_bmod), &
                       "def_var Bmod")
        call nc_assert(nf90_def_var(ncid, "num_field_periods", nf90_int, var_nfp), &
                       "def_var nfp")

        ! Global attributes
        call nc_assert(nf90_put_att(ncid, nf90_global, "rho_convention", "rho_tor"), &
                       "att rho_convention")
        call nc_assert(nf90_put_att(ncid, nf90_global, "zeta_convention", "boozer"), &
                       "att zeta_convention")
        call nc_assert(nf90_put_att(ncid, nf90_global, "rho_lcfs", rho_arr(ns_B)), &
                       "att rho_lcfs")
        call nc_assert(nf90_put_att(ncid, nf90_global, "boozer_field", 1), &
                       "att boozer_field")
        call nc_assert(nf90_put_att(ncid, nf90_global, "torflux", torflux), &
                       "att torflux")
        ! No rmajor attribute: the chartmap reader derives the major radius
        ! from the innermost-surface geometry (see boozer_chartmap_io).

        call nc_assert(nf90_enddef(ncid), "enddef")

        ! Write data
        call nc_assert(nf90_put_var(ncid, var_rho, rho_arr), "put rho")
        call nc_assert(nf90_put_var(ncid, var_s, s_arr), "put s")
        call nc_assert(nf90_put_var(ncid, var_theta, theta_arr), "put theta")
        call nc_assert(nf90_put_var(ncid, var_zeta, zeta_arr), "put zeta")
        call nc_assert(nf90_put_var(ncid, var_x, x_arr), "put x")
        call nc_assert(nf90_put_var(ncid, var_y, y_arr), "put y")
        call nc_assert(nf90_put_var(ncid, var_z, z_arr), "put z")
        call nc_assert(nf90_put_var(ncid, var_aphi, A_phi_arr), "put A_phi")
        call nc_assert(nf90_put_var(ncid, var_btheta, B_theta_arr), "put B_theta")
        call nc_assert(nf90_put_var(ncid, var_bphi, B_phi_arr), "put B_phi")
        call nc_assert(nf90_put_var(ncid, var_bmod, bmod_arr), "put Bmod")
        call nc_assert(nf90_put_var(ncid, var_nfp, nper), "put nfp")

        call nc_assert(nf90_close(ncid), "close")

        print *, 'Exported Boozer chartmap to ', trim(filename)
        print *, '  nfp=', nper, ' ns=', ns_B, ' ntheta=', n_theta_out, &
            ' nphi=', n_phi_out
        print *, '  torflux=', torflux

    contains

        subroutine nc_assert(stat, loc)
            integer, intent(in) :: stat
            character(len=*), intent(in) :: loc
            if (stat /= nf90_noerr) then
                print *, "export_boozer_chartmap: NetCDF error at ", trim(loc), &
                    ": ", trim(nf90_strerror(stat))
                error stop
            end if
        end subroutine nc_assert

        !> Exact VMEC magnetic axis (m=0) at the VMEC toroidal angle phi_V.
        !> Stellarator-symmetric: R = sum_n raxis_cc(n) cos(n nfp phi_V),
        !> Z = sum_n zaxis_cs(n) sin(n nfp phi_V). Sign convention is checked
        !> against the s=0 limit of splint_vmec_data (axis_geom_err).
        subroutine eval_vmec_axis(phi_V, R_axis_out, Z_axis_out)
            real(dp), intent(in) :: phi_V
            real(dp), intent(out) :: R_axis_out, Z_axis_out
            integer :: nm
            real(dp) :: ang
            R_axis_out = 0.0_dp
            Z_axis_out = 0.0_dp
            ! VMEC convention R = sum rmnc cos(m u - n v), Z = sum zmns sin(m u - n v).
            ! At m=0: cos(-n v) = cos(n v), sin(-n v) = -sin(n v), hence the Z sign.
            do nm = 0, ntor_axis
                ang = real(nm, dp)*real(nper, dp)*phi_V
                R_axis_out = R_axis_out + raxis_cc(nm)*cos(ang)
                Z_axis_out = Z_axis_out - zaxis_cs(nm)*sin(ang)
            end do
        end subroutine eval_vmec_axis

    end subroutine export_boozer_chartmap

end module boozer_chartmap
