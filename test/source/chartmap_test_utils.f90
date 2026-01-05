module chartmap_test_utils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: make_chartmap_coordinate_system, &
                                  coordinate_system_t, chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    use new_vmec_stuff_mod, only: netcdffile
    use spline_vmec_sub, only: spline_vmec_data
    use libneo_coordinates, only: vmec_coordinate_system_t
    use nctools_module, only: nc_open, nc_close, nc_get, nc_inq_dim
    implicit none

    private
    public :: chartmap_roundtrip_check
    public :: chartmap_boundary_matches_vmec_check
    public :: chartmap_boundary_offset_check

contains

    subroutine chartmap_roundtrip_check(filename, tol_x, nerrors)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: tol_x
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        integer :: ierr_local
        real(dp) :: zeta_period

        call make_chartmap_coordinate_system(cs, filename)
        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            zeta_period = TWOPI/real(ccs%num_field_periods, dp)
            u = [0.55_dp, 0.25_dp*TWOPI, 0.4_dp*zeta_period]
            call ccs%evaluate_cart(u, x)
            xcyl(1) = sqrt(x(1)**2 + x(2)**2)
            xcyl(2) = atan2(x(2), x(1))
            xcyl(3) = x(3)

            call ccs%from_cyl(xcyl, u_back, ierr_local)
            if (ierr_local /= 0) then
                print *, "  FAIL: from_cyl ierr=", ierr_local
                nerrors = nerrors + 1
                return
            end if

            call ccs%evaluate_cart(u_back, x_round)
            if (maxval(abs(x_round - x)) > tol_x) then
                print *, "  FAIL: x roundtrip mismatch max|dx|=", &
                    maxval(abs(x_round - x))
                nerrors = nerrors + 1
                return
            end if

            print *, "  PASS: generated chartmap roundtrip u->x->u_back"
        class default
            print *, &
                "  FAIL: make_chartmap_coordinate_system did not return chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine chartmap_roundtrip_check

    subroutine chartmap_boundary_matches_vmec_check(wout_file, chartmap_file, tol_r, tol_z, nerrors)
        character(len=*), intent(in) :: wout_file
        character(len=*), intent(in) :: chartmap_file
        real(dp), intent(in) :: tol_r
        real(dp), intent(in) :: tol_z
        integer, intent(inout) :: nerrors

        type(vmec_coordinate_system_t) :: vmec
        integer :: ncid
        integer :: nrho, ntheta, nzeta
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        real(dp) :: u_vmec(3), x_vmec_cyl(3)
        real(dp) :: r_chart, z_chart
        real(dp) :: dr, dz
        real(dp) :: max_dr, max_dz
        integer :: ir, it, iz
        integer :: iz_list(5)
        integer :: nz, n_used
        integer :: k

        max_dr = 0.0_dp
        max_dz = 0.0_dp

        netcdffile = wout_file
        call spline_vmec_data

        call nc_open(chartmap_file, ncid)
        call nc_inq_dim(ncid, "rho", nrho)
        call nc_inq_dim(ncid, "theta", ntheta)
        call nc_inq_dim(ncid, "zeta", nzeta)

        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        call nc_get(ncid, "rho", rho)
        call nc_get(ncid, "theta", theta)
        call nc_get(ncid, "zeta", zeta)

        ! NetCDF stores x/y/z in file order (zeta, theta, rho).
        ! Fortran is column-major; read into reversed dimension order so that
        ! indexing is natural: x(ir, it, iz) with (rho, theta, zeta).
        allocate (x(nrho, ntheta, nzeta))
        allocate (y(nrho, ntheta, nzeta))
        allocate (z(nrho, ntheta, nzeta))
        call nc_get(ncid, "x", x)
        call nc_get(ncid, "y", y)
        call nc_get(ncid, "z", z)
        call nc_close(ncid)

        nz = size(zeta)
        iz_list = [1, 2, max(1, nz/2), max(1, nz - 1), nz]
        n_used = 0
        do k = 1, size(iz_list)
            iz = iz_list(k)
            if (iz < 1 .or. iz > nz) cycle
            if (k > 1) then
                if (any(iz_list(1:k-1) == iz)) cycle
            end if
            n_used = n_used + 1

            ir = size(rho)
            do it = 1, size(theta)
                r_chart = sqrt(x(ir, it, iz)**2 + y(ir, it, iz)**2)
                z_chart = z(ir, it, iz)

                u_vmec = [1.0_dp, theta(it), zeta(iz)]
                call vmec%evaluate_cyl(u_vmec, x_vmec_cyl)
                dr = abs(r_chart - x_vmec_cyl(1))
                dz = abs(z_chart - x_vmec_cyl(3))
                max_dr = max(max_dr, dr)
                max_dz = max(max_dz, dz)
            end do
        end do

        if (n_used < 1) then
            print *, "  FAIL: no usable zeta slices to compare"
            nerrors = nerrors + 1
            return
        end if

        if (max_dr > tol_r .or. max_dz > tol_z) then
            print *, "  FAIL: chartmap boundary differs from VMEC"
            print *, "    max|dR|=", max_dr, " tol=", tol_r
            print *, "    max|dZ|=", max_dz, " tol=", tol_z
            nerrors = nerrors + 1
            return
        end if

        print *, "  PASS: chartmap boundary matches VMEC max|dR|=", max_dr, &
            " max|dZ|=", max_dz
    end subroutine chartmap_boundary_matches_vmec_check

    subroutine chartmap_boundary_offset_check(chartmap_base, chartmap_offset, offset_cm, tol_cm, nerrors)
        character(len=*), intent(in) :: chartmap_base
        character(len=*), intent(in) :: chartmap_offset
        real(dp), intent(in) :: offset_cm
        real(dp), intent(in) :: tol_cm
        integer, intent(inout) :: nerrors

        integer :: ncid0, ncid1
        integer :: nrho0, ntheta0, nzeta0
        integer :: nrho1, ntheta1, nzeta1
        real(dp), allocatable :: x0(:, :, :), y0(:, :, :), z0(:, :, :)
        real(dp), allocatable :: x1(:, :, :), y1(:, :, :), z1(:, :, :)
        real(dp) :: r0, zc0, r1, zc1
        real(dp) :: dist_cm
        real(dp) :: max_abs_err
        integer :: ir, it, iz
        integer :: iz_list(5)
        integer :: nz, n_used
        integer :: k

        max_abs_err = 0.0_dp

        call nc_open(chartmap_base, ncid0)
        call nc_open(chartmap_offset, ncid1)

        call nc_inq_dim(ncid0, "rho", nrho0)
        call nc_inq_dim(ncid0, "theta", ntheta0)
        call nc_inq_dim(ncid0, "zeta", nzeta0)

        call nc_inq_dim(ncid1, "rho", nrho1)
        call nc_inq_dim(ncid1, "theta", ntheta1)
        call nc_inq_dim(ncid1, "zeta", nzeta1)

        if (nrho0 /= nrho1 .or. ntheta0 /= ntheta1 .or. nzeta0 /= nzeta1) then
            print *, "  FAIL: chartmap offset check requires matching grids"
            nerrors = nerrors + 1
            call nc_close(ncid0)
            call nc_close(ncid1)
            return
        end if

        allocate (x0(nrho0, ntheta0, nzeta0))
        allocate (y0(nrho0, ntheta0, nzeta0))
        allocate (z0(nrho0, ntheta0, nzeta0))
        allocate (x1(nrho1, ntheta1, nzeta1))
        allocate (y1(nrho1, ntheta1, nzeta1))
        allocate (z1(nrho1, ntheta1, nzeta1))

        call nc_get(ncid0, "x", x0)
        call nc_get(ncid0, "y", y0)
        call nc_get(ncid0, "z", z0)
        call nc_get(ncid1, "x", x1)
        call nc_get(ncid1, "y", y1)
        call nc_get(ncid1, "z", z1)

        call nc_close(ncid0)
        call nc_close(ncid1)

        nz = nzeta0
        iz_list = [1, 2, max(1, nz/2), max(1, nz - 1), nz]
        n_used = 0

        ir = nrho0
        do k = 1, size(iz_list)
            iz = iz_list(k)
            if (iz < 1 .or. iz > nz) cycle
            if (k > 1) then
                if (any(iz_list(1:k-1) == iz)) cycle
            end if
            n_used = n_used + 1

            do it = 1, ntheta0
                r0 = sqrt(x0(ir, it, iz)**2 + y0(ir, it, iz)**2)
                zc0 = z0(ir, it, iz)
                r1 = sqrt(x1(ir, it, iz)**2 + y1(ir, it, iz)**2)
                zc1 = z1(ir, it, iz)
                dist_cm = sqrt((r1 - r0)**2 + (zc1 - zc0)**2)
                max_abs_err = max(max_abs_err, abs(dist_cm - offset_cm))
            end do
        end do

        if (n_used < 1) then
            print *, "  FAIL: no usable zeta slices for offset check"
            nerrors = nerrors + 1
            return
        end if

        if (max_abs_err > tol_cm) then
            print *, "  FAIL: chartmap boundary offset mismatch"
            print *, "    expected offset [cm]=", offset_cm, " tol [cm]=", tol_cm
            print *, "    max|d-offset| [cm]=", max_abs_err
            nerrors = nerrors + 1
            return
        end if

        print *, "  PASS: chartmap boundary offset max|d-offset| [cm]=", max_abs_err
    end subroutine chartmap_boundary_offset_check

end module chartmap_test_utils
