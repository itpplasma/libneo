program test_chartmap_coordinates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
                                  make_chartmap_coordinate_system, &
                                  chartmap_coordinate_system_t
    use math_constants, only: TWOPI
    use nctools_module, only: nc_close, nc_get, nc_inq_dim, nc_open
    use netcdf
    implicit none

    integer :: nerrors
    logical :: all_passed
    class(coordinate_system_t), allocatable :: cs
    class(coordinate_system_t), allocatable :: cs_unknown
    character(len=*), parameter :: volume_file = "chartmap.nc"
    character(len=*), parameter :: unknown_file = "chartmap_unknown_convention.nc"

    nerrors = 0
    all_passed = .true.

    print *, "Testing chartmap coordinate system..."

    call make_chartmap_coordinate_system(cs, volume_file)

    if (.not. allocated(cs)) then
        print *, "  FAIL: make_chartmap_coordinate_system did not allocate cs"
        nerrors = nerrors + 1
    end if

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        call run_roundtrip_u_check(ccs, nerrors)
        call run_roundtrip_check(ccs, nerrors)
        call run_boundary_check(ccs, nerrors)
        call run_metric_check(ccs, nerrors)
    class default
        print *, "  FAIL: coordinate system is not chartmap type"
        nerrors = nerrors + 1
    end select

    call write_chartmap_with_convention(volume_file, unknown_file, "unknown")
    call make_chartmap_coordinate_system(cs_unknown, unknown_file)
    select type (ccs2 => cs_unknown)
    type is (chartmap_coordinate_system_t)
        call run_unknown_convention_check(ccs2, nerrors)
    class default
        print *, "  FAIL: unknown convention coordinate system is not chartmap type"
        nerrors = nerrors + 1
    end select

    if (nerrors > 0) then
        all_passed = .false.
        print *, "FAILED: ", nerrors, " error(s) detected in chartmap tests"
    else
        print *, "All chartmap coordinate tests passed!"
    end if

    if (.not. all_passed) error stop 1

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    subroutine write_chartmap_with_convention(srcfile, dstfile, convention)
        character(len=*), intent(in) :: srcfile
        character(len=*), intent(in) :: dstfile
        character(len=*), intent(in) :: convention

        integer :: ncid_src, ncid_dst
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta
        integer :: var_x, var_y, var_z
        integer :: var_num_field_periods
        integer :: nrho, ntheta, nzeta
        integer :: num_field_periods
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)

        call nc_open(srcfile, ncid_src)
        call nc_inq_dim(ncid_src, "rho", nrho)
        call nc_inq_dim(ncid_src, "theta", ntheta)
        call nc_inq_dim(ncid_src, "zeta", nzeta)
        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        allocate (x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                    ntheta, nzeta))
        call nc_get(ncid_src, "rho", rho)
        call nc_get(ncid_src, "theta", theta)
        call nc_get(ncid_src, "zeta", zeta)
        call nc_get(ncid_src, "x", x)
        call nc_get(ncid_src, "y", y)
        call nc_get(ncid_src, "z", z)
        num_field_periods = 1
        call nc_get(ncid_src, "num_field_periods", num_field_periods)
        call nc_close(ncid_src)

        call nc_check(nf90_create(trim(dstfile), NF90_NETCDF4, ncid_dst))
        call nc_check(nf90_put_att(ncid_dst, NF90_GLOBAL, "zeta_convention", &
                                   trim(convention)))
        call nc_check(nf90_def_dim(ncid_dst, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid_dst, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid_dst, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid_dst, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid_dst, "theta", NF90_DOUBLE, [dim_theta], &
                                   var_theta))
        call nc_check(nf90_def_var(ncid_dst, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid_dst, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid_dst, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid_dst, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid_dst, "num_field_periods", NF90_INT, &
                                   var_num_field_periods))
        call nc_check(nf90_put_att(ncid_dst, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid_dst, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid_dst, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid_dst))

        call nc_check(nf90_put_var(ncid_dst, var_rho, rho))
        call nc_check(nf90_put_var(ncid_dst, var_theta, theta))
        call nc_check(nf90_put_var(ncid_dst, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid_dst, var_x, x))
        call nc_check(nf90_put_var(ncid_dst, var_y, y))
        call nc_check(nf90_put_var(ncid_dst, var_z, z))
        call nc_check(nf90_put_var(ncid_dst, var_num_field_periods, num_field_periods))
        call nc_check(nf90_close(ncid_dst))
    end subroutine write_chartmap_with_convention

    subroutine run_roundtrip_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: x_ref(:, :, :), y_ref(:, :, :), z_ref(:, :, :)
        real(dp) :: rho_val, theta_val
        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        real(dp) :: diff_x
        integer :: ierr, ncid
        integer :: i_rho, i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol_x = 1.0e-10_dp

        call nc_open(volume_file, ncid)
        allocate (x_ref(nrho, ntheta, nzeta))
        allocate (y_ref(nrho, ntheta, nzeta))
        allocate (z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "x", x_ref)
        call nc_get(ncid, "y", y_ref)
        call nc_get(ncid, "z", z_ref)
        call nc_close(ncid)

        i_rho = 12
        i_theta = 17
        rho_val = real(i_rho - 1, dp)/real(nrho - 1, dp)
        theta_val = TWOPI*real(i_theta - 1, dp)/real(ntheta, dp)

        u = [rho_val, theta_val, 0.0_dp]
        call ccs%evaluate_point(u, x)

        if (abs(x(1) - x_ref(i_rho, i_theta, 1)) > tol_x .or. &
            abs(x(2) - y_ref(i_rho, i_theta, 1)) > tol_x .or. &
            abs(x(3) - z_ref(i_rho, i_theta, 1)) > tol_x) then
            print *, "  FAIL: forward map does not reproduce reference X,Y,Z"
            nerrors = nerrors + 1
        end if

        xcyl(1) = sqrt(x(1)**2 + x(2)**2)
        xcyl(2) = atan2(x(2), x(1))
        xcyl(3) = x(3)

        call ccs%from_cyl(xcyl, u_back, ierr)

        if (ierr /= 0) then
            print *, "  FAIL: inverse mapping reported error code ", ierr
            nerrors = nerrors + 1
        else
            call ccs%evaluate_point(u_back, x_round)
            diff_x = maxval(abs(x_round - x))

            if (diff_x > tol_x) then
                print *, "  FAIL: roundtrip mismatch in Cartesian coordinates, "// &
                    "max |x_round - x| = ", diff_x
                nerrors = nerrors + 1
            else
                print *, "  PASS: forward/inverse roundtrip against chartmap volume"
            end if
        end if

    end subroutine run_roundtrip_check

    subroutine run_roundtrip_u_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), x(3), u_back(3), xcyl(3), x_round(3)
        real(dp) :: dtheta, dzeta
        integer :: ierr
        integer :: ir, it, iz
        real(dp), parameter :: rho_vals(3) = [0.05_dp, 0.55_dp, 0.95_dp]
        real(dp), parameter :: theta_vals(4) = &
                               [0.0_dp, TWOPI*0.5_dp, TWOPI*0.25_dp, TWOPI*0.75_dp]
        real(dp), parameter :: zeta_vals(4) = &
                               [TWOPI*0.5_dp, 0.0_dp, TWOPI*0.25_dp, TWOPI*0.75_dp]
        real(dp), parameter :: tol_u = 1.0e-8_dp
        real(dp), parameter :: tol_x = 2.0e-8_dp
        integer :: nerrors_start

        nerrors_start = nerrors

        do iz = 1, 4
            do it = 1, 4
                do ir = 0, 2
                    u(1) = rho_vals(ir + 1)
                    u(2) = theta_vals(it)
                    u(3) = zeta_vals(iz)

                    call ccs%evaluate_point(u, x)

                    xcyl(1) = sqrt(x(1)**2 + x(2)**2)
                    xcyl(2) = atan2(x(2), x(1))
                    xcyl(3) = x(3)

                    call ccs%from_cyl(xcyl, u_back, ierr)
                    if (ierr /= 0) then
                        print *, "  FAIL: from_cyl ierr=", ierr, " for u=", u
                        nerrors = nerrors + 1
                        cycle
                    end if

                    dtheta = abs(modulo(u_back(2) - u(2) + 0.5_dp*TWOPI, TWOPI) - &
                                 0.5_dp*TWOPI)
                    dzeta = abs(modulo(u_back(3) - u(3) + 0.5_dp*TWOPI, TWOPI) - &
                                0.5_dp*TWOPI)

                    if (abs(u_back(1) - u(1)) > tol_u .or. dtheta > tol_u .or. &
                        dzeta > tol_u) then
                        print *, "  FAIL: u roundtrip mismatch u_back-u=", u_back - u
                        nerrors = nerrors + 1
                        cycle
                    end if

                    call ccs%evaluate_point(u_back, x_round)
                    if (maxval(abs(x_round - x)) > tol_x) then
                        print *, "  FAIL: x(u_back) mismatch max|dx|=", &
                            maxval(abs(x_round - x))
                        nerrors = nerrors + 1
                        cycle
                    end if
                end do
            end do
        end do

        if (nerrors == nerrors_start) then
            print *, &
                "  PASS: u->x->u_back roundtrip across zeta slices and near-axis points"
        end if
    end subroutine run_roundtrip_u_check

    subroutine run_unknown_convention_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp) :: u(3), x(3)
        real(dp) :: u_cyl(3), u_cart(3)
        real(dp) :: xcyl(3)
        real(dp) :: x_cyl(3), x_cart(3)
        real(dp) :: dtheta, dzeta
        integer :: ierr1, ierr2
        integer :: ir, it, iz
        integer :: nerrors_start
        real(dp), parameter :: rho_vals(2) = [0.05_dp, 0.75_dp]
        real(dp), parameter :: theta_vals(2) = [0.0_dp, TWOPI/2.0_dp]
        real(dp), parameter :: zeta_vals(2) = [0.0_dp, TWOPI/2.0_dp]
        real(dp), parameter :: tol_x = 2.0e-8_dp
        real(dp), parameter :: tol_u = 2.0e-8_dp

        nerrors_start = nerrors

        do iz = 1, size(zeta_vals)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u = [rho_vals(ir), theta_vals(it), zeta_vals(iz)]
                    call ccs%evaluate_point(u, x)

                    xcyl(1) = sqrt(x(1)**2 + x(2)**2)
                    xcyl(2) = atan2(x(2), x(1))
                    xcyl(3) = x(3)

                    call ccs%from_cyl(xcyl, u_cyl, ierr1)
                    if (ierr1 /= 0) then
                        print *, "  FAIL: unknown convention from_cyl ierr=", ierr1
                        nerrors = nerrors + 1
                        cycle
                    end if

                    call ccs%from_cart(x, u_cart, ierr2)
                    if (ierr2 /= 0) then
                        print *, "  FAIL: from_cart ierr=", ierr2
                        nerrors = nerrors + 1
                        cycle
                    end if

                    dtheta = modulo(u_cyl(2) - u(2) + TWOPI/2.0_dp, TWOPI) - &
                             TWOPI/2.0_dp
                    dzeta = modulo(u_cyl(3) - u(3) + TWOPI/2.0_dp, TWOPI) - TWOPI/2.0_dp
                    if (abs(u_cyl(1) - u(1)) > tol_u .or. abs(dtheta) > tol_u .or. &
                        abs(dzeta) > tol_u) then
                        print *, "  FAIL: unknown convention u mismatch"
                        nerrors = nerrors + 1
                        cycle
                    end if

                    call ccs%evaluate_point(u_cyl, x_cyl)
                    call ccs%evaluate_point(u_cart, x_cart)
                    if (maxval(abs(x_cyl - x)) > tol_x .or. maxval(abs(x_cart - x)) > &
                        tol_x) then
                        print *, "  FAIL: unknown convention x roundtrip mismatch"
                        nerrors = nerrors + 1
                        cycle
                    end if
                end do
            end do
        end do

        if (nerrors == nerrors_start) then
            print *, "  PASS: unknown zeta_convention uses cartesian inversion path"
        end if
    end subroutine run_unknown_convention_check

    subroutine run_boundary_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp), allocatable :: x_ref(:, :, :), y_ref(:, :, :), z_ref(:, :, :)
        real(dp) :: u(3), x(3)
        integer :: ncid
        integer :: i_theta
        integer, parameter :: nrho = 63, ntheta = 64, nzeta = 65
        real(dp), parameter :: tol = 1.0e-8_dp
        integer :: nerrors_local

        nerrors_local = 0

        call nc_open(volume_file, ncid)
        allocate (x_ref(nrho, ntheta, nzeta))
        allocate (y_ref(nrho, ntheta, nzeta))
        allocate (z_ref(nrho, ntheta, nzeta))
        call nc_get(ncid, "x", x_ref)
        call nc_get(ncid, "y", y_ref)
        call nc_get(ncid, "z", z_ref)
        call nc_close(ncid)

        do i_theta = 1, 4
            u(1) = 1.0_dp
            u(2) = TWOPI*real(i_theta - 1, dp)/real(ntheta, dp)
            u(3) = 0.0_dp

            call ccs%evaluate_point(u, x)

            if (abs(x(1) - x_ref(nrho, i_theta, 1)) > tol .or. &
                abs(x(2) - y_ref(nrho, i_theta, 1)) > tol .or. &
                abs(x(3) - z_ref(nrho, i_theta, 1)) > tol) then
                print *, "  FAIL: boundary point mismatch at theta index=", &
                    i_theta
                nerrors_local = nerrors_local + 1
            end if
        end do

        nerrors = nerrors + nerrors_local
        if (nerrors_local == 0) then
            print *, "  PASS: boundary recovery matches chartmap volume"
        end if
    end subroutine run_boundary_check

    subroutine run_metric_check(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors
        real(dp) :: u(3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: identity(3, 3), prod(3, 3)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: i, j, k
        logical :: sym_ok, id_ok

        u = [0.2_dp, 0.7_dp, 0.9_dp]
        call ccs%metric_tensor(u, g, ginv, sqrtg)

        if (sqrtg <= 0.0_dp) then
            print *, "  FAIL: Jacobian determinant not positive"
            nerrors = nerrors + 1
        else
            print *, "  PASS: Jacobian determinant positive"
        end if

        sym_ok = .true.
        do i = 1, 3
            do j = i + 1, 3
                if (abs(g(i, j) - g(j, i)) > tol) sym_ok = .false.
                if (abs(ginv(i, j) - ginv(j, i)) > tol) sym_ok = .false.
            end do
        end do
        if (.not. sym_ok) then
            print *, "  FAIL: metric tensor not symmetric"
            nerrors = nerrors + 1
        else
            print *, "  PASS: metric tensor symmetry"
        end if

        identity = 0.0_dp
        do i = 1, 3
            identity(i, i) = 1.0_dp
        end do

        prod = 0.0_dp
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    prod(i, j) = prod(i, j) + g(i, k)*ginv(k, j)
                end do
            end do
        end do

        id_ok = .true.
        do i = 1, 3
            do j = 1, 3
                if (abs(prod(i, j) - identity(i, j)) > tol) id_ok = .false.
            end do
        end do
        if (.not. id_ok) then
            print *, "  FAIL: g * ginv /= identity"
            nerrors = nerrors + 1
        else
            print *, "  PASS: g * ginv = identity"
        end if
    end subroutine run_metric_check

end program test_chartmap_coordinates
